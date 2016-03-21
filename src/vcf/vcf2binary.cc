/**
 * Copyright (c) 2016 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifdef HTSDIR

#include "vcf2binary.h"
#include "htslib/bgzf.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VCF2BinaryException(#X);

//Binary gets priority
#if defined(PRODUCE_BINARY_CELLS) and defined(PRODUCE_CSV_CELLS)
#undef PRODUCE_CSV_CELLS
#endif

//If none defined, produce binary cells
#if !defined(PRODUCE_BINARY_CELLS) and !defined(PRODUCE_CSV_CELLS)
#define PRODUCE_BINARY_CELLS 1
#endif

//VCFReader functions
VCFReader::VCFReader()
{
  m_indexed_reader = 0;
  m_fptr = 0;
  m_hdr = 0;
  m_line = bcf_init();
  m_is_line_valid = false;
  m_buffer.l = 0;
  m_buffer.m = 4096;    //4KB
  m_buffer.s = (char*)malloc(m_buffer.m*sizeof(char));
}

VCFReader::~VCFReader()
{
  if(m_indexed_reader)
    bcf_sr_destroy(m_indexed_reader);
  if(m_fptr)
    bcf_close(m_fptr);
  if(m_hdr)
    bcf_hdr_destroy(m_hdr);
  if(m_line)
    bcf_destroy(m_line);
  free(m_buffer.s);
}

void VCFReader::initialize(const char* filename, const char* regions,
    const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, bool open_file)
{
  assert(m_indexed_reader == 0);
  m_indexed_reader = bcf_sr_init();
  bcf_sr_set_regions(m_indexed_reader, regions, 0);   //random value
  m_filename = std::move(std::string(filename));
  if(open_file)
  {
    add_reader();
    m_hdr = bcf_hdr_dup(bcf_sr_get_header(m_indexed_reader, 0));
  }
  else
  {
    m_fptr = bcf_open(m_filename.c_str(), "r");
    VERIFY_OR_THROW(m_fptr && ("Cannot open VCF/BCF file "+m_filename).c_str());
    m_hdr = bcf_hdr_read(m_fptr);
    bcf_close(m_fptr);
    m_fptr = 0;
  }
  //Add lines in the header for fields that are missing in the VCF, but requested in the input config JSON file
  //Only look at INFO and FORMAT fields
  for(auto field_type_idx=BCF_HL_INFO;field_type_idx<=BCF_HL_FMT;++field_type_idx)
  {
    assert(static_cast<size_t>(field_type_idx) < vcf_field_names.size());
    for(auto j=0u;j<vcf_field_names[field_type_idx].size();++j)
    {
      const auto& field_name = vcf_field_names[field_type_idx][j];
      auto field_idx = bcf_hdr_id2int(m_hdr, BCF_DT_ID, field_name.c_str());
      //bcf_hdr_idinfo_exists handles negative field idx
      auto idinfo_exists = bcf_hdr_idinfo_exists(m_hdr, field_type_idx, field_idx);
      //Field not found
      if(idinfo_exists == 0)
      {
        std::string header_line = "##";
        header_line += (field_type_idx == BCF_HL_INFO ? "INFO" : "FORMAT");
        header_line += "=<ID="+field_name+",Type=";
        //GT is weird
        if(field_type_idx == BCF_HL_FMT && field_name == "GT")
          header_line += "String,Number=1,Description=\"Genotype\"";
        else
        {
          assert(id_mapper->get_field_info(field_name));
          auto field_info = *(id_mapper->get_field_info(field_name));
          switch(field_info.m_bcf_ht_type)
          {
            case BCF_HT_INT:
              header_line += "Integer";
              break;
            case BCF_HT_REAL:
              header_line += "Float";
              break;
            case BCF_HT_CHAR:
            case BCF_HT_STR:
              header_line += "String";
              break;
            default:
              throw VCF2BinaryException("Field type "+std::to_string(field_info.m_bcf_ht_type)+" not handled");
              break;
          }
          header_line += ",Number=";
          switch(field_info.m_length_descriptor)
          {
            case BCF_VL_FIXED:
              header_line += std::to_string(field_info.m_num_elements);
              break;
            case BCF_VL_VAR:
              header_line += ".";
              break;
            case BCF_VL_A:
              header_line += "A";
              break;
            case BCF_VL_R:
              header_line += "R";
              break;
            case BCF_VL_G:
              header_line += "G";
              break;
            default:
              throw VCF2BinaryException("Unhandled field length descriptor "+std::to_string(field_info.m_length_descriptor));
              break;
          }
          header_line += ",Description=\""+field_name+"\"";
        }
        header_line += ">";
        int line_length = 0;
        auto hrec = bcf_hdr_parse_line(m_hdr, header_line.c_str(), &line_length);
        bcf_hdr_add_hrec(m_hdr, hrec);
        bcf_hdr_sync(m_hdr);
      }
    }
  }
}

void VCFReader::add_reader()
{
  assert(m_indexed_reader->nreaders == 0);      //no existing files are open
  assert(m_fptr == 0);  //normal file handle should be NULL
  bcf_sr_add_reader(m_indexed_reader, m_filename.c_str());
}

void VCFReader::remove_reader()
{
  if(m_fptr)    //file handle moved to m_fptr after discarding index
  {
    assert(m_indexed_reader->nreaders == 0);
    bcf_close(m_fptr);
    m_fptr = 0;
  }
  else
    bcf_sr_remove_reader(m_indexed_reader, 0);
}

void VCFReader::seek_read_advance(const char* contig, const int pos, bool discard_index)
{
  //Close file handle if open
  if(m_fptr)
  {
    bcf_close(m_fptr);
    m_fptr = 0;
  }
  if(m_indexed_reader->nreaders == 0)        //index not loaded
    bcf_sr_add_reader(m_indexed_reader, m_filename.c_str());
  assert(m_indexed_reader->nreaders == 1);
  bcf_sr_seek(m_indexed_reader, contig, pos);
  //Only read 1 record at a time
  if(discard_index)
    m_indexed_reader->readers[0].read_one_record_only = 1;
  read_and_advance();
  if(discard_index)
  {
    std::swap<htsFile*>(m_fptr, m_indexed_reader->readers[0].file);
    assert(m_indexed_reader->readers[0].file == 0);
    bcf_sr_remove_reader(m_indexed_reader, 0);
  }
}

void VCFReader::read_and_advance()
{
  if(m_fptr)    //normal file handle - no index
  {
    //Handle VCFs and BCFs differently since the indexed reader handles file pointers differently
    if(m_fptr->format.format == htsExactFormat::vcf)
    {
      //Since m_fptr is obtained from an indexed reader, use bgzf_getline function
      auto status = bgzf_getline(hts_get_bgzfp(m_fptr), '\n', &m_buffer);
      m_is_line_valid = (status <= 0) ? false : true;
      if(m_is_line_valid)
        vcf_parse(&m_buffer, m_hdr, m_line);
    }
    else        //BCF
    {
      m_line->errcode = 0;
      //simple bcf_read
      auto status = bcf_read(m_fptr, m_hdr, m_line);
      m_is_line_valid = (status < 0) ? false : true;
      assert(m_line->errcode == 0);
    }
  }
  else  //indexed reader
  {
    bcf_sr_next_line(m_indexed_reader);
    auto line = bcf_sr_get_line(m_indexed_reader, 0);
    if(line)
    {
      std::swap<bcf1_t*>(m_indexed_reader->readers[0].buffer[0], m_line);
      m_is_line_valid = true;
    }
    else
      m_is_line_valid = false;
  }
}

//VCFColumnPartition functions
void VCFColumnPartition::copy_simple_members(const VCFColumnPartition& other)
{
  m_column_interval_begin = other.m_column_interval_begin;
  m_column_interval_end = other.m_column_interval_end;
  m_local_contig_idx = other.m_local_contig_idx;
  m_contig_position = other.m_contig_position;
  m_contig_tiledb_column_offset = other.m_contig_tiledb_column_offset;
  m_vcf_get_buffer_size = other.m_vcf_get_buffer_size;
}

//Move constructor
VCFColumnPartition::VCFColumnPartition(VCFColumnPartition&& other)
{
  copy_simple_members(other);
  m_begin_buffer_offset_for_local_callset = std::move(other.m_begin_buffer_offset_for_local_callset);
  m_last_full_line_end_buffer_offset_for_local_callset = 
    std::move(other.m_last_full_line_end_buffer_offset_for_local_callset);
  m_buffer_offset_for_local_callset = std::move(other.m_buffer_offset_for_local_callset);
  m_vcf_get_buffer = other.m_vcf_get_buffer;
  other.m_vcf_get_buffer = 0;
  other.m_vcf_get_buffer_size = 0;
  m_reader = other.m_reader;
  other.m_reader = 0;
}


VCF2Binary::VCF2Binary(const std::string& vcf_filename, const std::vector<std::vector<std::string>>& vcf_fields,
    unsigned file_idx, VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
    size_t max_size_per_callset,
    bool treat_deletions_as_intervals,
    bool parallel_partitions, bool noupdates, bool close_file, bool discard_index)
{
  m_reader = 0;
  m_histogram = 0;
  clear();
  m_vcf_filename = vcf_filename;
  m_vcf_fields = &vcf_fields;
  m_file_idx = file_idx;
  m_vid_mapper = &vid_mapper;
  m_max_size_per_callset = max_size_per_callset;
  m_treat_deletions_as_intervals = treat_deletions_as_intervals;
  m_parallel_partitions = parallel_partitions;
  m_noupdates = noupdates;
  m_discard_index = discard_index;
  m_close_file = close_file || discard_index;   //close file if index has to be discarded
  initialize(partition_bounds);
}

void VCF2Binary::copy_simple_members(const VCF2Binary& other)
{
  m_parallel_partitions = other.m_parallel_partitions;
  m_noupdates = other.m_noupdates;
  m_close_file = other.m_close_file;
  m_treat_deletions_as_intervals = other.m_treat_deletions_as_intervals;
  m_vid_mapper = other.m_vid_mapper;
  m_file_idx = other.m_file_idx;
  m_max_size_per_callset = other.m_max_size_per_callset;
  m_vcf_fields = other.m_vcf_fields;
  m_discard_index = other.m_discard_index;
}

//Move constructor
VCF2Binary::VCF2Binary(VCF2Binary&& other)
{
  copy_simple_members(other);
  m_vcf_filename = std::move(other.m_vcf_filename);
  m_regions = std::move(other.m_regions);
  m_local_callset_idx_to_tiledb_row_idx = std::move(other.m_local_callset_idx_to_tiledb_row_idx);
  m_enabled_local_callset_idx_vec = std::move(other.m_enabled_local_callset_idx_vec);
  m_local_contig_idx_to_global_contig_idx = std::move(other.m_local_contig_idx_to_global_contig_idx);
  m_local_field_idx_to_global_field_idx = std::move(other.m_local_field_idx_to_global_field_idx);
  m_partitions = std::move(other.m_partitions);
  m_reader = other.m_reader;
  other.m_reader = 0;
  m_histogram = other.m_histogram;
  other.m_histogram = 0;
}

VCF2Binary::~VCF2Binary()
{
  if(!m_parallel_partitions)
    for(auto& x : m_partitions)
      x.m_reader = 0;
  clear();
  if(m_reader)
    delete m_reader;
  m_reader = 0;
  if(m_histogram)
    delete m_histogram;
  m_histogram = 0;
}

void VCF2Binary::clear()
{
  m_vid_mapper = 0;
  m_vcf_filename.clear();
  m_regions.clear();
  m_local_callset_idx_to_tiledb_row_idx.clear();
  m_enabled_local_callset_idx_vec.clear();
  m_local_contig_idx_to_global_contig_idx.clear();
  m_local_field_idx_to_global_field_idx.clear();
  m_partitions.clear();
}

void VCF2Binary::initialize(const std::vector<ColumnRange>& partition_bounds)
{
  assert(m_vid_mapper);
  //Setup local-global mappings
  auto* fptr = bcf_open(m_vcf_filename.c_str(), "r");
  auto* hdr = bcf_hdr_read(fptr);
  //Callset mapping
  m_local_callset_idx_to_tiledb_row_idx = std::move(std::vector<int64_t>(bcf_hdr_nsamples(hdr), -1ll));
  m_vid_mapper->get_local_tiledb_row_idx_vec(m_vcf_filename, m_local_callset_idx_to_tiledb_row_idx);
  for(auto i=0ull;i<m_local_callset_idx_to_tiledb_row_idx.size();++i)
    if(m_local_callset_idx_to_tiledb_row_idx[i] >= 0)
      m_enabled_local_callset_idx_vec.push_back(i);
  //Contig mapping
  m_local_contig_idx_to_global_contig_idx = std::move(std::vector<int>(hdr->n[BCF_DT_CTG], -1ll));
  for(auto i=0;i<hdr->n[BCF_DT_CTG];++i)
    m_vid_mapper->get_global_contig_idx(bcf_hdr_id2name(hdr, i), m_local_contig_idx_to_global_contig_idx[i]);
  //Build regions list - comma separated contigs - in increasing order of column values
  m_regions = "";
  int64_t column_value = -1ll;
  std::string contig_name;
  bool first_valid_contig = true;
  while(m_vid_mapper->get_next_contig_location(column_value, contig_name, column_value))
  {
    auto local_contig_idx = bcf_hdr_name2id(hdr, contig_name.c_str());
    if(local_contig_idx < 0)
      continue;
    if(!first_valid_contig)
      m_regions += ",";
    m_regions += contig_name;
    first_valid_contig = false;
  }
  //Field mapping
  m_local_field_idx_to_global_field_idx = std::move(std::vector<int>(hdr->n[BCF_DT_ID], -1));
  for(auto i=0;i<hdr->n[BCF_DT_ID];++i)
    m_vid_mapper->get_global_field_idx(bcf_hdr_int2id(hdr, BCF_DT_ID, i), m_local_field_idx_to_global_field_idx[i]);
  bcf_hdr_destroy(hdr);
  bcf_close(fptr);
  //Initialize reader, if needed
  if(!m_parallel_partitions)
  {
    m_reader = new VCFReader();
    m_reader->initialize(m_vcf_filename.c_str(), m_regions.c_str(), *m_vcf_fields, m_vid_mapper, !m_close_file);
  }
  //Initialize partition info
  m_partitions.resize(partition_bounds.size());
  for(auto i=0u;i<partition_bounds.size();++i)
    initialize_partition(i, partition_bounds);
}

void VCF2Binary::initialize_partition(unsigned idx, const std::vector<ColumnRange>& partition_bounds)
{
  auto& column_interval_info = m_partitions[idx];
  column_interval_info.m_column_interval_begin = partition_bounds[idx].first;
  column_interval_info.m_column_interval_end = partition_bounds[idx].second;
  //If parallel partitions, each interval gets its own reader
  if(m_parallel_partitions)
  {
    column_interval_info.m_reader = new VCFReader();
    column_interval_info.m_reader->initialize(m_vcf_filename.c_str(), m_regions.c_str(), *m_vcf_fields, m_vid_mapper, !m_close_file);
  }
  else
    column_interval_info.m_reader = m_reader;
  //buffer offsets for each callset in this partition
  column_interval_info.m_buffer_offset_for_local_callset.resize(m_enabled_local_callset_idx_vec.size());
  column_interval_info.m_begin_buffer_offset_for_local_callset.resize(m_enabled_local_callset_idx_vec.size());
  column_interval_info.m_last_full_line_end_buffer_offset_for_local_callset.resize(m_enabled_local_callset_idx_vec.size());
  //Indicates that nothing has been read for this interval
  column_interval_info.m_local_contig_idx = -1;
  column_interval_info.m_contig_position = -1;
}

void VCF2Binary::set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const
{
  for(auto local_callset_idx : m_enabled_local_callset_idx_vec)
  {
    assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size());
    auto row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
    assert(row_idx >= 0);
    assert(static_cast<size_t>(row_idx) < tiledb_row_idx_to_order.size());
    tiledb_row_idx_to_order[row_idx] = order_value++;
  }
}

void VCF2Binary::list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const
{
  auto& partition_file_batch = partition_batch.get_partition_file_batch(m_file_idx);
  if(partition_file_batch.m_fetch && !partition_file_batch.m_completed)
  {
    for(auto local_callset_idx : m_enabled_local_callset_idx_vec)
    {
      assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size());
      auto row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
      assert(row_idx >= 0);
      row_idx_vec[row_idx_offset++] = row_idx;
    }
  }
}

void VCF2Binary::read_next_batch(std::vector<std::vector<uint8_t>*>& buffer_vec,
    std::vector<ColumnPartitionBatch>& partition_batches, bool close_file)
{
  if(m_parallel_partitions)
  {
#pragma omp parallel for
    for(auto partition_idx=0u;partition_idx<partition_batches.size();++partition_idx)
    {
      auto& curr_file_batch = partition_batches[partition_idx].get_partition_file_batch(m_file_idx);
      assert(static_cast<size_t>(curr_file_batch.get_buffer_idx()) < buffer_vec.size());
      read_next_batch(*(buffer_vec[curr_file_batch.get_buffer_idx()]), m_partitions[partition_idx], curr_file_batch, close_file);
    }
  }
  else
  {
    //Open file handles if needed
    if(m_close_file)
      m_reader->add_reader();
    for(auto partition_idx=0u;partition_idx<partition_batches.size();++partition_idx)
    {
      auto& curr_file_batch = partition_batches[partition_idx].get_partition_file_batch(m_file_idx);
      assert(static_cast<size_t>(curr_file_batch.get_buffer_idx()) < buffer_vec.size());
      read_next_batch(*(buffer_vec[curr_file_batch.get_buffer_idx()]), m_partitions[partition_idx], curr_file_batch, close_file);
    }
    //Close file handles if needed
    if(close_file)
      m_reader->remove_reader();
  }
  m_close_file = close_file;
}

void VCF2Binary::read_next_batch(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
    ColumnPartitionFileBatch& partition_file_batch, bool close_file)
{
  //Nothing to do
  if(!partition_file_batch.m_fetch || partition_file_batch.m_completed)
    return;
  //Open file handles if needed
  if(m_parallel_partitions && m_close_file)
    vcf_partition.m_reader->add_reader();
  //Setup buffer offsets first 
  for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
  {
    auto curr_offset = partition_file_batch.get_offset_for_local_callset_idx(i, m_max_size_per_callset);
    vcf_partition.m_begin_buffer_offset_for_local_callset[i] = curr_offset;
    vcf_partition.m_buffer_offset_for_local_callset[i] = curr_offset;
    vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[i] = curr_offset;
  }
  //If file is re-opened, seek to position from which to begin reading, but do not advance iterator from previous position
  //The second parameter is useful if the file handler is open, but the buffer was full in a previous call
  auto has_data = seek_and_fetch_position(vcf_partition, m_close_file, false);
  auto buffer_full = false;
  auto read_one_line_fully = false;
  while(has_data && !buffer_full)
  {
    auto* line = vcf_partition.m_reader->get_line();
    bcf_unpack(line, BCF_UN_ALL);
    for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
    {
      buffer_full = buffer_full || convert_VCF_to_binary_for_callset(buffer, vcf_partition, m_max_size_per_callset, i);
      if(buffer_full)
        break;
    }
    vcf_partition.m_contig_position = line->pos;      //keep track of position from which to seek next time
    if(!buffer_full)
    {
      has_data = seek_and_fetch_position(vcf_partition, false, true);  //no need to re-seek, use next_line() directly, advance file pointer
      //Store buffer offsets at the beginning of the line
      for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
        vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[i] = vcf_partition.m_buffer_offset_for_local_callset[i];
      read_one_line_fully = true;
    } 
    else
      VERIFY_OR_THROW(read_one_line_fully && "Buffer did not have space to hold a line fully - increase buffer size")
  }
  //put Tiledb NULL for row_idx as end-of-batch marker 
  for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
  {
#ifdef PRODUCE_BINARY_CELLS
    tiledb_buffer_print_null<int64_t>(buffer, vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[i], 
        vcf_partition.m_begin_buffer_offset_for_local_callset[i] + m_max_size_per_callset);
#endif
#ifdef PRODUCE_CSV_CELLS
    tiledb_buffer_print<char>(buffer, vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[i], 
        vcf_partition.m_begin_buffer_offset_for_local_callset[i] + m_max_size_per_callset, '\0', false);
#endif
  }
  if(!has_data)
    partition_file_batch.m_completed = true;
  //Close file handles if needed
  if(m_parallel_partitions && close_file)
    vcf_partition.m_reader->remove_reader();
  //Advance write idx for partition_file_batch
  partition_file_batch.advance_write_idx();
}

inline void VCF2Binary::update_local_contig_idx(VCFColumnPartition& vcf_partition, const bcf1_t* line)
{
  //Different contig, update offset value
  if(line->rid != vcf_partition.m_local_contig_idx)
  {
    auto global_contig_idx = m_local_contig_idx_to_global_contig_idx[line->rid];
    //FIXME: contigs in this file that are unknown globally
    VERIFY_OR_THROW(global_contig_idx >= 0);
    vcf_partition.m_contig_tiledb_column_offset = m_vid_mapper->get_contig_info(global_contig_idx).m_tiledb_column_offset;
    vcf_partition.m_local_contig_idx = line->rid;
  }
}

bool VCF2Binary::seek_and_fetch_position(VCFColumnPartition& vcf_partition, bool force_seek, bool advance_reader)
{
  auto hdr = vcf_partition.m_reader->get_header();
  //If valid contig, i.e., continuing from a valid previous position
  //Common case: m_local_contig_idx >=0 and !force_seek and advance_reader  
  if(vcf_partition.m_local_contig_idx >= 0)
  {
    if(force_seek)
      vcf_partition.m_reader->seek_read_advance(bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx), vcf_partition.m_contig_position,
          m_discard_index);
    else
      if(advance_reader)
        vcf_partition.m_reader->read_and_advance();
    auto* line = vcf_partition.m_reader->get_line();
    //Next line exists && (either index is stored in memory or continuing in the same contig) - common case
    if(line && (!m_discard_index || line->rid == vcf_partition.m_local_contig_idx))
    {
      update_local_contig_idx(vcf_partition, line);
      if(vcf_partition.m_contig_tiledb_column_offset+static_cast<int64_t>(line->pos) <= vcf_partition.m_column_interval_end)
        return true;
      else
        return false;   //past column end
    }
    else
      if(!m_discard_index)  //index is in memory - implies there are no more lines found in file
        return false;
      else   //Index NOT in memory - either no line found or switched to new contig, either way need to re-seek
        vcf_partition.m_local_contig_idx = -1;
  }
  else //un-initialized
  {
    std::string contig_name;
    int64_t contig_position = -1;
    auto exists_contig = m_vid_mapper->get_contig_location(vcf_partition.m_column_interval_begin, 
        contig_name, contig_position);
    if(exists_contig) //decrement below the current contig so that next_contig gets the correct contig subsequently
      vcf_partition.m_contig_tiledb_column_offset = vcf_partition.m_column_interval_begin - contig_position - 1;
    else
      vcf_partition.m_contig_tiledb_column_offset = vcf_partition.m_column_interval_begin;
  }
  //while valid contig not found
  //Exits when valid contig found or next_contig() scans through all known contigs
  while(vcf_partition.m_local_contig_idx < 0)
  {
    std::string contig_name;
    auto exists_contig = m_vid_mapper->get_next_contig_location(vcf_partition.m_contig_tiledb_column_offset,
        contig_name, vcf_partition.m_contig_tiledb_column_offset);
    if(!exists_contig)  //no more valid contigs found, exit
      return false;
    //past column end
    if(vcf_partition.m_contig_tiledb_column_offset > vcf_partition.m_column_interval_end)
      return false;
    vcf_partition.m_local_contig_idx = bcf_hdr_name2id(hdr, contig_name.c_str());
    if(vcf_partition.m_local_contig_idx >= 0)
    {
      //Contig pos - if interval begin is in the middle of a contig, then compute diff, else start of contig 
      auto contig_pos = vcf_partition.m_column_interval_begin > vcf_partition.m_contig_tiledb_column_offset
        ? vcf_partition.m_column_interval_begin - vcf_partition.m_contig_tiledb_column_offset : 0ll;
      //Seek to start position
      vcf_partition.m_reader->seek_read_advance(bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx), contig_pos, m_discard_index);
      auto* line = vcf_partition.m_reader->get_line();
      //no more valid lines found, exit
      if(line == 0)
        return false;
      else      //valid VCF record
      {
        update_local_contig_idx(vcf_partition, line);
        //and is within this column interval, return true
        if(vcf_partition.m_contig_tiledb_column_offset+static_cast<int64_t>(line->pos) <= vcf_partition.m_column_interval_end)
          return true;
        else
          return false;
      }
    }
  }
  return false;
}

#ifdef PRODUCE_BINARY_CELLS

template<class FieldType>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const FieldType val, bool print_sep)
{
  int64_t add_size = sizeof(FieldType);
  //Do not write anything past the limit
  if(buffer_offset + add_size > buffer_offset_limit)
    return true;
  FieldType* ptr = reinterpret_cast<FieldType*>(&(buffer[buffer_offset]));
  *ptr = val;
  buffer_offset += add_size;
  return false;
}

//specialization for char* type
template<>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const char* val, bool print_sep)
{
  int64_t add_size = strlen(val)*sizeof(char);
  //Do not write anything past the limit
  if(buffer_offset + add_size > buffer_offset_limit)
    return true;
  memcpy(&(buffer[buffer_offset]), val, add_size);
  buffer_offset += add_size;
  return false;
}

template<>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, char* val, bool print_sep)
{
  return tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
}

//specialization for std::string type
template<>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const std::string& val, bool print_sep)
{
  int64_t add_size = val.length()*sizeof(char);
  //Do not write anything past the limit
  if(buffer_offset + add_size > buffer_offset_limit)
    return true;
  memcpy(&(buffer[buffer_offset]), val.c_str(), add_size);
  buffer_offset += add_size;
  return false;
}

template<>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, std::string& val, bool print_sep)
{
  return tiledb_buffer_print<const std::string&>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
}

template<class FieldType>
inline bool VCF2Binary::tiledb_buffer_print_null(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit)
{
  return tiledb_buffer_print<FieldType>(buffer, buffer_offset, buffer_offset_limit, get_tiledb_null_value<FieldType>());
}

#endif //ifdef PRODUCE_BINARY_CELLS

#ifdef PRODUCE_CSV_CELLS
template<class FieldType>
inline bool VCF2Binary::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const FieldType val, bool print_sep)
{
  std::stringstream ss;
  if(print_sep)
    ss << "," << val;
  else
    ss << val;
  std::string x = ss.str();
  size_t add_size = x.length()*sizeof(char);
  //Do not write anything past the limit
  if(buffer_offset + add_size > buffer_offset_limit)
    return true;
  memcpy(&(buffer[buffer_offset]), x.c_str(), add_size);
  buffer_offset += add_size;
  return false;
}

template<class FieldType>
inline bool VCF2Binary::tiledb_buffer_print_null(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit)
{
  return tiledb_buffer_print<char>(buffer, buffer_offset, buffer_offset_limit, NULL_VALUE);
}
#endif

template<class FieldType>
bool VCF2Binary::convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition, 
    int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
    const std::string& field_name, unsigned field_type_idx)
{
  auto* hdr = vcf_partition.m_reader->get_header();
  auto* line = vcf_partition.m_reader->get_line();
  //FIXME: avoid strings
  auto is_GT_field = (field_type_idx == BCF_HL_FMT && field_name == "GT");
  assert(line);
  assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size()
      && local_callset_idx < bcf_hdr_nsamples(hdr));
  auto field_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str());
  //This should always pass as missing fields are added to the header during initialization
  //Check left in for safety
  VERIFY_OR_THROW(field_idx >= 0 && bcf_hdr_idinfo_exists(hdr, field_type_idx, field_idx));
  //FIXME: special length descriptors
  auto length_descriptor = is_GT_field ? BCF_VL_P : bcf_hdr_id2length(hdr, field_type_idx, field_idx);
  auto field_length = bcf_hdr_id2number(hdr, field_type_idx, field_idx);
  auto bcf_ht_type = is_GT_field ? BCF_HT_INT : bcf_hdr_id2type(hdr, field_type_idx, field_idx);
  //The weirdness of VCF - string fields are marked as fixed length fields of size 1 (*facepalm*)
  auto is_vcf_str_type = (bcf_ht_type == BCF_HT_CHAR || bcf_ht_type == BCF_HT_STR) && !is_GT_field;
  length_descriptor =  is_vcf_str_type ? BCF_VL_VAR : length_descriptor;
  int max_num_values = vcf_partition.m_vcf_get_buffer_size/sizeof(FieldType);
  auto num_values = (field_type_idx == BCF_HL_INFO) ?
    bcf_get_info_values(hdr, line, field_name.c_str(), reinterpret_cast<void**>(&(vcf_partition.m_vcf_get_buffer)), &max_num_values, bcf_ht_type)
    : bcf_get_format_values(hdr, line, field_name.c_str(), reinterpret_cast<void**>(&(vcf_partition.m_vcf_get_buffer)), &max_num_values, bcf_ht_type);
  if(static_cast<uint64_t>(max_num_values)*sizeof(FieldType) >  vcf_partition.m_vcf_get_buffer_size)
    vcf_partition.m_vcf_get_buffer_size = static_cast<uint64_t>(max_num_values)*sizeof(FieldType);
  auto buffer_full = false;
  if(num_values < 0) //Curr line does not have this field
  {
    //variable length field, print #elements = 0
    if(length_descriptor != BCF_VL_FIXED)
    {
#ifdef PRODUCE_CSV_CELLS
      if(!is_vcf_str_type)
#endif
      {
        buffer_full = tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, 0);
      }
#ifdef PRODUCE_CSV_CELLS
      else
      {
        buffer_full = tiledb_buffer_print_null<FieldType>(buffer, buffer_offset, buffer_offset_limit);
      }
#endif
    }
    else        //fixed length field - fill with NULL values
      for(auto i=0u;i<field_length;++i)
      {
        buffer_full = buffer_full || tiledb_buffer_print_null<FieldType>(buffer, buffer_offset, buffer_offset_limit);
        if(buffer_full) return true;
      }
  }
  else
  {
    //variable length field, print #elements  first
    if(length_descriptor != BCF_VL_FIXED)
    {
#ifdef PRODUCE_CSV_CELLS
      if(!is_vcf_str_type)
#endif
      {
        buffer_full = tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, num_values);
        if(buffer_full) return true;
      }
    }
    auto* ptr = reinterpret_cast<const FieldType*>(vcf_partition.m_vcf_get_buffer);
    //For format fields, the ptr should point to where data for the current callset begins
    if(field_type_idx == BCF_HL_FMT)
    {
      assert(num_values%bcf_hdr_nsamples(hdr) == 0);
      num_values = num_values/bcf_hdr_nsamples(hdr);
      ptr += (local_callset_idx*num_values);
    }
    auto print_sep = true;
    for(auto k=0;k<num_values;++k)
    {
      auto val = ptr[k];
      if(is_GT_field)
        val = bcf_gt_allele(static_cast<int>(val));
      buffer_full = buffer_full || tiledb_buffer_print<FieldType>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
      if(buffer_full) return true;
      print_sep  = !is_vcf_str_type;
    }
  }
  return buffer_full;
}

bool VCF2Binary::convert_VCF_to_binary_for_callset(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
    size_t size_per_callset, uint64_t enabled_callsets_idx)
{
  auto* hdr = vcf_partition.m_reader->get_header();
  auto* line = vcf_partition.m_reader->get_line();
  assert(line);
  assert(enabled_callsets_idx < m_enabled_local_callset_idx_vec.size());
  auto local_callset_idx = m_enabled_local_callset_idx_vec[enabled_callsets_idx];
  assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size()
      && local_callset_idx < bcf_hdr_nsamples(hdr));
  assert(vcf_partition.m_local_contig_idx >= 0 && vcf_partition.m_local_contig_idx == line->rid);
  assert(vcf_partition.m_contig_tiledb_column_offset >= 0);
  //Buffer offsets tracking
  const int64_t begin_buffer_offset = vcf_partition.m_begin_buffer_offset_for_local_callset[enabled_callsets_idx];
  const int64_t line_begin_buffer_offset = vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[enabled_callsets_idx];
  int64_t& buffer_offset = vcf_partition.m_buffer_offset_for_local_callset[enabled_callsets_idx];
  assert(line_begin_buffer_offset >= begin_buffer_offset && line_begin_buffer_offset <= static_cast<int64_t>(begin_buffer_offset + size_per_callset));
  assert(buffer_offset >= begin_buffer_offset && buffer_offset <= static_cast<int64_t>(begin_buffer_offset + size_per_callset));
  assert(buffer_offset >= line_begin_buffer_offset);
  const int64_t buffer_offset_limit = begin_buffer_offset + size_per_callset;
  bool buffer_full = false;
  //Row
  int64_t row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
  if(row_idx < 0) return false;
  buffer_full = buffer_full || tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, row_idx, false);
  if(buffer_full) return true;
  //Column
  int64_t column_idx = vcf_partition.m_contig_tiledb_column_offset + line->pos;
  buffer_full = buffer_full || tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, column_idx);
  if(buffer_full) return true;
#ifdef PRODUCE_BINARY_CELLS
  //For binary cells, must print size of cell here. Keep track of offset and update later
  auto cell_size_offset = buffer_offset;
  buffer_offset += sizeof(size_t);
#endif
  //END position
  int max_num_values = vcf_partition.m_vcf_get_buffer_size/sizeof(int);
  //FIXME: avoid strings
  auto num_values = bcf_get_info_int32(hdr, line, "END", &(vcf_partition.m_vcf_get_buffer), &max_num_values);
  assert(num_values == 1 || num_values == -3);
  auto end_column_idx = column_idx;
  if(num_values < 0)    //missing end value
  {
    //handle spanning deletions
    if(m_treat_deletions_as_intervals)
    {
      auto alleles = line->d.allele;
      auto ref_length = strlen(alleles[0]);
      for(auto j=1;j<line->n_allele;++j)
      {    
        if(bcf_get_variant_type(line, j) == VCF_INDEL && ref_length > strlen(alleles[j]))
        {    
          end_column_idx = column_idx + ref_length - 1;
          break;
        }    
      }    
    }
  }
  else  //valid END found
    end_column_idx = vcf_partition.m_contig_tiledb_column_offset + *(reinterpret_cast<int*>(vcf_partition.m_vcf_get_buffer)) - 1; //convert 1-based END to 0-based
  buffer_full = buffer_full || tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, end_column_idx);
  if(buffer_full) return true;
  //REF
#ifdef PRODUCE_BINARY_CELLS
  buffer_full = buffer_full || tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, strlen(line->d.allele[0]));
  if(buffer_full) return true;
#endif
  buffer_full = buffer_full || tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, line->d.allele[0]);
  if(buffer_full) return true;
  //ALT
#ifdef PRODUCE_BINARY_CELLS
  //Must write length of ALT string in buffer, keep track of offset and update later
  auto alt_length_offset = buffer_offset;
  buffer_offset += sizeof(int);
#endif
  std::string alt_allele_serialized = std::move("");
  for(auto i=1;i<line->n_allele;++i)
  {
    if(i > 1)
      alt_allele_serialized += TILEDB_ALT_ALLELE_SEPARATOR;
    alt_allele_serialized += (bcf_get_variant_type(line, i) == VCF_NON_REF) ? TILEDB_NON_REF_VARIANT_REPRESENTATION
      : line->d.allele[i];
  }
  buffer_full = buffer_full || tiledb_buffer_print<const std::string&>(buffer, buffer_offset, buffer_offset_limit, alt_allele_serialized);
  if(buffer_full) return true;
#ifdef PRODUCE_BINARY_CELLS
  //write length of alt alleles string
  buffer_full = buffer_full ||  tiledb_buffer_print<int>(buffer, alt_length_offset, buffer_offset_limit, alt_allele_serialized.length());
  if(buffer_full) return true;
#endif
  //QUAL
  buffer_full = buffer_full || ( is_bcf_missing_value<float>(line->qual)
      ? tiledb_buffer_print_null<float>(buffer, buffer_offset, buffer_offset_limit) 
      : tiledb_buffer_print<float>(buffer, buffer_offset, buffer_offset_limit, line->qual) );
  if(buffer_full) return true;
  //Filter
  buffer_full = buffer_full ||  tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, line->d.n_flt);
  if(buffer_full) return true;
  for(auto i=0;i<line->d.n_flt;++i)
  {
    assert(line->d.flt[i] < static_cast<int64_t>(m_local_field_idx_to_global_field_idx.size()));
    buffer_full = buffer_full ||  tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit,
        m_local_field_idx_to_global_field_idx[line->d.flt[i]]);
    if(buffer_full) return true;
  }
  //Get INFO and FORMAT fields
  for(auto field_type_idx=BCF_HL_INFO;field_type_idx<=BCF_HL_FMT;++field_type_idx)
  {
    assert(static_cast<size_t>(field_type_idx) < m_vcf_fields->size());
    for(auto j=0u;j<(*m_vcf_fields)[field_type_idx].size();++j)
    {
      const auto& field_name = (*m_vcf_fields)[field_type_idx][j];
      //FIXME: avoid strings
      if(field_type_idx == BCF_HL_INFO && field_name == "END")   //ignore END field
        continue;
      auto field_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str());
      //Should always pass - left in for safety
      VERIFY_OR_THROW(field_idx >= 0 && bcf_hdr_idinfo_exists(hdr, field_type_idx, field_idx));
      auto field_ht_type = bcf_hdr_id2type(hdr, field_type_idx, field_idx);
      //Because GT is encoded type string in VCF - total nonsense
      //FIXME: avoid strings
      field_ht_type = (field_type_idx == BCF_HL_FMT && field_name == "GT") ? BCF_HT_INT : field_ht_type;
      switch(field_ht_type)
      {
        case BCF_HT_INT:
          buffer_full = buffer_full || convert_field_to_tiledb<int>(buffer, vcf_partition, buffer_offset, buffer_offset_limit, local_callset_idx,
              field_name, field_type_idx);
          if(buffer_full) return true;
          break;
        case BCF_HT_REAL:
          buffer_full = buffer_full || convert_field_to_tiledb<float>(buffer, vcf_partition, buffer_offset, buffer_offset_limit, local_callset_idx,
              field_name, field_type_idx);
          if(buffer_full) return true;
          break;
        case BCF_HT_STR:
        case BCF_HT_CHAR:
          buffer_full = buffer_full || convert_field_to_tiledb<char>(buffer, vcf_partition, buffer_offset, buffer_offset_limit, local_callset_idx,
              field_name, field_type_idx);
          if(buffer_full) return true;
          break;
        default: //FIXME: handle other types
          throw VCF2BinaryException("Unhandled VCF data type "+bcf_hdr_id2type(hdr, BCF_DT_ID, field_idx));
          break;
      }
    }
  }
#ifdef PRODUCE_BINARY_CELLS
  //Update total size
  buffer_full = buffer_full ||  tiledb_buffer_print<size_t>(buffer, cell_size_offset, buffer_offset_limit, buffer_offset-line_begin_buffer_offset);
  if(buffer_full) return true;
#endif
#ifdef PRODUCE_CSV_CELLS
  //Add newline
  if(buffer_offset+sizeof(char) <= buffer_offset_limit)
  {
    buffer[buffer_offset] = '\n';
    buffer_offset += sizeof(char);
  }
  else
    buffer_full = true;
  if(buffer_full) return true;
#endif
  return buffer_full;
}

void VCF2Binary::create_histogram_for_vcf(uint64_t max_histogram_range, unsigned num_bins)
{
  if(m_histogram)
    delete m_histogram;
  m_histogram = new UniformHistogram(0, max_histogram_range, num_bins);
  //Open file handles if needed
  if(m_close_file)
    m_reader->add_reader();
  //Should really have 1 partition only
  for(auto& vcf_partition : m_partitions)
  {
    auto has_data = seek_and_fetch_position(vcf_partition, m_close_file, false);
    while(has_data)
    {
      auto line = m_reader->get_line();
      int64_t column_idx = vcf_partition.m_contig_tiledb_column_offset + line->pos;
      for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
        m_histogram->add_interval(column_idx, column_idx);
      has_data = seek_and_fetch_position(vcf_partition, false, true);
    }
  }
  if(m_close_file)
    m_reader->remove_reader();
}

#endif //ifdef HTSDIR
