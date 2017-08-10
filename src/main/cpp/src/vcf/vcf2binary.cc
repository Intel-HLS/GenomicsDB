/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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
#include "vcf_adapter.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VCF2BinaryException(#X);

void VCFReaderBase::initialize(const char* filename,
    const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, const bool open_file)
{
  assert(m_hdr);
  m_name = std::move(std::string(filename));
  //Add lines in the header for fields that are missing in the VCF, but requested in the input config JSON file
  for(auto field_type_idx=BCF_HL_FLT;field_type_idx<=BCF_HL_FMT;++field_type_idx)
  {
    assert(static_cast<size_t>(field_type_idx) < vcf_field_names.size());
    try
    {
      for(auto j=0u;j<vcf_field_names[field_type_idx].size();++j)
        VCFAdapter::add_field_to_hdr_if_missing(m_hdr, id_mapper, vcf_field_names[field_type_idx][j], field_type_idx);
    }
    catch(const VCFAdapterException& e)
    {
      std::cerr << "ERROR: conflicting field description in the vid JSON and the VCF header of file: "<<filename<<"\n";
      throw e;
    }
  }
}

//VCFBufferReader functions
void VCFBufferReader::initialize(const char* stream_name,
    const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, const bool open_file)
{
  //Buffer MUST contain header information when initialize is called
  VERIFY_OR_THROW(contains_unread_data());
  m_hdr = bcf_hdr_init("r");
  VERIFY_OR_THROW(m_hdr);
  size_t hdr_length = 0ull;
  auto new_offset = bcf_hdr_deserialize(m_hdr, &(BufferReaderBase::m_buffer[0]), 0u, BufferReaderBase::m_num_valid_bytes_in_buffer, m_is_bcf ? 1 : 0);
  //Header sample line parsed correctly - might have ignored other incorrect lines 
  if(new_offset == 0u)
    throw VCF2BinaryException(std::string("Could not parse ")+(m_is_bcf ? "BCF" : "VCF")+" header for stream "+stream_name);
  advance_offset_by(hdr_length);
  VCFReaderBase::initialize(stream_name, vcf_field_names, id_mapper, open_file);
}

void VCFBufferReader::read_and_advance()
{
  m_is_record_valid = false;
  if(BufferReaderBase::contains_unread_data())
  {
    auto new_offset = bcf_deserialize(m_line, &(BufferReaderBase::m_buffer[0]), BufferReaderBase::m_offset, BufferReaderBase::m_num_valid_bytes_in_buffer,
        m_is_bcf ? 1 : 0, m_hdr);
    //Parsed or made progress
    assert(new_offset > BufferReaderBase::m_offset);
    BufferReaderBase::m_offset = new_offset;
    m_is_record_valid = true;
  }
}

//VCFReader functions
VCFReader::VCFReader()
  : GenomicsDBImportReaderBase(true), FileReaderBase(), VCFReaderBase(true)
{
  m_indexed_reader = 0;
  m_fptr = 0;
  m_vcf_file_buffer.l = 0;
  m_vcf_file_buffer.m = 4096;    //4KB
  m_vcf_file_buffer.s = (char*)malloc(m_vcf_file_buffer.m*sizeof(char));
}

VCFReader::~VCFReader()
{
  if(m_indexed_reader)
    bcf_sr_destroy(m_indexed_reader);
  m_indexed_reader = 0;
  if(m_fptr)
    bcf_close(m_fptr);
  m_fptr = 0;
  if(m_vcf_file_buffer.s && m_vcf_file_buffer.m)
    free(m_vcf_file_buffer.s);
  m_vcf_file_buffer.s = 0;
  m_vcf_file_buffer.m = 0;
}

void VCFReader::initialize(const char* filename,
    const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, const bool open_file)
{
  //Build regions list - comma separated contigs - in increasing order of column values
  //Needed to fix traversal order of indexed reader
  //So parse the header first
  m_fptr = bcf_open(filename, "r");
  VERIFY_OR_THROW(m_fptr && (std::string("Cannot open VCF/BCF file ")+filename).c_str());
  m_hdr = bcf_hdr_read(m_fptr);
  bcf_close(m_fptr);
  m_fptr = 0;
  std::string regions = "";
  int64_t column_value = -1ll;
  std::string contig_name;
  bool first_valid_contig = true;
  while(id_mapper->get_next_contig_location(column_value, contig_name, column_value))
  {
    auto local_contig_idx = bcf_hdr_name2id(m_hdr, contig_name.c_str());
    if(local_contig_idx < 0)
      continue;
    if(!first_valid_contig)
      regions += ",";
    //enclose contig name with quotes to deal with weird contig names containing :,- etc - messes up synced reader
    regions += ('"' + contig_name + '"');
    first_valid_contig = false;
  }
  assert(m_indexed_reader == 0);
  m_indexed_reader = bcf_sr_init();
  bcf_sr_set_regions(m_indexed_reader, regions.c_str(), 0);
  VCFReaderBase::initialize(filename, vcf_field_names, id_mapper, open_file);
  if(open_file)
    add_reader();
}

void VCFReader::add_reader()
{
  assert(m_indexed_reader->nreaders == 0);      //no existing files are open
  assert(m_fptr == 0);  //normal file handle should be NULL
  if(bcf_sr_add_reader(m_indexed_reader, m_name.c_str()) != 1)
    throw VCF2BinaryException(std::string("Could not open file ")+m_name+" : " + bcf_sr_strerror(m_indexed_reader->errnum) + " (VCF/BCF files must be block compressed and indexed)");
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
    if(bcf_sr_add_reader(m_indexed_reader, m_name.c_str()) != 1)
      throw VCF2BinaryException(std::string("Could not open file ")+m_name+" or its index doesn't exist - VCF/BCF files must be block compressed and indexed");
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
      auto status = bgzf_getline(hts_get_bgzfp(m_fptr), '\n', &m_vcf_file_buffer);
      m_is_record_valid = (status <= 0) ? false : true;
      if(m_is_record_valid)
        vcf_parse(&m_vcf_file_buffer, m_hdr, m_line);
    }
    else        //BCF
    {
      m_line->errcode = 0;
      //simple bcf_read
      auto status = bcf_read(m_fptr, m_hdr, m_line);
      m_is_record_valid = (status < 0) ? false : true;
      assert(m_line->errcode == 0);
    }
  }
  else  //indexed reader
  {
    auto next_line_exists = bcf_sr_next_line(m_indexed_reader);
    auto line = next_line_exists ? bcf_sr_get_line(m_indexed_reader, 0) : 0;
    if(line)
    {
      std::swap<bcf1_t*>(m_indexed_reader->readers[0].buffer[0], m_line);
      m_is_record_valid = true;
    }
    else
      m_is_record_valid = false;
  }
}

//Move constructor
VCFColumnPartition::VCFColumnPartition(VCFColumnPartition&& other)
  : File2TileDBBinaryColumnPartitionBase(std::move(other))
{
  m_local_contig_idx = other.m_local_contig_idx;
  m_contig_tiledb_column_offset = other.m_contig_tiledb_column_offset;
  m_vcf_get_buffer_size = other.m_vcf_get_buffer_size;
  m_vcf_get_buffer = other.m_vcf_get_buffer;
  m_split_output_fptr = other.m_split_output_fptr;
  other.m_vcf_get_buffer = 0;
  other.m_vcf_get_buffer_size = 0;
  other.m_split_output_fptr = 0;
}

VCFColumnPartition::~VCFColumnPartition()
{ 
  if(m_vcf_get_buffer && m_vcf_get_buffer_size)
    free(m_vcf_get_buffer);
  m_vcf_get_buffer = 0;
  m_vcf_get_buffer_size = 0;
  if(m_split_output_fptr)
    bcf_close(m_split_output_fptr);
  m_split_output_fptr = 0;
}

//VCF2Binary functions
//For VCF files
VCF2Binary::VCF2Binary(const std::string& vcf_filename, const std::vector<std::vector<std::string>>& vcf_fields,
    unsigned file_idx, VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
    size_t max_size_per_callset,
    bool treat_deletions_as_intervals,
    bool parallel_partitions, bool noupdates, bool close_file, bool discard_index)
  : File2TileDBBinaryBase(vcf_filename, file_idx, vid_mapper,
        max_size_per_callset,
        treat_deletions_as_intervals,
        parallel_partitions, noupdates, close_file)
{
  clear();
  m_vcf_fields = &vcf_fields;
  m_discard_index = discard_index;
  m_import_ID_field = false;
  m_close_file = close_file || discard_index;   //close file if index has to be discarded
  m_vcf_buffer_reader_buffer_size = 0;
  m_vcf_buffer_reader_is_bcf = false;
  m_vcf_buffer_reader_init_buffer = 0;
  m_vcf_buffer_reader_init_num_valid_bytes = 0;
  initialize(partition_bounds);
}

//For VCFBufferReader
VCF2Binary::VCF2Binary(const std::string& stream_name, const std::vector<std::vector<std::string>>& vcf_fields,
        unsigned file_idx, const int64_t buffer_stream_idx,
        VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
        const size_t vcf_buffer_reader_buffer_size, const bool vcf_buffer_reader_is_bcf,
        const uint8_t* vcf_buffer_reader_init_buffer, const size_t vcf_buffer_reader_init_num_valid_bytes,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals)
  : File2TileDBBinaryBase(stream_name,
      file_idx, buffer_stream_idx,
      vid_mapper,
      max_size_per_callset,
      treat_deletions_as_intervals,
      true, true, false) //must use parallel partitions for buffered reader
{
  clear();
  m_vcf_fields = &vcf_fields;
  //The next parameter is irrelevant for buffered readers
  m_discard_index = false;
  m_import_ID_field = false;
  //VCFBufferReader relevant params
  m_vcf_buffer_reader_buffer_size = vcf_buffer_reader_buffer_size;
  m_vcf_buffer_reader_is_bcf = vcf_buffer_reader_is_bcf;
  m_vcf_buffer_reader_init_buffer = vcf_buffer_reader_init_buffer;
  m_vcf_buffer_reader_init_num_valid_bytes = vcf_buffer_reader_init_num_valid_bytes;
  initialize(partition_bounds);
}

//Move constructor
VCF2Binary::VCF2Binary(VCF2Binary&& other)
  : File2TileDBBinaryBase(std::move(other))
{
  m_vcf_fields = other.m_vcf_fields;
  m_discard_index = other.m_discard_index;
  m_import_ID_field = other.m_import_ID_field;
  m_local_contig_idx_to_global_contig_idx = std::move(other.m_local_contig_idx_to_global_contig_idx);
  m_local_field_idx_to_global_field_idx = std::move(other.m_local_field_idx_to_global_field_idx);
  m_vcf_buffer_reader_buffer_size = other.m_vcf_buffer_reader_buffer_size;
  m_vcf_buffer_reader_is_bcf = other.m_vcf_buffer_reader_is_bcf;
  //Not useful, but copying to be safe
  m_vcf_buffer_reader_init_buffer = other.m_vcf_buffer_reader_init_buffer;
  m_vcf_buffer_reader_init_num_valid_bytes = other.m_vcf_buffer_reader_init_num_valid_bytes;
}

VCF2Binary::~VCF2Binary()
{
  clear();
}

void VCF2Binary::clear()
{
  m_local_contig_idx_to_global_contig_idx.clear();
  m_local_field_idx_to_global_field_idx.clear();
}

GenomicsDBImportReaderBase* VCF2Binary::create_new_reader_object(const std::string& filename, bool open_file) const
{
  //either reading from file or buffer parameters initialized
  assert(m_get_data_from_file || (m_vcf_buffer_reader_init_buffer && m_vcf_buffer_reader_init_num_valid_bytes && m_vcf_buffer_reader_buffer_size));
  return (m_get_data_from_file ? dynamic_cast<GenomicsDBImportReaderBase*>(new VCFReader())
      : dynamic_cast<GenomicsDBImportReaderBase*>(new VCFBufferReader(m_vcf_buffer_reader_buffer_size, m_vcf_buffer_reader_is_bcf,
       m_vcf_buffer_reader_init_buffer,  m_vcf_buffer_reader_init_num_valid_bytes))
      );
}

void VCF2Binary::initialize(const std::vector<ColumnRange>& partition_bounds)
{
  assert(m_vid_mapper);
  //Initialize partition info
  initialize_base_column_partitions(partition_bounds);
  //Setup local-global mappings
  //Get reader from column partition struct if parallel partitions, else use common base reader ptr
  auto base_reader_ptr = (m_parallel_partitions && m_base_partition_ptrs.size()) ? m_base_partition_ptrs[0]->get_base_reader_ptr()
    : m_base_reader_ptr;
  VERIFY_OR_THROW(base_reader_ptr && (std::string("Could not find valid VCF reader for ")+m_filename).c_str());
  assert(dynamic_cast<VCFReader*>(base_reader_ptr) || dynamic_cast<VCFBufferReader*>(base_reader_ptr));
  auto hdr = dynamic_cast<VCFReaderBase*>(base_reader_ptr)->get_header();
  VERIFY_OR_THROW(hdr && (std::string("Could not find valid VCF header for ")+m_filename).c_str());
  //Callset mapping
  //Length might be more than what's available in hdr due to JSON error
  m_local_callset_idx_to_tiledb_row_idx.resize(bcf_hdr_nsamples(hdr), -1ll);
  //Contig mapping
  m_local_contig_idx_to_global_contig_idx = std::move(std::vector<int>(hdr->n[BCF_DT_CTG], -1ll));
  for(auto i=0;i<hdr->n[BCF_DT_CTG];++i)
    m_vid_mapper->get_global_contig_idx(bcf_hdr_id2name(hdr, i), m_local_contig_idx_to_global_contig_idx[i]); 
  //Field mapping
  m_local_field_idx_to_global_field_idx = std::move(std::vector<int>(hdr->n[BCF_DT_ID], -1));
  for(auto i=0;i<hdr->n[BCF_DT_ID];++i)
    m_vid_mapper->get_global_field_idx(bcf_hdr_int2id(hdr, BCF_DT_ID, i), m_local_field_idx_to_global_field_idx[i]);
  int ID_field_idx = -1;
  m_import_ID_field = m_vid_mapper->get_global_field_idx("ID", ID_field_idx);
}

void VCF2Binary::initialize_column_partitions(const std::vector<ColumnRange>& partition_bounds)
{
  //Initialize reader, if needed
  if(!m_parallel_partitions)
  {
    auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(m_base_reader_ptr);
    assert(vcf_reader_ptr);
    vcf_reader_ptr->initialize(m_filename.c_str(), *m_vcf_fields, m_vid_mapper, !m_close_file);
  }
  for(auto i=0u;i<partition_bounds.size();++i)
  {
    auto vcf_column_partition_ptr = dynamic_cast<VCFColumnPartition*>(m_base_partition_ptrs[i]);
    assert(vcf_column_partition_ptr);
    //If parallel partitions, each interval gets its own reader
    if(m_parallel_partitions)
    {
      auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(vcf_column_partition_ptr->m_base_reader_ptr);
      assert(vcf_reader_ptr);
      vcf_reader_ptr->initialize(m_filename.c_str(), *m_vcf_fields, m_vid_mapper, !m_close_file);
    }
    //Indicates that nothing has been read for this interval
    vcf_column_partition_ptr->m_local_contig_idx = -1;
    //Since every line in the VCF contains data from all samples, use the sample with the smallest tiledb row idx
    //for filling the m_min_current_tiledb_row_idx field
    auto min_row_idx = INT64_MAX;
    //All enabled samples are present in every record
    vcf_column_partition_ptr->m_current_enabled_local_callset_idx_vec = m_enabled_local_callset_idx_vec;
    for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
    {
      auto tiledb_row_idx = m_local_callset_idx_to_tiledb_row_idx[m_enabled_local_callset_idx_vec[i]];
      if(tiledb_row_idx < min_row_idx)
      {
        vcf_column_partition_ptr->m_min_current_tiledb_row_idx = tiledb_row_idx;
        min_row_idx = tiledb_row_idx;
      }
    }
  }
}

bool VCF2Binary::convert_record_to_binary(std::vector<uint8_t>& buffer, File2TileDBBinaryColumnPartitionBase& partition_info)
{
  auto buffer_full = false;
  auto& vcf_partition = dynamic_cast<VCFColumnPartition&>(partition_info);
  //Cast to VCFReaderBase
  auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(partition_info.get_base_reader_ptr());
  assert(vcf_reader_ptr);
  auto* line = vcf_reader_ptr->get_line();
  assert(line);
  bcf_unpack(line, BCF_UN_ALL);
  for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
  {
    buffer_full = buffer_full || convert_VCF_to_binary_for_callset(buffer, vcf_partition, m_max_size_per_callset, i);
    if(buffer_full)
      break;
  }
  return buffer_full;
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

inline void VCF2Binary::update_local_contig_idx(VCFColumnPartition& vcf_partition, const bcf_hdr_t* hdr, bcf1_t* line)
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
  vcf_partition.m_current_column_position = vcf_partition.m_contig_tiledb_column_offset + static_cast<int64_t>(line->pos);
  bcf_set_end_point_from_info(hdr, line);
  vcf_partition.m_current_end_position = vcf_partition.m_contig_tiledb_column_offset + static_cast<int64_t>(line->m_end_point);
  //By default, bcf_set_end_point_from_info treats deletions as intervals
  if(!m_treat_deletions_as_intervals && (line->d.var_type & VCF_INDEL))
    vcf_partition.m_current_end_position = vcf_partition.m_current_column_position;
}

bool VCF2Binary::seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool& is_read_buffer_exhausted,
    bool force_seek, bool advance_reader)
{
  auto& vcf_partition = static_cast<VCFColumnPartition&>(partition_info);
  bcf_hdr_t* hdr = 0;
  if(!m_get_data_from_file) //handle VCFBufferReader
  {
    //Cast to VCFBufferReader
    auto vcf_reader_ptr = dynamic_cast<VCFBufferReader*>(partition_info.get_base_reader_ptr());
    assert(vcf_reader_ptr);
    hdr = vcf_reader_ptr->get_header();
    //advance or nothing in the buffer has been deserialized yet
    auto advance_flag = (advance_reader || vcf_reader_ptr->get_offset() == 0u);
    if(advance_flag)
    {
      vcf_reader_ptr->read_and_advance();
      auto line = vcf_reader_ptr->get_line();
      if(line)
        update_local_contig_idx(vcf_partition, hdr, line);
    }
    //read buffer done
    is_read_buffer_exhausted = !(vcf_reader_ptr->contains_unread_data());
    if(vcf_reader_ptr->get_line())
      return true;
    else
      return false; //no valid line and the buffer had no valid data at all, this stream is done
  }
  else //VCF file
  {
    is_read_buffer_exhausted = false;
    //Cast to VCFReader
    auto vcf_reader_ptr = dynamic_cast<VCFReader*>(partition_info.get_base_reader_ptr());
    assert(vcf_reader_ptr);
    hdr = vcf_reader_ptr->get_header();
    //If valid contig, i.e., continuing from a valid previous position
    //Common case: m_local_contig_idx >=0 and !force_seek and advance_reader  
    if(vcf_partition.m_local_contig_idx >= 0)
    {
      if(force_seek)
        vcf_reader_ptr->seek_read_advance(bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx),
            vcf_partition.get_column_position_in_record()-vcf_partition.m_contig_tiledb_column_offset,
            m_discard_index);
      else
        if(advance_reader)
          vcf_reader_ptr->read_and_advance();
      auto* line = vcf_reader_ptr->get_line();
      //Next line exists && (either index is stored in memory or continuing in the same contig) - common case
      if(line && (!m_discard_index || line->rid == vcf_partition.m_local_contig_idx))
      {
        update_local_contig_idx(vcf_partition, hdr, line);
        if(vcf_partition.get_column_position_in_record() <= vcf_partition.m_column_interval_end)
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
        vcf_reader_ptr->seek_read_advance(bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx), contig_pos, m_discard_index);
        auto* line = vcf_reader_ptr->get_line();
        //no more valid lines found, exit
        if(line == 0)
          return false;
        else      //valid VCF record
        {
          update_local_contig_idx(vcf_partition, hdr, line);
          //and is within this column interval, return true
          if(vcf_partition.get_column_position_in_record() <= vcf_partition.m_column_interval_end)
            return true;
          else
            return false;
        }
      }
    }
    return false;
  }
}

template<class FieldType>
bool VCF2Binary::convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition, 
    int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
    const std::string& field_name, unsigned field_type_idx)
{
  //Cast to VCFReader
  auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(vcf_partition.get_base_reader_ptr());
  assert(vcf_reader_ptr);
  auto* hdr = vcf_reader_ptr->get_header();
  auto* line = vcf_reader_ptr->get_line();
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
  auto bcf_ht_type = is_GT_field ? BCF_HT_INT : bcf_hdr_id2type(hdr, field_type_idx, field_idx);
  auto field_length = bcf_hdr_id2number(hdr, field_type_idx, field_idx);
  //Flag field lengths are set to 0 in the header :(
  field_length = (bcf_ht_type == BCF_HT_FLAG && field_length == 0) ? 1 : field_length;
  //The weirdness of VCF - string fields are marked as fixed length fields of size 1 (*facepalm*)
  auto is_vcf_str_type = ((bcf_ht_type == BCF_HT_CHAR && length_descriptor != BCF_VL_FIXED) 
      || bcf_ht_type == BCF_HT_STR) && !is_GT_field;
  length_descriptor =  is_vcf_str_type ? BCF_VL_VAR : length_descriptor;
  int max_num_values = vcf_partition.m_vcf_get_buffer_size/sizeof(FieldType);
  auto num_values = (field_type_idx == BCF_HL_INFO) ?
    bcf_get_info_values(hdr, line, field_name.c_str(), reinterpret_cast<void**>(&(vcf_partition.m_vcf_get_buffer)), &max_num_values, bcf_ht_type)
    : bcf_get_format_values(hdr, line, field_name.c_str(), reinterpret_cast<void**>(&(vcf_partition.m_vcf_get_buffer)), &max_num_values, bcf_ht_type);
  if(static_cast<uint64_t>(max_num_values)*sizeof(FieldType) >  vcf_partition.m_vcf_get_buffer_size)
    vcf_partition.m_vcf_get_buffer_size = static_cast<uint64_t>(max_num_values)*sizeof(FieldType);
  auto buffer_full = false;
  if(num_values <= 0) //Curr line does not have this field, or flag is not set
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
    auto* ptr = reinterpret_cast<const FieldType*>(vcf_partition.m_vcf_get_buffer);
    //For format fields, the ptr should point to where data for the current callset begins
    if(field_type_idx == BCF_HL_FMT)
    {
      assert(num_values%bcf_hdr_nsamples(hdr) == 0);
      num_values = num_values/bcf_hdr_nsamples(hdr);
      ptr += (local_callset_idx*num_values);
    }
    //Exclude trailing null characters for strings
    if(is_vcf_str_type)
      num_values = strnlen(reinterpret_cast<const char*>(ptr), num_values);
    auto field_length_offset = buffer_offset;
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
    auto print_sep = true;
    for(auto k=0;k<num_values;++k)
    {
      auto val = (bcf_ht_type == BCF_HT_FLAG) ? static_cast<char>(1) : ptr[k];
      //For variable length fields, terminate loop if vector_end seen
      if(is_bcf_vector_end_value<FieldType>(val) && length_descriptor != BCF_VL_FIXED)
      {
        //Update field length
#ifdef PRODUCE_CSV_CELLS
        if(!is_vcf_str_type)
#endif
        {
          buffer_full = tiledb_buffer_print<int>(buffer, field_length_offset, buffer_offset_limit, k);
          assert(!buffer_full);
        }
        break;
      }
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
  //Cast to VCFReader
  auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(vcf_partition.get_base_reader_ptr());
  assert(vcf_reader_ptr);
  auto* hdr = vcf_reader_ptr->get_header();
  auto* line = vcf_reader_ptr->get_line();
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
  int64_t column_idx = vcf_partition.get_column_position_in_record();
  buffer_full = buffer_full || tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, column_idx);
  if(buffer_full) return true;
#ifdef PRODUCE_BINARY_CELLS
  //For binary cells, must print size of cell here. Keep track of offset and update later
  auto cell_size_offset = buffer_offset;
  buffer_offset += sizeof(size_t);
#endif
  //END position
  auto end_column_idx = vcf_partition.get_end_position_in_record();
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
  //ID (if needed)
  if(m_import_ID_field)
  {
    auto ID_length = strlen(line->d.id);
    if(ID_length > 0 && (ID_length != 1 || line->d.id[0] != '.'))
    {
#ifdef PRODUCE_BINARY_CELLS
      buffer_full = buffer_full || tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, ID_length);
      if(buffer_full) return true;
#endif
      buffer_full = buffer_full || tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, line->d.id);
      if(buffer_full) return true;
    }
#ifdef PRODUCE_BINARY_CELLS
    else
      buffer_full = buffer_full || tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, 0);
#endif
  }
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
        case BCF_HT_FLAG:
          buffer_full = buffer_full || convert_field_to_tiledb<char>(buffer, vcf_partition, buffer_offset, buffer_offset_limit, local_callset_idx,
              field_name, field_type_idx);
          if(buffer_full) return true;
          break;
        default: //FIXME: handle other types
          throw VCF2BinaryException(std::string("Unhandled VCF data type ")+std::to_string(bcf_hdr_id2type(hdr, BCF_DT_ID, field_idx)));
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
  if(static_cast<size_t>(buffer_offset)+sizeof(char) <= static_cast<size_t>(buffer_offset_limit))
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

//Print partitions of the file - useful when splitting files into partitions

bool VCF2Binary::open_partition_output_file(const std::string& results_directory, std::string& output_filename,
    const std::string& output_type, File2TileDBBinaryColumnPartitionBase& partition_info, const unsigned partition_idx)
{
  auto copy_output_type = output_type;
  output_filename = m_vid_mapper->get_split_file_path(m_filename, results_directory, copy_output_type, partition_idx);
  auto valid_output_types = std::unordered_set<std::string>({ "b", "z" });
  if(valid_output_types.find(copy_output_type) == valid_output_types.end())
    throw VCF2BinaryException(std::string("Unknown output file type ")+copy_output_type+" for VCF/BCF file "+m_filename
        +"; only compressed BCF or VCF supported");
  auto& vcf_partition = static_cast<VCFColumnPartition&>(partition_info);
  vcf_partition.m_split_filename = output_filename;
  vcf_partition.m_split_output_fptr = bcf_open(output_filename.c_str(), (std::string("w")+copy_output_type).c_str());
  if(vcf_partition.m_split_output_fptr == 0)
    return false;
  auto status = bcf_hdr_write(vcf_partition.m_split_output_fptr, vcf_partition.get_header());
  if(status != 0)
    throw VCF2BinaryException(std::string("Error writing VCF header to output split file ")+output_filename+" for partition "+std::to_string(partition_idx));
  return true;
}

void VCF2Binary::write_partition_data(File2TileDBBinaryColumnPartitionBase& partition_info)
{
  auto& vcf_partition = static_cast<VCFColumnPartition&>(partition_info);
  auto hdr = vcf_partition.get_header();
  auto vcf_reader_ptr = dynamic_cast<VCFReaderBase*>(vcf_partition.m_base_reader_ptr);
  auto is_read_buffer_exhausted = false;
  //If file is re-opened, seek to position from which to begin reading, but do not advance iterator from previous position
  //The second parameter is useful if the file handler is open, but the buffer was full in a previous call
  auto has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, m_close_file, false);
  while(has_data)
  {
    bcf_write(vcf_partition.m_split_output_fptr, hdr, vcf_reader_ptr->get_line());
    has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, false, true);
  }
}

void VCF2Binary::close_partition_output_file(File2TileDBBinaryColumnPartitionBase& partition_info)
{
  auto& vcf_partition = static_cast<VCFColumnPartition&>(partition_info);
  if(vcf_partition.m_split_output_fptr == 0)
    return;
  auto format_value = vcf_partition.m_split_output_fptr->format.format;
  bcf_close(vcf_partition.m_split_output_fptr);
  vcf_partition.m_split_output_fptr = 0;
  auto status = -1;
  //Index the split file
  switch(format_value)
  {
    case htsExactFormat::bcf:
    case htsExactFormat::binary_format:
      status = bcf_index_build(vcf_partition.m_split_filename.c_str(), 14); //CSI index, min_shift value is what bcftools provides
      break;
    case htsExactFormat::vcf:
    case htsExactFormat::text_format:
      status = tbx_index_build(vcf_partition.m_split_filename.c_str(), 0, &tbx_conf_vcf); //tabix index
      break;
    default:
      break; //do nothing
  }
  if(status != 0)
    std::cerr << "WARNING: indexing of partition file "<< vcf_partition.m_split_filename <<" failed\n";
}

#endif //ifdef HTSDIR
