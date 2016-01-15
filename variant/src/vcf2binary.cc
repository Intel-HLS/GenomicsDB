#include "vcf2binary.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw VCF2BinaryException(#X);

//Binary gets priority
#if defined(PRODUCE_BINARY_CELLS) and defined(PRODUCE_CSV_CELLS)
#undef PRODUCE_CSV_CELLS
#endif

//If none defined, produce binary cells
#if !defined(PRODUCE_BINARY_CELLS) and !defined(PRODUCE_CSV_CELLS)
#define PRODUCE_BINARY_CELLS 1
#endif

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
    bool treat_deletions_as_intervals, bool parallel_partitions, bool noupdates, bool close_file)
{
  m_reader = 0;
  clear();
  m_vcf_filename = vcf_filename;
  m_vcf_fields = &vcf_fields;
  m_file_idx = file_idx;
  m_vid_mapper = &vid_mapper;
  m_max_size_per_callset = max_size_per_callset;
  m_treat_deletions_as_intervals = treat_deletions_as_intervals;
  m_parallel_partitions = parallel_partitions;
  m_noupdates = noupdates;
  m_close_file = close_file;
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
}

VCF2Binary::~VCF2Binary()
{
  if(!m_parallel_partitions)
    for(auto& x : m_partitions)
      x.m_reader = 0;
  clear();
  if(m_reader)
    bcf_sr_destroy(m_reader);
  m_reader = 0;
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

bcf_srs_t* VCF2Binary::initialize_reader(bool open_file)
{
  auto* reader = bcf_sr_init();
  bcf_sr_set_regions(reader, m_regions.c_str(), 0);   //random value
  if(open_file)
    bcf_sr_add_reader(reader, m_vcf_filename.c_str());
  return reader;
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
    m_reader = initialize_reader(!m_close_file);
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
    column_interval_info.m_reader = initialize_reader(!m_close_file);
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
      bcf_sr_add_reader(m_reader, m_vcf_filename.c_str());
    for(auto partition_idx=0u;partition_idx<partition_batches.size();++partition_idx)
    {
      auto& curr_file_batch = partition_batches[partition_idx].get_partition_file_batch(m_file_idx);
      assert(static_cast<size_t>(curr_file_batch.get_buffer_idx()) < buffer_vec.size());
      read_next_batch(*(buffer_vec[curr_file_batch.get_buffer_idx()]), m_partitions[partition_idx], curr_file_batch, close_file);
    }
    //Close file handles if needed
    if(close_file)
      bcf_sr_remove_reader(m_reader, 0);
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
    bcf_sr_add_reader(vcf_partition.m_reader, m_vcf_filename.c_str());
  auto* hdr = bcf_sr_get_header(vcf_partition.m_reader, 0);
  //Setup buffer offsets first 
  for(auto i=0ull;i<m_enabled_local_callset_idx_vec.size();++i)
  {
    auto curr_offset = partition_file_batch.get_offset_for_local_callset_idx(i, m_max_size_per_callset);
    vcf_partition.m_begin_buffer_offset_for_local_callset[i] = curr_offset;
    vcf_partition.m_buffer_offset_for_local_callset[i] = curr_offset;
    vcf_partition.m_last_full_line_end_buffer_offset_for_local_callset[i] = curr_offset;
  }
  //If file is re-opened, seek to position from which to begin reading, but do not advance iterator from previous position
  //The second parameter is useful if the file handler is open, buffer was full in a previous call
  auto has_data = seek_and_fetch_position(vcf_partition, m_close_file, false);
  auto buffer_full = false;
  auto read_one_line_fully = false;
  while(has_data && !buffer_full)
  {
    auto* line = bcf_sr_get_line(vcf_partition.m_reader, 0);
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
    bcf_sr_remove_reader(vcf_partition.m_reader, 0);
  //Caller must set to true to force fetch for this partition-file batch again
  partition_file_batch.m_fetch = false;
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
  auto hdr = bcf_sr_get_header(vcf_partition.m_reader, 0);
  //If valid contig, i.e., continuing from a valid previous position
  //Common case: m_local_contig_idx >=0 and !force_seek 
  if(vcf_partition.m_local_contig_idx >= 0)
  {
    if(force_seek)
      bcf_sr_seek(vcf_partition.m_reader, bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx), vcf_partition.m_contig_position);
    if(advance_reader || force_seek)
      bcf_sr_next_line(vcf_partition.m_reader);
    auto* line = bcf_sr_get_line(vcf_partition.m_reader, 0);
    //Next line exists
    if(line)
    {
      update_local_contig_idx(vcf_partition, line);
      if(vcf_partition.m_contig_tiledb_column_offset+static_cast<int64_t>(line->pos) <= vcf_partition.m_column_interval_end)
        return true;
      else
        return false;   //past column end
    }
    else        //no more lines found in file
      return false;
    //Either no line found or switched to new contig, either way need to seek
    //vcf_partition.m_local_contig_idx = -1;
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
      bcf_sr_seek(vcf_partition.m_reader, bcf_hdr_id2name(hdr, vcf_partition.m_local_contig_idx), contig_pos);
      bcf_sr_next_line(vcf_partition.m_reader);
      auto* line = bcf_sr_get_line(vcf_partition.m_reader, 0);
      //no valid line found, keep iterating
      if(line == 0)
        vcf_partition.m_local_contig_idx = -1;
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
  tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
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
  tiledb_buffer_print<const std::string&>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
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
  tiledb_buffer_print<char>(buffer, buffer_offset, buffer_offset_limit, NULL_VALUE);
}
#endif

template<class FieldType>
bool VCF2Binary::convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition, 
    int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
    const std::string& field_name, unsigned field_type_idx)
{
  auto* hdr = bcf_sr_get_header(vcf_partition.m_reader, 0);
  auto* line = bcf_sr_get_line(vcf_partition.m_reader, 0);
  //FIXME: avoid strings
  auto is_GT_field = (field_type_idx == BCF_HL_FMT && field_name == "GT");
  assert(line);
  assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size()
      && local_callset_idx < bcf_hdr_nsamples(hdr));
  auto field_idx = bcf_hdr_id2int(hdr, BCF_DT_ID, field_name.c_str());
  //FIXME: handle missing fields gracefully
  VERIFY_OR_THROW(field_idx >= 0);
  //FIXME: special length descriptors
  auto length_descriptor = is_GT_field ? BCF_VL_P : bcf_hdr_id2length(hdr, field_type_idx, field_idx);
  auto field_length = bcf_hdr_id2number(hdr, field_type_idx, field_idx);
  auto bcf_ht_type = is_GT_field ? BCF_HT_INT : bcf_hdr_id2type(hdr, field_type_idx, field_idx);
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
      buffer_full = tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, 0);
    else        //fixed length field - fill with NULL values
      for(auto i=0;i<field_length;++i)
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
      buffer_full = tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, num_values);
      if(buffer_full) return true;
    }
    auto* ptr = reinterpret_cast<const FieldType*>(vcf_partition.m_vcf_get_buffer);
    //For format fields, the ptr should point to where data for the current callset begins
    if(field_type_idx == BCF_HL_FMT)
    {
      assert(num_values%bcf_hdr_nsamples(hdr) == 0);
      num_values = num_values/bcf_hdr_nsamples(hdr);
      ptr += (local_callset_idx*num_values);
    }
    for(auto k=0;k<num_values;++k)
    {
      auto val = ptr[k];
      if(is_GT_field)
        val = bcf_gt_allele(static_cast<int>(val));
      buffer_full = buffer_full || tiledb_buffer_print<FieldType>(buffer, buffer_offset, buffer_offset_limit, val);
      if(buffer_full) return true;
    }
  }
  return buffer_full;
}

bool VCF2Binary::convert_VCF_to_binary_for_callset(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
    size_t size_per_callset, uint64_t enabled_callsets_idx)
{
  auto* hdr = bcf_sr_get_header(vcf_partition.m_reader, 0);
  auto* line = bcf_sr_get_line(vcf_partition.m_reader, 0);
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
  assert(line_begin_buffer_offset >= begin_buffer_offset && line_begin_buffer_offset <= begin_buffer_offset + size_per_callset);
  assert(buffer_offset >= begin_buffer_offset && buffer_offset <= begin_buffer_offset + size_per_callset);
  assert(buffer_offset >= line_begin_buffer_offset);
  const int64_t buffer_offset_limit = begin_buffer_offset + size_per_callset;
  bool buffer_full = false;
  auto curr_cell_begin_offset = buffer_offset;
  //Row
  int64_t row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
  if(row_idx < 0) return false;
  buffer_full = buffer_full || tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, row_idx, false);
  if(buffer_full) return true;
//#ifdef PRODUCE_CSV_CELLS
  ////For CSV cells, no comma at the beginning of the cell - shift left and reduce offset
  //memmove(&(buffer[curr_cell_begin_offset]), &(buffer[curr_cell_begin_offset+sizeof(char)]), buffer_offset-curr_cell_begin_offset-sizeof(char));
  //buffer_offset -= sizeof(char);
//#endif
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
    assert(line->d.flt[i] < m_local_field_idx_to_global_field_idx.size());
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
      //FIXME: handle missing fields gracefully
      VERIFY_OR_THROW(field_idx >= 0);
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
