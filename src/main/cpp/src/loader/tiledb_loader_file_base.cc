#include "tiledb_loader_file_base.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw File2TileDBBinaryException(#X);

//Move constructor
File2TileDBBinaryColumnPartitionBase::File2TileDBBinaryColumnPartitionBase(File2TileDBBinaryColumnPartitionBase&& other)
{
  m_column_interval_begin = other.m_column_interval_begin;
  m_column_interval_end = other.m_column_interval_end;
  m_begin_buffer_offset_for_local_callset = std::move(other.m_begin_buffer_offset_for_local_callset);
  m_last_full_line_end_buffer_offset_for_local_callset =
    std::move(other.m_last_full_line_end_buffer_offset_for_local_callset);
  m_buffer_offset_for_local_callset = std::move(other.m_buffer_offset_for_local_callset);
  m_buffer_full_for_local_callset = std::move(other.m_buffer_full_for_local_callset);
  m_split_filename = std::move(other.m_split_filename);
  //Move and nullify other
  m_base_reader_ptr = other.m_base_reader_ptr;
  other.m_base_reader_ptr = 0;
  m_buffer_ptr = other.m_buffer_ptr;
  other.m_buffer_ptr = 0;
}

File2TileDBBinaryColumnPartitionBase::~File2TileDBBinaryColumnPartitionBase()
{
  if(m_base_reader_ptr)
    delete m_base_reader_ptr;
  m_base_reader_ptr = 0;
  clear();
}

void File2TileDBBinaryColumnPartitionBase::initialize_base_class_members(const int64_t begin, const int64_t end,
    const uint64_t num_enabled_callsets, GenomicsDBImportReaderBase* reader_ptr)
{
  m_column_interval_begin = begin;
  m_column_interval_end = end;
  //buffer offsets for each callset in this partition
  m_buffer_offset_for_local_callset.resize(num_enabled_callsets);
  m_begin_buffer_offset_for_local_callset.resize(num_enabled_callsets);
  m_last_full_line_end_buffer_offset_for_local_callset.resize(num_enabled_callsets);
  m_buffer_full_for_local_callset.resize(num_enabled_callsets, false);
  m_base_reader_ptr = reader_ptr;
}

//File2TileDBBinaryBase functions

#ifdef PRODUCE_BINARY_CELLS

template<class FieldType>
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const FieldType val, bool print_sep)
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
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const char* val, bool print_sep)
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
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, char* val, bool print_sep)
{
  return tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
}

//specialization for std::string type
template<>
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const std::string& val, bool print_sep)
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
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, std::string& val, bool print_sep)
{
  return tiledb_buffer_print<const std::string&>(buffer, buffer_offset, buffer_offset_limit, val, print_sep);
}

template<class FieldType>
bool File2TileDBBinaryBase::tiledb_buffer_print_null(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit)
{
  return tiledb_buffer_print<FieldType>(buffer, buffer_offset, buffer_offset_limit, get_tiledb_null_value<FieldType>());
}

#endif //ifdef PRODUCE_BINARY_CELLS

#ifdef PRODUCE_CSV_CELLS
template<class FieldType>
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const FieldType val, bool print_sep)
{
  std::stringstream ss;
  if(print_sep)
    ss << "," << val;
  else
    ss << val;
  std::string x = ss.str();
  size_t add_size = x.length()*sizeof(char);
  //Do not write anything past the limit
  if(static_cast<size_t>(buffer_offset) + add_size > static_cast<size_t>(buffer_offset_limit))
    return true;
  memcpy(&(buffer[buffer_offset]), x.c_str(), add_size);
  buffer_offset += add_size;
  return false;
}

template<class FieldType>
bool File2TileDBBinaryBase::tiledb_buffer_print_null(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit)
{
  //return tiledb_buffer_print<char>(buffer, buffer_offset, buffer_offset_limit, '\0');
  //Print nothing for null fields when producing a CSV - except for a separator
  return tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, "");
}
#endif

template<>
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit, const std::string val, bool print_sep)
{
  return tiledb_buffer_print<const std::string&>(buffer, buffer_offset, buffer_offset_limit, static_cast<const std::string&>(val), print_sep);
}

//Template instantiations
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const char val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const int val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const unsigned val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const int64_t val, bool print_sep);
#ifdef __MACH__
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const size_t val, bool print_sep);
#endif
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const uint64_t val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const float val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const double val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const char* val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit, const std::string& val, bool print_sep);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<char>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<int>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<unsigned>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<int64_t>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
#ifdef __MACH__
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<size_t>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
#endif
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<uint64_t>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<float>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);
template
bool File2TileDBBinaryBase::tiledb_buffer_print_null<double>(std::vector<uint8_t>& buffer, int64_t& buffer_offset,
    const int64_t buffer_offset_limit);

//Constructor
File2TileDBBinaryBase::File2TileDBBinaryBase(const std::string& filename,
    unsigned file_idx, const int64_t buffer_stream_idx,
    VidMapper& vid_mapper,
    size_t max_size_per_callset,
    bool treat_deletions_as_intervals,
    bool parallel_partitions, bool noupdates, bool close_file)
{
  m_filename = filename;
  m_file_idx = file_idx;
  m_buffer_stream_idx = buffer_stream_idx;
  m_vid_mapper = &(vid_mapper);
  m_max_size_per_callset = max_size_per_callset;
  m_treat_deletions_as_intervals = treat_deletions_as_intervals;
  m_parallel_partitions = parallel_partitions;
  m_noupdates = noupdates;
  m_close_file = close_file;
  m_get_data_from_file = (buffer_stream_idx < 0) ? true : false;
  m_no_mandatory_VCF_fields = false;
  //Callset mapping
  vid_mapper.get_local_tiledb_row_idx_vec(filename, m_local_callset_idx_to_tiledb_row_idx);
  m_local_callset_idx_to_enabled_idx.resize(m_local_callset_idx_to_tiledb_row_idx.size(), -1ll);
  m_enabled_local_callset_idx_vec.clear();
  for(auto i=0ull;i<m_local_callset_idx_to_tiledb_row_idx.size();++i)
    if(m_local_callset_idx_to_tiledb_row_idx[i] >= 0)
    {
      m_local_callset_idx_to_enabled_idx[i] = m_enabled_local_callset_idx_vec.size();
      m_enabled_local_callset_idx_vec.push_back(i);
    }
  std::sort(m_enabled_local_callset_idx_vec.begin(), m_enabled_local_callset_idx_vec.end(),
      LocalCallSetIdxCompareByTileDBRowIdx(m_local_callset_idx_to_tiledb_row_idx));
  m_base_reader_ptr = 0;
  m_histogram = 0;
  auto* GT_field_info_ptr = vid_mapper.get_field_info("GT");
  m_store_phase_information_for_GT = (GT_field_info_ptr && GT_field_info_ptr->m_length_descriptor.contains_phase_information());
}

void File2TileDBBinaryBase::initialize_base_column_partitions(const std::vector<ColumnRange>& partition_bounds)
{
  m_base_reader_ptr = !m_parallel_partitions ? create_new_reader_object(m_filename, !m_close_file) : 0;
  //Initialize base partition pointers
  m_base_partition_ptrs.resize(partition_bounds.size(), 0);
  for(auto i=0ull;i<m_base_partition_ptrs.size();++i)
  {
    m_base_partition_ptrs[i] = create_new_column_partition_object();
    const auto& column_partition = partition_bounds[i];
    //If parallel partitions, each interval gets its own reader
    auto base_reader_ptr = m_parallel_partitions ? create_new_reader_object(m_filename, !m_close_file) : m_base_reader_ptr;
    //Initialize base class members
    m_base_partition_ptrs[i]->initialize_base_class_members(
        column_partition.first, column_partition.second,
        m_enabled_local_callset_idx_vec.size(), base_reader_ptr);
  }
  //Sub-class virtual function
  initialize_column_partitions(partition_bounds);
}

void File2TileDBBinaryBase::copy_simple_members(const File2TileDBBinaryBase& other)
{
  m_parallel_partitions = other.m_parallel_partitions;
  m_noupdates = other.m_noupdates;
  m_close_file = other.m_close_file;
  m_treat_deletions_as_intervals = other.m_treat_deletions_as_intervals;
  m_get_data_from_file = other.m_get_data_from_file;
  m_no_mandatory_VCF_fields = other.m_no_mandatory_VCF_fields;
  m_store_phase_information_for_GT = other.m_store_phase_information_for_GT;
  m_file_idx = other.m_file_idx;
  m_buffer_stream_idx = other.m_buffer_stream_idx;
  m_max_size_per_callset = other.m_max_size_per_callset;
}

File2TileDBBinaryBase::File2TileDBBinaryBase(File2TileDBBinaryBase&& other)
{
  copy_simple_members(other);
  m_vid_mapper = other.m_vid_mapper;
  other.m_vid_mapper = 0;
  m_filename = std::move(other.m_filename);
  m_local_callset_idx_to_tiledb_row_idx = std::move(other.m_local_callset_idx_to_tiledb_row_idx);
  m_enabled_local_callset_idx_vec = std::move(other.m_enabled_local_callset_idx_vec);
  m_local_callset_idx_to_enabled_idx = std::move(other.m_local_callset_idx_to_enabled_idx);
  m_base_reader_ptr = other.m_base_reader_ptr;
  other.m_base_reader_ptr = 0;
  m_base_partition_ptrs = std::move(other.m_base_partition_ptrs);
  m_histogram = other.m_histogram;
  other.m_histogram = 0;
}

File2TileDBBinaryBase::~File2TileDBBinaryBase()
{
  for(auto& x : m_base_partition_ptrs)
  {
    if(!m_parallel_partitions)
      x->m_base_reader_ptr = 0;
    delete x;
    x = 0;
  }
  if(m_base_reader_ptr)
    delete m_base_reader_ptr;
  m_base_reader_ptr = 0;
  m_vid_mapper = 0;
  clear();
  if(m_histogram)
    delete m_histogram;
  m_histogram = 0;
}

void File2TileDBBinaryBase::clear()
{
  m_filename.clear();
  m_local_callset_idx_to_tiledb_row_idx.clear();
  m_enabled_local_callset_idx_vec.clear();
  m_local_callset_idx_to_enabled_idx.clear();
  m_base_partition_ptrs.clear();
}

void File2TileDBBinaryBase::read_next_batch(std::vector<std::vector<uint8_t>*>& buffer_vec,
    std::vector<ColumnPartitionBatch>& partition_batches,
    std::vector<BufferStreamIdentifier>& exhausted_buffer_stream_identifiers, size_t& num_exhausted_buffer_streams,
    bool close_file)
{
  if(m_parallel_partitions)
  {
#pragma omp parallel for
    for(auto partition_idx=0u;partition_idx<partition_batches.size();++partition_idx)
    {
      auto& curr_file_batch = partition_batches[partition_idx].get_partition_file_batch(m_file_idx);
      assert(static_cast<size_t>(curr_file_batch.get_buffer_idx()) < buffer_vec.size());
      read_next_batch(*(buffer_vec[curr_file_batch.get_buffer_idx()]), *(m_base_partition_ptrs[partition_idx]),
          curr_file_batch, partition_idx,
          exhausted_buffer_stream_identifiers, num_exhausted_buffer_streams,
          close_file);
    }
  }
  else
  {
    //Open file handles if needed
    if(m_close_file)
      m_base_reader_ptr->add_reader();
    for(auto partition_idx=0u;partition_idx<partition_batches.size();++partition_idx)
    {
      auto& curr_file_batch = partition_batches[partition_idx].get_partition_file_batch(m_file_idx);
      assert(static_cast<size_t>(curr_file_batch.get_buffer_idx()) < buffer_vec.size());
      read_next_batch(*(buffer_vec[curr_file_batch.get_buffer_idx()]), *(m_base_partition_ptrs[partition_idx]),
          curr_file_batch, partition_idx,
          exhausted_buffer_stream_identifiers, num_exhausted_buffer_streams,
          close_file);
    }
    //Close file handles if needed
    if(close_file)
      m_base_reader_ptr->remove_reader();
  }
  m_close_file = close_file;
}

void File2TileDBBinaryBase::read_next_batch(std::vector<uint8_t>& buffer,
    File2TileDBBinaryColumnPartitionBase& partition_info,
    ColumnPartitionFileBatch& partition_file_batch, const unsigned partition_idx,
    std::vector<BufferStreamIdentifier>& exhausted_buffer_stream_identifiers, size_t& num_exhausted_buffer_streams,
    bool close_file)
{
  //Nothing to do
  if(!partition_file_batch.m_fetch || partition_file_batch.m_completed)
    return;
  //Open file handles if needed
  if(m_parallel_partitions && m_close_file)
    partition_info.m_base_reader_ptr->add_reader();
  //Setup buffer offsets first
  for(auto i=0ull;i<static_cast<size_t>(partition_file_batch.get_num_orders());++i)
  {
    auto curr_offset = partition_file_batch.get_offset_for_local_callset_idx(i, m_max_size_per_callset);
    partition_info.m_begin_buffer_offset_for_local_callset[i] = curr_offset;
    partition_info.m_buffer_offset_for_local_callset[i] = curr_offset;
    partition_info.m_last_full_line_end_buffer_offset_for_local_callset[i] = curr_offset;
  }
  auto is_read_buffer_exhausted = false;
  //If file is re-opened, seek to position from which to begin reading, but do not advance iterator from previous position
  //The second parameter is useful if the file handler is open, but the buffer was full in a previous call
  auto has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, m_close_file, false);
  auto buffer_full = false;
  auto read_one_line_fully = false;
  //Set buffer full to false
  for(auto i=0ull;i<partition_info.m_buffer_full_for_local_callset.size();++i)
    partition_info.m_buffer_full_for_local_callset[i] = false;
  partition_info.m_buffer_ptr = &(buffer);
  while(has_data && !buffer_full)
  {
    buffer_full = convert_record_to_binary(buffer, partition_info);
    if(!buffer_full)
    {
      //Store buffer offsets at the beginning of the line
      for(auto i=0ull;i<static_cast<size_t>(partition_file_batch.get_num_orders());++i)
        partition_info.m_last_full_line_end_buffer_offset_for_local_callset[i] = partition_info.m_buffer_offset_for_local_callset[i];
      read_one_line_fully = true;
      //For buffered readers, if the read buffer is empty return control to caller
      if(is_read_buffer_exhausted)
      {
        assert(m_buffer_stream_idx >= 0); //must be valid buffer stream
        exhausted_buffer_stream_identifiers[num_exhausted_buffer_streams++] = std::move(BufferStreamIdentifier(m_buffer_stream_idx, partition_idx));
        break;
      }
      has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, false, true);  //no need to re-seek, use next_line() directly, advance file pointer
    }
    else
      VERIFY_OR_THROW(read_one_line_fully && "Buffer did not have space to hold a line fully - increase buffer size")
  }
  //put Tiledb NULL for row_idx as end-of-batch marker
  for(auto i=0ull;i<static_cast<size_t>(partition_file_batch.get_num_orders());++i)
  {
#ifdef PRODUCE_BINARY_CELLS
    tiledb_buffer_print_null<int64_t>(buffer, partition_info.m_last_full_line_end_buffer_offset_for_local_callset[i],
        partition_info.m_begin_buffer_offset_for_local_callset[i] + m_max_size_per_callset);
#endif
#ifdef PRODUCE_CSV_CELLS
    tiledb_buffer_print<char>(buffer, partition_info.m_last_full_line_end_buffer_offset_for_local_callset[i],
        partition_info.m_begin_buffer_offset_for_local_callset[i] + m_max_size_per_callset, '\0', false);
#endif
  }
  if(!has_data)
    partition_file_batch.m_completed = true;
  //Close file handles if needed
  if(m_parallel_partitions && close_file)
    partition_info.m_base_reader_ptr->remove_reader();
  //Advance write idx for partition_file_batch
  partition_file_batch.advance_write_idx();
}

void File2TileDBBinaryBase::create_histogram(uint64_t max_histogram_range, unsigned num_bins)
{
  if(m_histogram)
    delete m_histogram;
  m_histogram = new UniformHistogram(0, max_histogram_range, num_bins);
  //Open file handles if needed
  if(m_close_file && !m_parallel_partitions)
    m_base_reader_ptr->add_reader();
  //Should really have 1 partition only
  for(auto* base_partition_ptr : m_base_partition_ptrs)
  {
    auto& partition_info = *base_partition_ptr;
    if(m_close_file && m_parallel_partitions)
      partition_info.m_base_reader_ptr->add_reader();
    auto is_read_buffer_exhausted = false;
    auto has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, m_close_file, false);
    while(has_data && !is_read_buffer_exhausted)
    {
      auto column_idx = partition_info.get_column_position_in_record();
      auto num_callsets = get_num_callsets_in_record(partition_info);
      for(auto i=0ull;i<num_callsets;++i)
        m_histogram->add_interval(column_idx, column_idx);
      has_data = seek_and_fetch_position(partition_info, is_read_buffer_exhausted, false, true);
    }
    if(m_close_file && m_parallel_partitions)
      partition_info.m_base_reader_ptr->remove_reader();
  }
  if(m_close_file && !m_parallel_partitions)
    m_base_reader_ptr->remove_reader();
}

void File2TileDBBinaryBase::print_all_partitions(const std::string& results_directory, const std::string& output_type, const int rank,
    const bool close_file)
{
  //Open file handles if needed
  if(!m_parallel_partitions && m_close_file)
    m_base_reader_ptr->add_reader();
#pragma omp parallel for if(m_parallel_partitions)
  for(auto partition_idx=0u;partition_idx<m_base_partition_ptrs.size();++partition_idx)
    print_partition(*(m_base_partition_ptrs[partition_idx]), results_directory, output_type, rank < 0 ? partition_idx : rank, close_file);
  //Close file handles if needed
  if(!m_parallel_partitions && close_file)
    m_base_reader_ptr->remove_reader();
  m_close_file = close_file;
}

void File2TileDBBinaryBase::print_partition(File2TileDBBinaryColumnPartitionBase& partition_info,
    const std::string& results_directory, const std::string& output_type,
    const unsigned partition_idx, const bool close_file)
{
  //Open file handles if needed
  if(m_parallel_partitions && m_close_file)
    partition_info.m_base_reader_ptr->add_reader();
  std::string output_filename = "";
  auto status = open_partition_output_file(results_directory, output_filename, output_type, partition_info, partition_idx);
  if(!status)
    throw File2TileDBBinaryException(std::string("Could not open partition output file ")+output_filename);
  write_partition_data(partition_info);
  close_partition_output_file(partition_info);
  //Close file handles if needed
  if(m_parallel_partitions && close_file)
    partition_info.m_base_reader_ptr->remove_reader();
}
