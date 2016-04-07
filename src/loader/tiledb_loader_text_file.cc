/**
 * The MIT License (MIT)
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

#include "tiledb_loader_text_file.h"
#include "vcf.h"
#include "variant_field_data.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw LineBasedTextFileException(#X);

std::string g_tmp_scratch_dir = "/tmp";

LineBasedTextFileReader::LineBasedTextFileReader()
  : FileReaderBase()
{
  m_filename.clear();
  m_fptr = 0;
  m_line_buffer_size = 4096u;    //4KB
  m_line_buffer = new char[m_line_buffer_size];
  m_line_length = 0;
}

LineBasedTextFileReader::~LineBasedTextFileReader()
{
  m_filename.clear();
  if(m_fptr)
    fclose(m_fptr);
  m_fptr = 0;
  if(m_line_buffer && m_line_buffer_size)
    delete[] m_line_buffer;
  m_line_buffer = 0;
  m_line_buffer_size = 0;
  m_line_length = 0;
}

void LineBasedTextFileReader::initialize(const char* filename, bool open_file)
{
  m_filename = filename;
  add_reader();
  if(!open_file)
  {
    fclose(m_fptr);
    m_fptr = 0;
  }
}

void LineBasedTextFileReader::add_reader()
{
  if(m_fptr)
    return;
  m_fptr = fopen(m_filename.c_str(), "r");
  if(m_fptr == 0)
    throw LineBasedTextFileException(std::string("Could not open file: ")+m_filename);
}

void LineBasedTextFileReader::remove_reader()
{
  if(m_fptr)
    fclose(m_fptr);
  m_fptr = 0; 
}

void LineBasedTextFileReader::read_and_advance()
{
  assert(m_fptr);
  if(!feof(m_fptr))
  {
    auto num_bytes_read = getline(&m_line_buffer, &m_line_buffer_size, m_fptr);
    m_is_record_valid = (num_bytes_read < 0) ? false : true;
    //m_line_length = (num_bytes_read > 0) ? ((m_line_buffer[num_bytes_read-1] == static_cast<uint8_t>('\n')) ? num_bytes_read-1 : num_bytes_read)
    //includes newline
    m_line_length = m_is_record_valid ? num_bytes_read : 0ull; 
  }
  else
  {
    m_is_record_valid = false;
    m_line_length = 0;
  }
}

CSV2TileDBBinary::CSV2TileDBBinary(const std::string& filename,
        unsigned file_idx, VidMapper& vid_mapper,
        size_t max_size_per_callset,
        const std::vector<ColumnRange>& partition_bounds,
        bool treat_deletions_as_intervals,
        bool parallel_partitions, bool noupdates, bool close_file)
      : LineBasedTextFile2TileDBBinary(filename, file_idx, vid_mapper,
          max_size_per_callset,
          treat_deletions_as_intervals,
          parallel_partitions, noupdates, close_file)
{
  m_cleanup_file = false;
  auto file_type = 0u;
  auto status = vid_mapper.get_file_type(filename, file_type);
  if(!status)
    throw LineBasedTextFileException(std::string("Could not find an entry for file ")+filename);
  //Sort CSV
  if(file_type == VidFileTypeEnum::UNSORTED_CSV_FILE_TYPE)
  {
    auto unsorted_fptr = fopen(filename.c_str(), "r");
    //File doesn't exist
    if(unsorted_fptr == 0)
      throw LineBasedTextFileException(std::string("Could not open file ")+filename);
    fclose(unsorted_fptr);
    auto sorted_filename = strdup((g_tmp_scratch_dir+"/sorted_csv_XXXXXX").c_str());
    auto fd = mkstemp(sorted_filename);
    if(fd == -1)
      throw LineBasedTextFileException(std::string("Could not create sorted file in temporary directory ")+g_tmp_scratch_dir+" for file "+filename);
    close(fd);
    auto cmd_string = std::string("sort -T ")+g_tmp_scratch_dir+" -t, -k2,2n -k1,1n -o "+sorted_filename+" "+filename;
    auto fptr = popen(cmd_string.c_str(), "r");
    if(fptr == 0)
      throw LineBasedTextFileException(std::string("Sort failed for file ")+filename);
    auto status = pclose(fptr);
    if(status != 0)
     throw LineBasedTextFileException(std::string("Sort failed for file ")+filename);
    m_filename = sorted_filename;
    free(sorted_filename);
    m_cleanup_file = true;
  }
  //Initialize partition info
  initialize_base_column_partitions(partition_bounds);
}

CSV2TileDBBinary::~CSV2TileDBBinary()
{
  //Cleanup if unsorted csv file
  if(m_cleanup_file)
    remove(m_filename.c_str());
  m_cleanup_file = false;
}

void CSV2TileDBBinary::initialize_column_partitions(const std::vector<ColumnRange>& partition_bounds)
{
  //Initialize reader, if needed
  if(!m_parallel_partitions)
  {
    auto csv_reader_ptr = dynamic_cast<LineBasedTextFileReader*>(m_base_reader_ptr);
    assert(csv_reader_ptr);
    csv_reader_ptr->initialize(m_filename.c_str(), !m_close_file);
  }
  for(auto i=0u;i<partition_bounds.size();++i)
  {
    auto csv_column_partition_ptr = dynamic_cast<CSV2TileDBBinaryColumnPartition*>(m_base_partition_ptrs[i]);
    assert(csv_column_partition_ptr);
    //If parallel partitions, each interval gets its own reader
    if(m_parallel_partitions)
    {
      auto csv_reader_ptr = dynamic_cast<LineBasedTextFileReader*>(csv_column_partition_ptr->get_base_reader_ptr());
      assert(csv_reader_ptr);
      csv_reader_ptr->initialize(m_filename.c_str(), !m_close_file);
    }
  }
}

void CSV2TileDBBinary::set_order_of_enabled_callsets(int64_t& order_value, std::vector<int64_t>& tiledb_row_idx_to_order) const
{
  //Point all callsets to the same value of order i.e. the order value for the first callset
  //This ensures that effectively, a single buffer space is allocated for all callsets in this 
  //file
  if(m_enabled_local_callset_idx_vec.size())
  {
    for(auto local_callset_idx : m_enabled_local_callset_idx_vec)
    {
      assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size());
      auto row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
      assert(row_idx >= 0);
      assert(static_cast<size_t>(row_idx) < tiledb_row_idx_to_order.size());
      tiledb_row_idx_to_order[row_idx] = order_value;
    }
    order_value++;
  }
}

void CSV2TileDBBinary::list_active_row_idxs(const ColumnPartitionBatch& partition_batch, int64_t& row_idx_offset, std::vector<int64_t>& row_idx_vec) const
{
  auto& partition_file_batch = partition_batch.get_partition_file_batch(m_file_idx);
  if(partition_file_batch.m_fetch && !partition_file_batch.m_completed)
  {
    //Effectively inform loader that only 1 callset in this file is ready
    //Since all callsets in this file use the same buffer space, it doesn't really matter
    if(m_enabled_local_callset_idx_vec.size())
    {
      auto local_callset_idx = m_enabled_local_callset_idx_vec[0];
      assert(static_cast<size_t>(local_callset_idx) < m_local_callset_idx_to_tiledb_row_idx.size());
      auto row_idx = m_local_callset_idx_to_tiledb_row_idx[local_callset_idx];
      assert(row_idx >= 0);
      row_idx_vec[row_idx_offset++] = row_idx;
    }
  }
}

bool CSV2TileDBBinary::seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool force_seek, bool advance_reader)
{
  auto& csv_partition_info = dynamic_cast<CSV2TileDBBinaryColumnPartition&>(partition_info);
  auto csv_reader_ptr = dynamic_cast<LineBasedTextFileReader*>(partition_info.get_base_reader_ptr());
  assert(csv_reader_ptr);
  if(force_seek || !(csv_partition_info.is_initialized_file_position_to_partition_begin()))
  {
    //Had previously sought file ptr to the column partition begin - now just use fpos_t
    if(csv_partition_info.is_initialized_file_position_to_partition_begin())
    {
      csv_reader_ptr->seek(csv_partition_info.m_file_position);
      csv_reader_ptr->read_and_advance();
    }
    else
    {
      //First read
      csv_reader_ptr->read_and_advance();
      auto line = csv_reader_ptr->get_line();
      while(line != 0)
      {
        parse_line(line, csv_partition_info, TileDBCSVFieldPosIdxEnum::TILEDB_CSV_COLUMN_POS_IDX, false);        //parse only till column idx
        if(csv_partition_info.m_current_column_position >= csv_partition_info.m_column_interval_begin)
          break;
        csv_reader_ptr->read_and_advance();
        line = csv_reader_ptr->get_line();
      }
      csv_partition_info.set_initialized_file_position_to_partition_begin(true);
    }
  }
  else
    if(advance_reader)
      csv_reader_ptr->read_and_advance();
  //Store file position
  csv_reader_ptr->get_position(csv_partition_info.m_file_position);
  auto line = csv_reader_ptr->get_line();
  if(line)
  {
    //Parse line for column idx to check whether file pointer has gone beyond column_partition_end
    parse_line(line, csv_partition_info, TileDBCSVFieldPosIdxEnum::TILEDB_CSV_COLUMN_POS_IDX, false);        //parse only till column idx
    return (csv_partition_info.m_current_column_position <= csv_partition_info.m_column_interval_end);
  }
  else
    return false;
}

template<class FieldType>
void CSV2TileDBBinary::handle_field_token(const char* token_ptr,
    CSVLineParseStruct* csv_line_parse_ptr, CSV2TileDBBinaryColumnPartition& csv_partition_info,
    std::vector<uint8_t>& buffer, int64_t& buffer_offset, const int64_t buffer_offset_limit,
    VariantFieldTypeEnum variant_field_type_enum)
{
  auto field_idx = csv_line_parse_ptr->get_field_idx(); 
  auto num_elements = csv_line_parse_ptr->get_num_elements();
  auto enabled_idx_in_file = csv_line_parse_ptr->get_enabled_idx_in_file();
  //Direct all data to the buffer corresponding to the first enabled callset in this file
  auto buffer_idx = 0u;
  //For fixed length fields, set length first
  if(!(csv_line_parse_ptr->read_num_elements()) && !(m_array_schema->is_variable_length_field(field_idx)))
  {
    num_elements = m_array_schema->val_num(field_idx);
    csv_line_parse_ptr->set_num_elements(num_elements);
  }
  //Read #elements if not already read
  if(!(csv_line_parse_ptr->read_num_elements()) && m_array_schema->is_variable_length_field(field_idx))
  {
    //string fields need to be handled weirdly
    if(variant_field_type_enum == VariantFieldTypeEnum::VARIANT_FIELD_STRING ||
        variant_field_type_enum == VariantFieldTypeEnum::VARIANT_FIELD_CHAR)
    {
#ifdef PRODUCE_BINARY_CELLS
      //Copy length into buffer
      csv_partition_info.set_buffer_full_if_true(buffer_idx,
          tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, strlen(token_ptr)));
#endif
      csv_partition_info.set_buffer_full_if_true(buffer_idx,
          tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit, token_ptr));
      //since field is already handled, set num elements to 0
      //code at the end of the function will increment field idx
      num_elements = 0u;
    }
    else
    {
      num_elements = from_string_to_tiledb<int>(token_ptr);
      csv_partition_info.set_buffer_full_if_true(buffer_idx,
          tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit, num_elements));
    }
    csv_line_parse_ptr->set_num_elements(num_elements);
  }
  else
  {
    csv_partition_info.set_buffer_full_if_true(buffer_idx,
      tiledb_buffer_print<FieldType>(buffer, buffer_offset, buffer_offset_limit,
          from_string_to_tiledb<FieldType>(token_ptr)));
    csv_line_parse_ptr->increment_field_element_idx();
  }
  if(csv_line_parse_ptr->get_field_element_idx() >= num_elements)
  {
    csv_line_parse_ptr->reset_field_element_idx();
    csv_line_parse_ptr->increment_field_idx();
  }
}

void CSV2TileDBBinary::handle_token(CSVLineParseStruct* csv_line_parse_ptr, const char* token_ptr, const size_t field_size)
{
  if(csv_line_parse_ptr->is_past_max_token_idx())
    return;
  auto& csv_partition_info = *(csv_line_parse_ptr->get_csv_column_partition_ptr());
  char* endptr = 0;
  //Special fields
  switch(csv_line_parse_ptr->get_token_idx())
  {
    case TileDBCSVFieldPosIdxEnum::TILEDB_CSV_ROW_POS_IDX:
      {
        auto row_idx = strtoll(token_ptr, &endptr, 0);
        VERIFY_OR_THROW((endptr != token_ptr) && "Could not parse row field");
        csv_line_parse_ptr->set_row_idx(row_idx);
        auto local_callset_idx = m_vid_mapper->get_idx_in_file_for_row_idx(row_idx);
        auto enabled_idx_in_file = get_enabled_idx_for_local_callset_idx(local_callset_idx);
        csv_line_parse_ptr->set_enabled_idx_in_file(enabled_idx_in_file);
        break;
      }
    case TileDBCSVFieldPosIdxEnum::TILEDB_CSV_COLUMN_POS_IDX:
      {
        csv_partition_info.m_current_column_position = strtoll(token_ptr, &endptr, 0);
        VERIFY_OR_THROW((endptr != token_ptr) && "Could not parse column field");
        break;
      }
    default:
      break;
  }
  auto enabled_idx_in_file = csv_line_parse_ptr->get_enabled_idx_in_file();
  //Direct all data to the buffer corresponding to the first enabled callset in this file
  auto buffer_idx = 0u;
  //Should store this row in buffer?
  if(csv_line_parse_ptr->should_store_in_buffer() && enabled_idx_in_file >= 0 && !(csv_partition_info.is_buffer_full(buffer_idx)))
  {
    const int64_t begin_buffer_offset = csv_partition_info.m_begin_buffer_offset_for_local_callset[buffer_idx];
    const int64_t line_begin_buffer_offset = csv_partition_info.m_last_full_line_end_buffer_offset_for_local_callset[buffer_idx];
    int64_t& buffer_offset = csv_partition_info.m_buffer_offset_for_local_callset[buffer_idx];
    assert(line_begin_buffer_offset >= begin_buffer_offset && line_begin_buffer_offset <= static_cast<int64_t>(begin_buffer_offset + m_max_size_per_callset));
    assert(buffer_offset >= begin_buffer_offset && buffer_offset <= static_cast<int64_t>(begin_buffer_offset + m_max_size_per_callset));
    assert(buffer_offset >= line_begin_buffer_offset);
    const int64_t buffer_offset_limit = begin_buffer_offset + m_max_size_per_callset;
    auto& buffer = csv_partition_info.get_buffer();
    switch(csv_line_parse_ptr->get_token_idx())
    {
      case TileDBCSVFieldPosIdxEnum::TILEDB_CSV_ROW_POS_IDX:
        {
          //Do not print sep when printing row idx
          csv_partition_info.set_buffer_full_if_true(buffer_idx,
            tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, from_string_to_tiledb<int64_t>(token_ptr),
                false));
          break;
        }
      case TileDBCSVFieldPosIdxEnum::TILEDB_CSV_COLUMN_POS_IDX:
        {
          csv_partition_info.set_buffer_full_if_true(buffer_idx,
            tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit, from_string_to_tiledb<int64_t>(token_ptr)));
#ifdef PRODUCE_BINARY_CELLS
          //reserve space for cell size
          csv_line_parse_ptr->set_cell_size_offset(buffer_offset);
          buffer_offset += sizeof(size_t);
          csv_partition_info.set_buffer_full_if_true(buffer_idx, buffer_offset > buffer_offset_limit);
#endif
          break;
        }
      default:
        {
          //Depends on field idx
          auto field_idx = csv_line_parse_ptr->get_field_idx();
          switch(field_idx)
          {
            case VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_END_IDX:
              {
                csv_partition_info.set_buffer_full_if_true(buffer_idx,
                  tiledb_buffer_print<int64_t>(buffer, buffer_offset, buffer_offset_limit,
                      from_string_to_tiledb<int64_t>(token_ptr)));
                break;
              }
            case VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_REF_IDX:
            case VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_ALT_IDX:
              {
#ifdef PRODUCE_BINARY_CELLS
                csv_partition_info.set_buffer_full_if_true(buffer_idx,
                  tiledb_buffer_print<int>(buffer, buffer_offset, buffer_offset_limit,
                      strlen(token_ptr)));
#endif
                csv_partition_info.set_buffer_full_if_true(buffer_idx,
                  tiledb_buffer_print<const char*>(buffer, buffer_offset, buffer_offset_limit,
                      token_ptr));
                break;
              }
            case VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_QUAL_IDX:
              {
                csv_partition_info.set_buffer_full_if_true(buffer_idx,
                  tiledb_buffer_print<float>(buffer, buffer_offset, buffer_offset_limit,
                      from_string_to_tiledb<float>(token_ptr)));
                break;
              }
            case VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_FILTER_IDX:
              {
                handle_field_token<int>(token_ptr,
                    csv_line_parse_ptr, csv_partition_info,
                    buffer, buffer_offset, buffer_offset_limit,
                    VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(std::type_index(typeid(int))));
                break;
              }
            default:    //other optional fields
              {
                auto variant_field_type_enum = VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(
                    m_array_schema->type(field_idx));
                switch(variant_field_type_enum)
                {
                  case VariantFieldTypeEnum::VARIANT_FIELD_INT:
                    {
                      handle_field_token<int>(token_ptr,
                          csv_line_parse_ptr, csv_partition_info,
                          buffer, buffer_offset, buffer_offset_limit,
                          variant_field_type_enum);
                      break;
                    }
                  case VariantFieldTypeEnum::VARIANT_FIELD_FLOAT:
                    {
                      handle_field_token<float>(token_ptr,
                          csv_line_parse_ptr, csv_partition_info,
                          buffer, buffer_offset, buffer_offset_limit,
                          variant_field_type_enum);
                      break;
                    }
                  case VariantFieldTypeEnum::VARIANT_FIELD_CHAR:
                  case VariantFieldTypeEnum::VARIANT_FIELD_STRING:
                    {
                      handle_field_token<std::string>(token_ptr,
                          csv_line_parse_ptr, csv_partition_info,
                          buffer, buffer_offset, buffer_offset_limit,
                          variant_field_type_enum);
                      break;
                    }
                  default:
                    throw LineBasedTextFileException(std::string("Type ")+m_array_schema->type(field_idx).name()+
                        " not handled by the CSV importer");
                    break;
                }
              }
          }
          //Increment field idx for fixed mandatory fields
          if(field_idx <= VariantArraySchemaFixedFieldsEnum::VARIANT_ARRAY_SCHEMA_QUAL_IDX)
            csv_line_parse_ptr->increment_field_idx();
          break;
        }
    }
  }
}

void CSV2TileDBBinary::handle_end_of_line(CSVLineParseStruct* csv_line_parse_ptr)
{
  auto& csv_partition_info = *(csv_line_parse_ptr->get_csv_column_partition_ptr());
  auto enabled_idx_in_file = csv_line_parse_ptr->get_enabled_idx_in_file();
  //Direct all data to the buffer corresponding to the first enabled callset in this file
  auto buffer_idx = 0u;
  if(csv_line_parse_ptr->should_store_in_buffer() && enabled_idx_in_file >= 0 && !(csv_partition_info.is_buffer_full(buffer_idx)))
  {
    //If a fixed length field is at the end of the schema, then it's possible that the last field or last element
    //of the fixed length field is omitted
    if(csv_line_parse_ptr->get_field_idx() < m_array_schema->attribute_num())
      handle_token(csv_line_parse_ptr, "", 0);
    assert(csv_line_parse_ptr->get_field_idx() == m_array_schema->attribute_num());
    const int64_t begin_buffer_offset = csv_partition_info.m_begin_buffer_offset_for_local_callset[buffer_idx];
    const int64_t line_begin_buffer_offset =
      csv_partition_info.m_last_full_line_end_buffer_offset_for_local_callset[buffer_idx];
    int64_t& buffer_offset = csv_partition_info.m_buffer_offset_for_local_callset[buffer_idx];
    const int64_t buffer_offset_limit = begin_buffer_offset + m_max_size_per_callset;
#ifdef PRODUCE_BINARY_CELLS
    //Write cell size into buffer
    auto cell_size_offset = csv_line_parse_ptr->get_cell_size_offset();
    assert(cell_size_offset >= 0);
    auto buffer_full = tiledb_buffer_print<size_t>(csv_partition_info.get_buffer(), cell_size_offset, buffer_offset_limit,
        buffer_offset-line_begin_buffer_offset);
    assert(!buffer_full);
#endif
#ifdef PRODUCE_CSV_CELLS
    //Add newline
    csv_partition_info.set_buffer_full_if_true(buffer_idx,
        tiledb_buffer_print<char>(csv_partition_info.get_buffer(), buffer_offset, buffer_offset_limit, '\n'));
#endif
  }
}

bool CSV2TileDBBinary::parse_line(const char* line, CSV2TileDBBinaryColumnPartition& csv_partition_info, const unsigned max_token_idx, const bool store_in_buffer)
{
  auto csv_reader_ptr = dynamic_cast<LineBasedTextFileReader*>(csv_partition_info.get_base_reader_ptr());
  assert(csv_reader_ptr);
  CSVLineParseStruct parse_obj(this, &csv_partition_info, max_token_idx, store_in_buffer);
#ifdef USE_LIBCSV
  csv_parse(&(csv_partition_info.m_csv_parser), line, csv_reader_ptr->get_line_length(),
      csv_parse_callback, csv_line_end_callback,
      reinterpret_cast<void*>(&(parse_obj)));
#endif
  //Direct all data to the buffer corresponding to the first enabled callset in this file
  auto buffer_idx = 0u;
  return (store_in_buffer && csv_partition_info.is_buffer_full(buffer_idx));
}

void csv_parse_callback(void* token_ptr, size_t field_size, void* parse_ptr)
{
  auto csv_line_parse_ptr = reinterpret_cast<CSVLineParseStruct*>(parse_ptr);
  if(!(csv_line_parse_ptr->is_past_max_token_idx()))
    csv_line_parse_ptr->get_csv2tiledb_binary_ptr()->handle_token(csv_line_parse_ptr, reinterpret_cast<const char*>(token_ptr), field_size);
  csv_line_parse_ptr->increment_token_idx();
}

void csv_line_end_callback(int terminating_token, void* parse_ptr)
{
  auto csv_line_parse_ptr = reinterpret_cast<CSVLineParseStruct*>(parse_ptr);
  csv_line_parse_ptr->get_csv2tiledb_binary_ptr()->handle_end_of_line(csv_line_parse_ptr);
}
