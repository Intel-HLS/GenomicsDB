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

#include "vid_mapper.h"
#include "genomicsdb_multid_vector_field.h"
#include "tiledb_loader_file_base.h"

template<class ElementType>
ElementType str_to_element(const char* str, const size_t element_begin_idx,
    const size_t current_element_length);

//Template specializations
template<>
int64_t str_to_element(const char* str, const size_t element_begin_idx,
    const size_t current_element_length)
{
  char* endptr = 0;
  //Why copy? No strntol function :(
  //Optimization - no heap operations if small enough
  const size_t array_size = 64u;
  char array[array_size];
  if(current_element_length+1u < array_size) //+1 for NULL char
  {
    memcpy_s(array, current_element_length, str+element_begin_idx, current_element_length);
    array[current_element_length] = '\0';
    return strtoll(array, &endptr, 0);
  }
  else
  {
    std::string copy(str+element_begin_idx, current_element_length);
    return strtoll(copy.c_str(), &endptr, 0);
  }
}

template<>
int str_to_element(const char* str, const size_t element_begin_idx,
    const size_t current_element_length)
{
  return (current_element_length > 0u
      && !(current_element_length == 3u && strncasecmp(str+element_begin_idx, "NaN", 3u) == 0))
    ? str_to_element<int64_t>(str, element_begin_idx, current_element_length)
    : get_bcf_missing_value<int>();
}

template<>
float str_to_element(const char* str, const size_t element_begin_idx,
    const size_t current_element_length)
{
  if(current_element_length == 0u
      || (current_element_length == 3u && strncasecmp(str+element_begin_idx, "NaN", 3u) == 0))
    return get_bcf_missing_value<float>();
  char* endptr = 0;
  //Why copy? No strntol function :(
  //Optimization - no heap operations if small enough
  const size_t array_size = 32u;
  char array[array_size];
  if(current_element_length+1u < array_size) //+1 for NULL char
  {
    memcpy_s(array, current_element_length, str+element_begin_idx, current_element_length);
    array[current_element_length] = '\0';
    return strtof(array, &endptr);
  }
  else
  {
    std::string copy(str+element_begin_idx, current_element_length);
    return strtof(copy.c_str(), &endptr);
  }
}

template<class T>
void cast_join_and_print(std::ostream& fptr, const uint8_t* ptr, const size_t idx,
    const size_t num_elements, const char sep)
{
  if(num_elements)
  {
    auto data_ptr = reinterpret_cast<const T*>(ptr);
    auto val = data_ptr[idx];
    if(!is_bcf_missing_value<T>(val))
      fptr << val;
    for(auto i=idx+1ull;i<num_elements;++i)
    {
      fptr << sep;
      auto val = data_ptr[i];
      if(!is_bcf_missing_value<T>(val))
        fptr << val;
    }
  }
}

void GenomicsDBMultiDVectorIdx::advance_to_index_in_next_dimension(const size_t idx)
{
  //<size_of_data starting at F> <F[0] vector data> <F[1] vector data> <#elements in F> <offset of F[0], offset of F[1], ...> ...
  auto& length_descriptor = m_field_info_ptr->m_length_descriptor;
  assert(m_current_dim_idx+1u < length_descriptor.get_num_dimensions());
  //the innermost dimension is simply a vector of elements - the size, #elements and offsets are NOT stored on disk,
  //Since the innermost dimension is simply a raw vector, hence the +2u check for dimensions which store 
  //size, #elements, offsets on disk
  if(m_current_dim_idx+2u < length_descriptor.get_num_dimensions())
  {
    auto size_of_data_at_current_index = *(reinterpret_cast<const uint64_t*>(m_ro_field_ptr));
    auto num_entries_and_offsets_ptr = reinterpret_cast<const uint64_t*>(m_ro_field_ptr
        + sizeof(uint64_t) //the size_at_current_ptr in the previous statement
        + size_of_data_at_current_index);
    m_num_entries_in_current_dim = *num_entries_and_offsets_ptr;
    //Below: Use +1, not sizeof(uint64_t), since num_entries_and_offsets_ptr is uint64_t*
    m_offsets_ptr = num_entries_and_offsets_ptr + 1u;
    assert(idx < m_num_entries_in_current_dim);
    m_ro_field_ptr = m_ro_field_ptr
      + sizeof(uint64_t) //the size_at_current_ptr position
      + m_offsets_ptr[idx];
  }
  else
  {
    //Update #entries based on size of innermost vector 
    assert(m_field_info_ptr->get_element_size() > 0u);
    m_num_entries_in_current_dim = get_size_of_current_index()/(m_field_info_ptr->get_element_size());
    m_offsets_ptr = 0;
    assert(idx < m_num_entries_in_current_dim);
    m_ro_field_ptr = m_ro_field_ptr + idx*(m_field_info_ptr->get_element_size());
  }
  m_current_index_in_current_dimension = idx;
  ++m_current_dim_idx;
}

void GenomicsDBMultiDVectorIdx::advance_index_in_current_dimension()
{
  //assert(static_cast<size_t>(m_current_dim_idx) < m_field_info_ptr->m_length_descriptor.get_num_dimensions());
  //assert(m_current_index_in_current_dimension < m_num_entries_in_current_dim);
  ////the innermost dimension is simply a vector of elements - the size, #elements and offsets are NOT stored on disk,
  ////Since the innermost dimension is simply a raw vector, hence the +2u check for dimensions which store
  ////size, #elements, offsets on disk
  //if(m_current_dim_idx+1u < m_field_info_ptr->m_length_descriptor.get_num_dimensions())
  //{
    //auto update_offset = m_offsets_ptr[m_current_index_in_current_dimension+1u]
      //- m_offsets_ptr[m_current_index_in_current_dimension];
    //m_ro_field_ptr = m_ro_field_ptr + update_offset;
  //}
  //else
    //m_ro_field_ptr = m_ro_field_ptr + (m_field_info_ptr->get_element_size());
  //++m_current_index_in_current_dimension;
  set_index_in_current_dimension(m_current_index_in_current_dimension+1u);
}

void GenomicsDBMultiDVectorIdx::set_index_in_current_dimension(const uint64_t idx)
{
  assert(static_cast<size_t>(m_current_dim_idx) < m_field_info_ptr->m_length_descriptor.get_num_dimensions());
  assert(idx <= m_num_entries_in_current_dim);
  auto add_to_ptr = (idx >= m_current_index_in_current_dimension);
  //the innermost dimension is simply a vector of elements - the size, #elements and offsets are NOT stored on disk,
  //Since the innermost dimension is simply a raw vector, hence the +2u check for dimensions which store 
  //size, #elements, offsets on disk
  auto update_offset = 0ull;
  if(m_current_dim_idx+1u < m_field_info_ptr->m_length_descriptor.get_num_dimensions())
    update_offset = add_to_ptr
      ?  (m_offsets_ptr[idx] - m_offsets_ptr[m_current_index_in_current_dimension])
      :  (m_offsets_ptr[m_current_index_in_current_dimension] - m_offsets_ptr[idx]);
  else
    update_offset = (
        add_to_ptr
        ? (idx-m_current_index_in_current_dimension)
        : (m_current_index_in_current_dimension-idx)
        )
      *(m_field_info_ptr->get_element_size());
  m_ro_field_ptr = add_to_ptr ? (m_ro_field_ptr + update_offset)
      : (m_ro_field_ptr - update_offset);
  m_current_index_in_current_dimension = idx;
}

template<class ElementType>
bool parse_and_store_tuple_element(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const char* str, const uint64_t begin, const uint64_t length)
{
  auto value = str_to_element<ElementType>(str, begin, length);
  File2TileDBBinaryBase::tiledb_buffer_resize_if_needed_and_print<ElementType>(
      buffer, reinterpret_cast<int64_t&>(write_offset), value);
  return is_bcf_valid_value<ElementType>(value);
}

bool GenomicsDBMultiDVectorFieldParseAndStoreOperator::parse_and_store_tuple_element_int(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const char* str, const uint64_t begin, const uint64_t length,
    const unsigned tuple_element_idx) const
{
  return parse_and_store_tuple_element<int>(buffer, write_offset, str, begin, length);
}

bool GenomicsDBMultiDVectorFieldParseAndStoreOperator::parse_and_store_tuple_element_float(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const char* str, const uint64_t begin, const uint64_t length,
    const unsigned tuple_element_idx) const
{
  return parse_and_store_tuple_element<float>(buffer, write_offset, str, begin, length);
}

bool GenomicsDBMultiDVectorFieldParseDivideUpAndStoreOperator::parse_and_store_tuple_element_int(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const char* str, const uint64_t begin, const uint64_t length,
    const unsigned tuple_element_idx) const
{
  auto status = parse_and_store_tuple_element<int>(buffer, write_offset, str, begin, length);
  auto offset_written_to = write_offset - sizeof(int);
  auto ptr = reinterpret_cast<int*>(&(buffer[offset_written_to]));
  assert(tuple_element_idx < m_tuple_indexes_to_divide_bitset.size());
  if(m_tuple_indexes_to_divide_bitset[tuple_element_idx])
  {
    auto val = *ptr;
    *ptr = (val/m_divisor + ((m_curr_idx < val%m_divisor) ? 1 : 0));
  }
  return status;
}

bool GenomicsDBMultiDVectorFieldParseDivideUpAndStoreOperator::parse_and_store_tuple_element_float(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const char* str, const uint64_t begin, const uint64_t length,
    const unsigned tuple_element_idx) const
{
  auto status = parse_and_store_tuple_element<float>(buffer, write_offset, str, begin, length);
  auto offset_written_to = write_offset - sizeof(int);
  auto ptr = reinterpret_cast<float*>(&(buffer[offset_written_to]));
  assert(tuple_element_idx < m_tuple_indexes_to_divide_bitset.size());
  if(m_tuple_indexes_to_divide_bitset[tuple_element_idx])
  {
    auto val = *ptr;
    *ptr = (val/m_divisor);
  }
  return status;
}

template<class ElementType>
void fill_with_bcf_missing_values(std::vector<uint8_t>& buffer,
    uint64_t& write_offset,
    const uint64_t curr_num, const uint64_t n)
{
  for(auto i=curr_num;i<n;++i)
    File2TileDBBinaryBase::tiledb_buffer_resize_if_needed_and_print<ElementType>(
        buffer, reinterpret_cast<int64_t&>(write_offset), get_bcf_missing_value<ElementType>());
}

/*
   Given a delimited string, parse and store the data in a binary buffer as described in
   the header.
   1,2,3,4|5,6|7,8$9,10|11|12
*/
std::vector<uint64_t> GenomicsDBMultiDVectorField::parse_and_store_numeric(
    std::vector<std::vector<uint8_t>>& buffer_vec, //outer vector - one for each element of tuple
    const FieldInfo& field_info,
    const char* str, const size_t str_length,
    const GenomicsDBMultiDVectorFieldParseAndStoreOperator& op)
{
  auto& length_descriptor = field_info.m_length_descriptor;
  auto num_elements_in_tuple = field_info.get_genomicsdb_type().get_num_elements_in_tuple();
  if(num_elements_in_tuple > buffer_vec.size())
    buffer_vec.resize(num_elements_in_tuple);
  for(auto& buffer : buffer_vec)
    buffer.resize(4096u); //4KiB
  auto r_idx = 0ull; //read idx
  //#bytes for curr data in dim i
  std::vector<std::vector<uint64_t>> dim_sizes_vec(num_elements_in_tuple,
      std::vector<uint64_t>(length_descriptor.get_num_dimensions(), 0ull));
  //Offsets for the data in each dim F[0] offset, F[1] offset etc
  //Each dimension begins with offset 0
  std::vector<std::vector<std::vector<uint64_t>>> dim_offsets_vec(
      num_elements_in_tuple,
      std::vector<std::vector<uint64_t>>(length_descriptor.get_num_dimensions(),
        std::vector<uint64_t>(1u, 0ull))
      );
  //Specifies the current offset in buffer at which the data for dim i should be written
  //the last 2 dimensions don't get any sizes
  std::vector<std::vector<uint64_t>> dim_write_begin_offsets_vec(
      num_elements_in_tuple,
      std::vector<uint64_t>(length_descriptor.get_num_dimensions())
      );
  //The first dimension begins writing at offset 0
  //A 64-bit uint64_t is allocated for each dimension's size
  for(auto tuple_element_idx=0u;tuple_element_idx<num_elements_in_tuple;++tuple_element_idx)
    for(auto i=0ull;i<dim_write_begin_offsets_vec[tuple_element_idx].size();++i)
      dim_write_begin_offsets_vec[tuple_element_idx][i] = i*sizeof(uint64_t);
  //Points to the begin of the current element - length
  auto current_element_begin_read_idx = 0ull; //r stands for read
  auto current_element_length = 0ull;
  std::vector<uint64_t> num_elements_in_innermost_dim_read_vec(
      num_elements_in_tuple, 0ull);
  auto max_num_elements_in_innermost_dim = 0ull;
  std::vector<uint64_t> total_size_of_multi_d_data_vec(num_elements_in_tuple, 0ull);
  //Idx of entry inside a tuple
  auto curr_element_idx_in_tuple_parsed = 0ull;
  //Whether there is at least one valid element in the tuple
  auto curr_tuple_contains_one_valid_element = false;
  //TODO: can be vectorized, probably not worth it since only used in loading
  while(r_idx <= str_length)
  {
    auto sep_dim_idx = 0;
    auto found_delim = false;
    if(r_idx < str_length)
    {
      //is this a delimiter char?
      for(sep_dim_idx=static_cast<int>(length_descriptor.get_num_dimensions()-1u);
          sep_dim_idx>=0;--sep_dim_idx)
      {
        if(str[r_idx] == length_descriptor.get_vcf_delimiter(sep_dim_idx))
        {
          found_delim = true;
          break;
        }
      }
    }
    else
    {
      found_delim = true;
      sep_dim_idx = -1; //flush out offsets for outermost dimension (dim_idx 0)
    }
    if(found_delim)
    {
      //Write out element
      auto is_valid_element = false;
      switch(field_info.get_genomicsdb_type().get_tuple_element_bcf_ht_type(curr_element_idx_in_tuple_parsed))
      {
        case BCF_HT_INT:
            is_valid_element = op.parse_and_store_tuple_element_int(buffer_vec[curr_element_idx_in_tuple_parsed],
                dim_write_begin_offsets_vec[curr_element_idx_in_tuple_parsed].back(),
                str, current_element_begin_read_idx, current_element_length,
                curr_element_idx_in_tuple_parsed);
          break;
        case BCF_HT_REAL:
            is_valid_element = op.parse_and_store_tuple_element_float(buffer_vec[curr_element_idx_in_tuple_parsed],
                dim_write_begin_offsets_vec[curr_element_idx_in_tuple_parsed].back(),
                str, current_element_begin_read_idx, current_element_length,
                curr_element_idx_in_tuple_parsed);
          break;
        default:
          throw GenomicsDBMultiDVectorFieldOperatorException(
              "Cannot parse multi-D fields whose elements are neither int nor float");
      }
      ++num_elements_in_innermost_dim_read_vec[curr_element_idx_in_tuple_parsed];
      max_num_elements_in_innermost_dim = std::max<uint64_t>(
          max_num_elements_in_innermost_dim,
          num_elements_in_innermost_dim_read_vec[curr_element_idx_in_tuple_parsed]);
      //Can this info be used somehow? unclear
      curr_tuple_contains_one_valid_element = curr_tuple_contains_one_valid_element
        || is_valid_element;
      //Delimiter not of the innermost dimension
      if(sep_dim_idx+1u < length_descriptor.get_num_dimensions())
      {
        //For each element of the tuple
        for(auto tuple_element_idx=0u;tuple_element_idx<num_elements_in_tuple;++tuple_element_idx)
        {
          //Fill out other tuple elements with bcf_missing_values to make equi-length
          auto tuple_element_size = 0ull;
          switch(field_info.get_genomicsdb_type().get_tuple_element_bcf_ht_type(tuple_element_idx))
          {
            case BCF_HT_INT:
              fill_with_bcf_missing_values<int>(buffer_vec[tuple_element_idx],
                  dim_write_begin_offsets_vec[tuple_element_idx].back(),
                  num_elements_in_innermost_dim_read_vec[tuple_element_idx],
                  max_num_elements_in_innermost_dim);
              tuple_element_size = sizeof(int);
              break;
            case BCF_HT_REAL:
              fill_with_bcf_missing_values<float>(buffer_vec[tuple_element_idx],
                  dim_write_begin_offsets_vec[tuple_element_idx].back(),
                  num_elements_in_innermost_dim_read_vec[tuple_element_idx],
                  max_num_elements_in_innermost_dim);
              tuple_element_size = sizeof(float);
              break;
            default:
              throw GenomicsDBMultiDVectorFieldOperatorException(
                  "Cannot parse multi-D fields whose elements are neither int nor float");
          }
          //Example: 1,2,3,4|5,6|7,8$9,10|11|12
          //Length of innermost dimension is the number of elements read
          auto last_dim_size =
            (max_num_elements_in_innermost_dim*tuple_element_size);
          //Convenience variables
          auto& buffer = buffer_vec[tuple_element_idx];
          auto& dim_sizes = dim_sizes_vec[tuple_element_idx];
          auto& dim_offsets = dim_offsets_vec[tuple_element_idx];
          auto& dim_write_begin_offsets = dim_write_begin_offsets_vec[tuple_element_idx];
          //Ignore innermost dimension - hence -2u instead of -1u
          for(auto i=static_cast<int>(length_descriptor.get_num_dimensions()-2u);
              i>=std::max<int>(sep_dim_idx,0);--i)  //sep_dim_idx can be -1, used to flush out offsets at string (str) end
          {
            dim_sizes[i] += last_dim_size;
            dim_offsets[i].push_back(dim_offsets[i].back()+last_dim_size);
            //Outer dim delimiter found - write out curr dim size and offsets
            //If sep_dim_idx == -1, flushes out the outermost dim
            if(i > sep_dim_idx)
            {
              int64_t write_offset = dim_write_begin_offsets[i];
              //write the dim size at dim_write_begin_offsets
              File2TileDBBinaryBase::tiledb_buffer_resize_if_needed_and_print<uint64_t>(
                  buffer,
                  write_offset,
                  dim_sizes[i]);
              write_offset = dim_write_begin_offsets[i]
                +sizeof(uint64_t)  //64-bit size - previous statement
                +dim_sizes[i];     //data bytes
              //write #entries in dim i
              File2TileDBBinaryBase::tiledb_buffer_resize_if_needed_and_print<uint64_t>(
                  buffer,
                  write_offset,
                  dim_offsets[i].size()-1u); //why -1? If there are N entries, there are N+1 offsets
              //write offsets
              auto offsets_start_offset = dim_write_begin_offsets[i]
                +sizeof(uint64_t)  //64-bit size
                +dim_sizes[i]      //data bytes
                +sizeof(uint64_t); //#elements
              auto offsets_size = dim_offsets[i].size()*sizeof(uint64_t);
              if(offsets_start_offset+offsets_size >= buffer.size())
                buffer.resize(2*(offsets_start_offset+offsets_size+1u));
              memcpy_s(&(buffer[offsets_start_offset]), offsets_size, &(dim_offsets[i][0u]), offsets_size);
              //Update last_dim_size
              last_dim_size = (dim_sizes[i]
                  + sizeof(uint64_t)    //for storing size of data at dim i
                  + sizeof(uint64_t)    //for storing #entries at dim i
                  + offsets_size        //for storing offsets at dim i
                  );
            }
          }
          if(sep_dim_idx == -1)
            total_size_of_multi_d_data_vec[tuple_element_idx] = last_dim_size;
          //sep_dim_idx can be -1
          //The outermost dimension whose delimiter was hit
          auto outermost_delim_dim_idx = static_cast<size_t>(std::max<int>(sep_dim_idx, 0));
          auto new_write_begin_offset =
            dim_write_begin_offsets[outermost_delim_dim_idx]
            + sizeof(uint64_t)                          //64-bit for storing size of outermost_delim_dim_idx
            + dim_sizes[outermost_delim_dim_idx];       //current size of outermost_delim_dim_idx
          //Update dim_write_begin_offsets for dimensions which 'completed' i.e. wrote out all the data
          //for one idx value for the dimension
          for(unsigned i=outermost_delim_dim_idx+1u,j=0ul;i<length_descriptor.get_num_dimensions();++i,++j)
          {
            dim_write_begin_offsets[i] = new_write_begin_offset + j*sizeof(uint64_t); //64-bit uint64_t for sizes
            //Reset dim sizes
            dim_sizes[i] = 0u;
            //Reset offsets vector - single element with value 0
            dim_offsets[i].resize(1u);
            dim_offsets[i][0u] = 0u;
          }
        }
        curr_element_idx_in_tuple_parsed = 0u; //reset index in tuple
        //Reset number of elements in innermost dim to 0
        num_elements_in_innermost_dim_read_vec.assign(num_elements_in_tuple, 0u);
        max_num_elements_in_innermost_dim = 0u;
      }
      else
        curr_element_idx_in_tuple_parsed = (curr_element_idx_in_tuple_parsed+1u)
          %num_elements_in_tuple; //move to next element in tuple
      current_element_begin_read_idx = r_idx+1u;
      current_element_length = 0u;
    }
    else
      ++current_element_length;
    ++r_idx;
  }
  return total_size_of_multi_d_data_vec;
}

std::vector<uint64_t> GenomicsDBMultiDVectorField::parse_and_store_numeric(const char* str, const size_t str_length)
{
  auto total_size_of_multi_d_data_vec = std::move(GenomicsDBMultiDVectorField::parse_and_store_numeric(
      m_rw_field_data, *m_field_info_ptr,
      str, str_length));
  for(auto tuple_element_idx=0u;tuple_element_idx<total_size_of_multi_d_data_vec.size();
      ++tuple_element_idx)
    m_rw_field_data[tuple_element_idx].resize(total_size_of_multi_d_data_vec[tuple_element_idx]);
  return total_size_of_multi_d_data_vec;
}


//Iterating over all dimensions and all index values can be a recursive process
//We can model this as a stack
void GenomicsDBMultiDVectorField::run_operation(GenomicsDBMultiDVectorFieldOperator& multid_vector_field_operator,
    const std::vector<const uint8_t*>& data_ptr_vec) const
{
  if(m_field_info_ptr->get_genomicsdb_type().get_num_elements_in_tuple() != data_ptr_vec.size())
    throw GenomicsDBMultiDVectorFieldOperatorException("Data ptr vec and genomicsdb_type do not have the same number of elements");
  auto& length_descriptor = m_field_info_ptr->m_length_descriptor;
  //Stack replacing the recursive function
  std::vector<std::vector<GenomicsDBMultiDVectorIdx>> idx_stack_vec;
  for(auto i=0u;i<data_ptr_vec.size();++i)
    idx_stack_vec.emplace_back(std::vector<GenomicsDBMultiDVectorIdx>(
          1u, GenomicsDBMultiDVectorIdx(data_ptr_vec[i], m_field_info_ptr, 0u)));
  //current index vector e.g. A[5][0][3] will have 5,0,3
  //changes as stack is traversed
  //Don't bother with the last dimension
  std::vector<uint64_t> curr_index_vector(length_descriptor.get_num_dimensions()-1u);
  //Optimization - tracks the outermost dimension index which changed since the
  //last call to operate(). Initialized to N-2
  auto outermost_dim_idx_changed_since_last_call_to_operate = length_descriptor.get_num_dimensions()-2u;
  //Vector to hold pointers and sizes for calling the operator
  std::vector<const uint8_t*> op_ptr_vec(data_ptr_vec.size(), 0);
  std::vector<size_t> op_size_vec(data_ptr_vec.size(), 0u);
  while(!idx_stack_vec[0].empty())
  {
    auto& top_idx = idx_stack_vec[0].back();
    auto curr_dim_idx_in_curr_index_vector = top_idx.get_current_dim_index();
    assert(static_cast<size_t>(curr_dim_idx_in_curr_index_vector) < curr_index_vector.size());
    curr_index_vector[curr_dim_idx_in_curr_index_vector++] = top_idx.get_current_index_in_current_dimension();
    //At the end of this block, the stack must look like the following:
    //[  [top+1], [top][0] ]
    //This ensures that subsequent iterations expand on lower dimensions starting with [top][0]
    //first followed by [top+1]
    auto is_top_idx_valid = (top_idx.get_current_index_in_current_dimension()
        < top_idx.get_num_entries_in_current_dimension());
    if(is_top_idx_valid)
    {
      //First index in the next dimension is added to the top of the stack
      //Why +2u? dimension N-2 corresponds to a pointer to a contiguous segment of values
      if(top_idx.get_current_dim_index()+2u < length_descriptor.get_num_dimensions())
      {
        for(auto i=0u;i<data_ptr_vec.size();++i)
        {
          auto& curr_tuple_element_top_idx = idx_stack_vec[i].back();
          auto copy_idx = curr_tuple_element_top_idx; //copy
          assert(copy_idx.get_current_index_in_current_dimension()
              < copy_idx.get_num_entries_in_current_dimension()
              && copy_idx.get_current_dim_index()+2u
              < length_descriptor.get_num_dimensions());
          //index 0 in the next dimension
          copy_idx.advance_to_index_in_next_dimension(0u);
          //Cannot move this to outside the if-else block because idx_stack memory can be reallocated
          //by the emplace_back() statement
          curr_tuple_element_top_idx.advance_index_in_current_dimension();
          idx_stack_vec[i].emplace_back(copy_idx);
        }
      }
      else
      {
        for(auto i=0u;i<data_ptr_vec.size();++i)
        {
          auto& curr_tuple_element_top_idx = idx_stack_vec[i].back();
          assert(curr_tuple_element_top_idx.get_current_dim_index()+2u
              >= length_descriptor.get_num_dimensions());
          op_ptr_vec[i] = curr_tuple_element_top_idx.get_ptr<uint8_t>();
          op_size_vec[i] = curr_tuple_element_top_idx.get_size_of_current_index();
        }
        multid_vector_field_operator.operate(op_ptr_vec,
            op_size_vec, curr_index_vector,
            outermost_dim_idx_changed_since_last_call_to_operate);
        //reset outermost_dim_idx_changed_since_last_call_to_operate 
        outermost_dim_idx_changed_since_last_call_to_operate = length_descriptor.get_num_dimensions()-2u; 
        //Cannot move this to outside the if-else block because idx_stack memory can be reallocated (if block)
        for(auto i=0u;i<data_ptr_vec.size();++i)
          idx_stack_vec[i].back().advance_index_in_current_dimension();
      }
      //DO NOT USE top_idx in this if block after this - emplace_back() might have reallocated the vector
    }
    else
    {
      for(auto i=0u;i<idx_stack_vec.size();++i)
        idx_stack_vec[i].pop_back();
      --outermost_dim_idx_changed_since_last_call_to_operate;
    }
  }
}

//GenomicsDBMultiDVectorFieldOperator functions

GenomicsDBMultiDVectorFieldVCFPrinter::GenomicsDBMultiDVectorFieldVCFPrinter(std::ostream& fptr,
    const FieldInfo& field_info)
{
  m_first_call = true;
  m_fptr = &fptr;
  m_field_info_ptr = &field_info;
}

void GenomicsDBMultiDVectorFieldVCFPrinter::operate(const std::vector<const uint8_t*>& ptr_vec,
    const std::vector<size_t>& size_of_data_vec,
    const std::vector<uint64_t>& idx_vector, int outermost_dim_idx_changed_since_last_call_to_operate)
{
  auto& length_descriptor = m_field_info_ptr->m_length_descriptor;
  //Should be called for a pointer to contiguous segment - next to last dimension
  assert(idx_vector.size() == length_descriptor.get_num_dimensions()-1u);
  auto sep = length_descriptor.get_vcf_delimiter(length_descriptor.get_num_dimensions()-1u);
  //Print outer dimension separator
  if(!m_first_call)
  {
    assert(static_cast<size_t>(outermost_dim_idx_changed_since_last_call_to_operate)
        < length_descriptor.get_num_dimensions());
    (*m_fptr) << length_descriptor.get_vcf_delimiter(outermost_dim_idx_changed_since_last_call_to_operate);
  }
  auto num_elements = size_of_data_vec[0]
    /(m_field_info_ptr->get_genomicsdb_type().get_tuple_element_size(0u));
#ifdef DEBUG
  for(auto i=0u;i<m_field_info_ptr->get_genomicsdb_type().get_num_elements_in_tuple();++i)
    assert(num_elements == (size_of_data_vec[i]
          /(m_field_info_ptr->get_genomicsdb_type().get_tuple_element_size(i))));
#endif
  auto first_element = true;
  for(auto j=0ull;j<num_elements;++j)
  {
    for(auto i=0u;i<m_field_info_ptr->get_genomicsdb_type().get_num_elements_in_tuple();++i)
    {
      if(!first_element)
        (*m_fptr) << sep;
      switch(m_field_info_ptr->get_genomicsdb_type().get_tuple_element_bcf_ht_type(i))
      {
#define CASE_STATEMENTS(T) \
        { \
          cast_join_and_print<T>(*m_fptr, ptr_vec[i], j, 1u, sep); \
          break; \
        }
        case BCF_HT_INT:
          CASE_STATEMENTS(int);
        case BCF_HT_UINT:
          CASE_STATEMENTS(unsigned);
        case BCF_HT_INT64:
          CASE_STATEMENTS(int64_t);
        case BCF_HT_UINT64:
          CASE_STATEMENTS(uint64_t);
        case BCF_HT_REAL:
          CASE_STATEMENTS(float);
        case BCF_HT_DOUBLE:
          CASE_STATEMENTS(double);
        default:
          throw GenomicsDBMultiDVectorFieldOperatorException(std::string("Unhandled type in GenomicsDBMultiDVectorFieldVCFPrinter ")
              +m_field_info_ptr->get_genomicsdb_type().get_tuple_element_type_index(0u).name());
          break;
#undef CASE_STATEMENTS
      }
      first_element = false;
    }
  }
  m_first_call = false;
}
