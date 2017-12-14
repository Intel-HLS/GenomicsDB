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
    memcpy(array, str+element_begin_idx, current_element_length);
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
    memcpy(array, str+element_begin_idx, current_element_length);
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
void cast_join_and_print(std::ostream& fptr, const uint8_t* ptr, const size_t num_elements, const char sep)
{
  if(num_elements)
  {
    auto data_ptr = reinterpret_cast<const T*>(ptr);
    auto val = data_ptr[0];
    if(!is_bcf_missing_value<T>(val))
      fptr << val;
    for(auto i=1ull;i<num_elements;++i)
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
  assert(static_cast<size_t>(m_current_dim_idx) < m_field_info_ptr->m_length_descriptor.get_num_dimensions());
  assert(m_current_index_in_current_dimension < m_num_entries_in_current_dim);
  //the innermost dimension is simply a vector of elements - the size, #elements and offsets are NOT stored on disk,
  //Since the innermost dimension is simply a raw vector, hence the +2u check for dimensions which store 
  //size, #elements, offsets on disk
  if(m_current_dim_idx+1u < m_field_info_ptr->m_length_descriptor.get_num_dimensions())
  {
    auto update_offset = m_offsets_ptr[m_current_index_in_current_dimension+1u]
      - m_offsets_ptr[m_current_index_in_current_dimension]; 
    m_ro_field_ptr = m_ro_field_ptr + update_offset;
  }
  else
    m_ro_field_ptr = m_ro_field_ptr + (m_field_info_ptr->get_element_size());
  ++m_current_index_in_current_dimension;
}

/*
   Given a delimited string, parse and store the data in a binary buffer as described in
   the header.
   1,2,3,4|5,6|7,8$9,10|11|12
*/
template<class ElementType>
uint64_t GenomicsDBMultiDVectorField::parse_and_store_numeric(std::vector<uint8_t>& buffer,
        const FieldInfo& field_info,
        const char* str, const size_t str_length)
{
  auto& length_descriptor = field_info.m_length_descriptor;
  buffer.resize(4096u); //4KiB
  auto w_idx = 0ull; //write idx
  auto r_idx = 0ull; //read idx
  //#bytes for curr data in dim i
  std::vector<uint64_t> dim_sizes;
  dim_sizes.resize(length_descriptor.get_num_dimensions(), 0ull);
  //Offsets for the data in each dim F[0] offset, F[1] offset etc
  std::vector<std::vector<uint64_t>> dim_offsets;
  //Each dimension begins with offset 0
  dim_offsets.resize(length_descriptor.get_num_dimensions(), std::vector<uint64_t>(1u, 0ull));
  //Specifies the current offset in buffer at which the data for dim i should be written
  //the last 2 dimensions don't get any sizes
  std::vector<uint64_t> dim_write_begin_offsets;
  dim_write_begin_offsets.resize(length_descriptor.get_num_dimensions());
  //The first dimension begins writing at offset 0
  //A 64-bit uint64_t is allocated for each dimension's size
  for(auto i=0ull;i<dim_write_begin_offsets.size();++i)
    dim_write_begin_offsets[i] = i*sizeof(uint64_t);
  //Points to the begin of the current element - length
  auto current_element_begin_read_idx = 0ull; //r stands for read
  auto current_element_length = 0ull;
  auto num_elements_in_innermost_dim_read = 0ull;
  auto total_size_of_multi_d_data = 0ull;
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
      auto value = str_to_element<ElementType>(str, current_element_begin_read_idx, current_element_length);
      int64_t write_offset = dim_write_begin_offsets.back();
      File2TileDBBinaryBase::tiledb_buffer_resize_if_needed_and_print<ElementType>(
          buffer, write_offset, value);
      dim_write_begin_offsets.back() = write_offset;
      ++num_elements_in_innermost_dim_read;
      //Delimiter not of the innermost dimension
      if(sep_dim_idx+1u < length_descriptor.get_num_dimensions())
      {
        //Single element - but actually empty string
        if(num_elements_in_innermost_dim_read == 1u
            && current_element_length == 0u)
          num_elements_in_innermost_dim_read = 0u;
        //Example: 1,2,3,4|5,6|7,8$9,10|11|12
        //Length of innermost dimension is the number of elements read
        auto last_dim_size =
          (num_elements_in_innermost_dim_read*sizeof(ElementType));
        //Reset number of elements in innermost dim to 0
        num_elements_in_innermost_dim_read = 0u;
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
            memcpy(&(buffer[offsets_start_offset]), &(dim_offsets[i][0u]), offsets_size);
            //Update last_dim_size
            last_dim_size = (dim_sizes[i]
                + sizeof(uint64_t)    //for storing size of data at dim i
                + sizeof(uint64_t)    //for storing #entries at dim i
                + offsets_size        //for storing offsets at dim i
                );
          }
        }
        if(sep_dim_idx == -1)
          total_size_of_multi_d_data = last_dim_size;
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
      current_element_begin_read_idx = r_idx+1u;
      current_element_length = 0u;
    }
    else
      ++current_element_length;
    ++r_idx;
  }
  return total_size_of_multi_d_data;
}

template<class ElementType>
uint64_t GenomicsDBMultiDVectorField::parse_and_store_numeric(const char* str, const size_t str_length)
{
  auto total_size_of_multi_d_data = GenomicsDBMultiDVectorField::parse_and_store_numeric<ElementType>(m_rw_field_data, *m_field_info_ptr,
      str, str_length);
  m_rw_field_data.resize(total_size_of_multi_d_data);
  return total_size_of_multi_d_data;
}


//Template instantiations
template
uint64_t GenomicsDBMultiDVectorField::parse_and_store_numeric<int>(const char* str, const size_t str_length);
template
uint64_t GenomicsDBMultiDVectorField::parse_and_store_numeric<int64_t>(const char* str, const size_t str_length);
template
uint64_t GenomicsDBMultiDVectorField::parse_and_store_numeric<float>(const char* str, const size_t str_length);

//Iterating over all dimensions and all index values can be a recursive process
//We can model this as a stack
void GenomicsDBMultiDVectorField::run_operation(GenomicsDBMultiDVectorFieldOperator& multid_vector_field_operator,
    const uint8_t* data_ptr) const
{
  auto& length_descriptor = m_field_info_ptr->m_length_descriptor;
  //Stack replacing the recursive function
  std::vector<GenomicsDBMultiDVectorIdx> idx_stack;
  idx_stack.emplace_back(data_ptr, m_field_info_ptr, 0u);
  //current index vector e.g. A[5][0][3] will have 5,0,3
  //changes as stack is traversed
  //Don't bother with the last dimension
  std::vector<uint64_t> curr_index_vector(length_descriptor.get_num_dimensions()-1u);
  //Optimization - tracks the outermost dimension index which changed since the
  //last call to operate(). Initialized to N-2
  auto outermost_dim_idx_changed_since_last_call_to_operate = length_descriptor.get_num_dimensions()-2u;
  while(!idx_stack.empty())
  {
    auto& top_idx = idx_stack.back();
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
        auto copy_idx = top_idx; //copy
        //index 0 in the next dimension
        copy_idx.advance_to_index_in_next_dimension(0u);
        //Cannot move this to outside the if-else block because idx_stack memory can be reallocated
        //by the emplace_back() statement
        top_idx.advance_index_in_current_dimension();
        idx_stack.emplace_back(copy_idx);
      }
      else
      {
        multid_vector_field_operator.operate(top_idx.get_ptr<uint8_t>(),
            top_idx.get_size_of_current_index(), curr_index_vector,
            outermost_dim_idx_changed_since_last_call_to_operate);
        //reset outermost_dim_idx_changed_since_last_call_to_operate 
        outermost_dim_idx_changed_since_last_call_to_operate = length_descriptor.get_num_dimensions()-2u; 
        //Cannot move this to outside the if-else block because idx_stack memory can be reallocated (if block)
        top_idx.advance_index_in_current_dimension();
      }
      //DO NOT USE top_idx in this if block after this - emplace_back() might have reallocated the vector
    }
    else
    {
      idx_stack.pop_back();
      --outermost_dim_idx_changed_since_last_call_to_operate;
    }
  }
}

//GenomicsDBMultiDVectorFieldOperator functions

GenomicsDBMultiDVectorFieldVCFPrinter::GenomicsDBMultiDVectorFieldVCFPrinter(std::ostream& fptr, const FieldInfo& field_info)
{
  m_first_call = true;
  m_type_enum = VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(field_info.get_genomicsdb_type_index());
  m_fptr = &fptr;
  m_field_info_ptr = &field_info;
}

void GenomicsDBMultiDVectorFieldVCFPrinter::operate(const uint8_t* ptr, const size_t size_of_data,
        const std::vector<uint64_t>& idx_vector, int outermost_dim_idx_changed_since_last_call_to_operate)
{
  auto& length_descriptor = m_field_info_ptr->m_length_descriptor;
  //Should be called for a pointer to contiguous segment - next to last dimension
  assert(idx_vector.size() == length_descriptor.get_num_dimensions()-1u);
  auto num_elements = size_of_data/m_field_info_ptr->get_element_size();
  auto sep = length_descriptor.get_vcf_delimiter(length_descriptor.get_num_dimensions()-1u);
  //Print outer dimension separator
  if(!m_first_call)
  {
    assert(static_cast<size_t>(outermost_dim_idx_changed_since_last_call_to_operate)
        < length_descriptor.get_num_dimensions());
    (*m_fptr) << length_descriptor.get_vcf_delimiter(outermost_dim_idx_changed_since_last_call_to_operate);
  }
  switch(m_type_enum)
  {
    case VariantFieldTypeEnum::VARIANT_FIELD_INT:
      cast_join_and_print<int>(*m_fptr, ptr, num_elements, sep);
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_FLOAT:
      cast_join_and_print<float>(*m_fptr, ptr, num_elements, sep);
      break;
    default:
      throw GenomicsDBMultiDVectorFieldOperatorException(std::string("Unhandled type in GenomicsDBMultiDVectorFieldVCFPrinter ")
          +m_field_info_ptr->get_genomicsdb_type_index().name());
      break;
  }
  m_first_call = false;
}
