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

#include "genomicsdb_columnar_field.h"
#include "variant_field_data.h"

GenomicsDBColumnarField::GenomicsDBColumnarField(const std::type_index element_type, const int length_descriptor,
    const unsigned fixed_length_field_length, const size_t num_bytes)
: m_element_type(element_type), m_length_descriptor(length_descriptor),
  m_fixed_length_field_num_elements(fixed_length_field_length)
{
  m_element_size = VariantFieldTypeUtil::size(m_element_type);
  m_fixed_length_field_size = m_fixed_length_field_num_elements*m_element_size;
  m_buffer_size = (length_descriptor == BCF_VL_FIXED)
    ? GET_ALIGNED_BUFFER_SIZE(num_bytes, m_fixed_length_field_size)
    : GET_ALIGNED_BUFFER_SIZE(num_bytes, m_element_size);
  assign_function_pointers();
  m_free_buffer_list_head_ptr = 0;
  m_live_buffer_list_head_ptr = 0;
  m_live_buffer_list_tail_ptr = 0;
  m_curr_index_in_live_buffer_list_tail = 0u;
  //Create 1 buffer
  add_new_buffer();
}

void GenomicsDBColumnarField::copy_simple_members(const GenomicsDBColumnarField& other)
{
  m_length_descriptor = other.m_length_descriptor;
  m_fixed_length_field_num_elements = other.m_fixed_length_field_num_elements;
  m_fixed_length_field_size = other.m_fixed_length_field_size;
  m_element_size = other.m_element_size;
  m_element_type = other.m_element_type;
  m_check_tiledb_valid_element = other.m_check_tiledb_valid_element;
  m_buffer_size = other.m_buffer_size;
  m_curr_index_in_live_buffer_list_tail = other.m_curr_index_in_live_buffer_list_tail;
}

GenomicsDBColumnarField::GenomicsDBColumnarField(GenomicsDBColumnarField&& other)
  : m_element_type(typeid(void))
{
  copy_simple_members(other);
  m_free_buffer_list_head_ptr = other.m_free_buffer_list_head_ptr;
  other.m_free_buffer_list_head_ptr = 0;
  m_live_buffer_list_head_ptr = other.m_live_buffer_list_head_ptr;
  other.m_live_buffer_list_head_ptr = 0;
  m_live_buffer_list_tail_ptr = other.m_live_buffer_list_tail_ptr;
  other.m_live_buffer_list_tail_ptr = 0;
}

GenomicsDBColumnarField::~GenomicsDBColumnarField()
{
  GenomicsDBBuffer* ptr = m_free_buffer_list_head_ptr;
  GenomicsDBBuffer* next_ptr = ptr;
  while(ptr)
  {
    next_ptr = ptr->get_next_buffer();
    delete ptr;
    ptr = next_ptr;
  }
  m_free_buffer_list_head_ptr = 0;
  ptr = m_live_buffer_list_head_ptr;
  while(ptr)
  {
    next_ptr = ptr->get_next_buffer();
    delete ptr;
    ptr = next_ptr;
  }
  m_live_buffer_list_head_ptr = 0;
  m_live_buffer_list_tail_ptr = 0;
  clear();
}

void GenomicsDBColumnarField::assign_function_pointers()
{
  switch(VariantFieldTypeUtil::get_variant_field_type_enum_for_variant_field_type(m_element_type))
  {
    case VariantFieldTypeEnum::VARIANT_FIELD_INT:
    case VariantFieldTypeEnum::VARIANT_FIELD_UNSIGNED:
      m_check_tiledb_valid_element = GenomicsDBColumnarField::check_tiledb_valid_element<int>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_INT64_T:
    case VariantFieldTypeEnum::VARIANT_FIELD_UINT64_T:
      m_check_tiledb_valid_element = GenomicsDBColumnarField::check_tiledb_valid_element<int64_t>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_FLOAT:
      m_check_tiledb_valid_element = GenomicsDBColumnarField::check_tiledb_valid_element<float>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_DOUBLE:
      m_check_tiledb_valid_element = GenomicsDBColumnarField::check_tiledb_valid_element<double>;
      break;
    case VariantFieldTypeEnum::VARIANT_FIELD_CHAR:
    case VariantFieldTypeEnum::VARIANT_FIELD_STRING:
      m_check_tiledb_valid_element = GenomicsDBColumnarField::check_tiledb_valid_element<char>;
      break;
    default:
      throw GenomicsDBColumnarFieldException(std::string("Unhandled type ")+m_element_type.name());
  }
}

void GenomicsDBColumnarField::move_buffer_to_live_list(GenomicsDBBuffer* buffer)
{
  assert(buffer && !(buffer->is_in_live_list()));
  //Store neighbours in free list in tmp variables
  auto next_in_free_list = buffer->get_next_buffer();
  auto previous_in_free_list = buffer->get_previous_buffer();
  //Insert into live list
  //Since this will be the last element in the live list
  buffer->set_next_buffer(0);
  if(m_live_buffer_list_head_ptr == 0)
  {
    assert(m_live_buffer_list_tail_ptr == 0);
    m_live_buffer_list_head_ptr = buffer;
    m_live_buffer_list_tail_ptr = buffer;
    buffer->set_previous_buffer(0);
  }
  else
  {
    assert(m_live_buffer_list_tail_ptr);
    assert(m_live_buffer_list_tail_ptr->get_next_buffer() == 0);
    m_live_buffer_list_tail_ptr->set_next_buffer(buffer);
    buffer->set_previous_buffer(m_live_buffer_list_tail_ptr);
    m_live_buffer_list_tail_ptr = buffer;
  }
  buffer->set_is_in_live_list(true);
  //Fix neighbours of buffer in free list
  if(next_in_free_list)
    next_in_free_list->set_previous_buffer(previous_in_free_list);
  if(previous_in_free_list)
    previous_in_free_list->set_next_buffer(next_in_free_list);
  //If this was the head pointer, advance it
  if(buffer == m_free_buffer_list_head_ptr)
    m_free_buffer_list_head_ptr = next_in_free_list;
  //Reset idx pointed to in tail to 0
  m_curr_index_in_live_buffer_list_tail = 0u;
}

void GenomicsDBColumnarField::move_buffer_to_free_list(GenomicsDBBuffer* buffer)
{
  assert(buffer && (buffer->is_in_live_list()));
  //Store neighbours in live list in tmp variables
  auto next_in_live_list = buffer->get_next_buffer();
  auto previous_in_live_list = buffer->get_previous_buffer();
  //Insert into free list
  //Since this will be the first element in the free list
  buffer->set_previous_buffer(0);
  buffer->set_next_buffer(m_free_buffer_list_head_ptr);
  if(m_free_buffer_list_head_ptr)
    m_free_buffer_list_head_ptr->set_previous_buffer(buffer);
  m_free_buffer_list_head_ptr = buffer;
  buffer->set_is_in_live_list(false);
  //Fix neighbours of buffer in live list
  if(next_in_live_list)
    next_in_live_list->set_previous_buffer(previous_in_live_list);
  if(previous_in_live_list)
    previous_in_live_list->set_next_buffer(next_in_live_list);
  //If this was the tail, update to previous
  if(buffer == m_live_buffer_list_tail_ptr)
    m_live_buffer_list_tail_ptr = previous_in_live_list;
  //If this was the head, update to next
  if(buffer == m_live_buffer_list_head_ptr)
    m_live_buffer_list_head_ptr = next_in_live_list;
}

void GenomicsDBColumnarField::set_valid_vector_in_live_buffer_list_tail_ptr()
{
  auto buffer_ptr = get_live_buffer_list_tail_ptr();
  assert(buffer_ptr);
  auto& valid_vector = buffer_ptr->get_valid_vector();
  if(m_length_descriptor == BCF_VL_FIXED)
    for(auto i=0ull;i<buffer_ptr->get_num_live_entries();++i)
    {
      assert(i < valid_vector.size());
      valid_vector[i] = m_check_tiledb_valid_element(buffer_ptr->get_buffer_pointer() + (m_fixed_length_field_size*i),
          m_fixed_length_field_num_elements);
    }
  else
    for(auto i=0ull;i<buffer_ptr->get_num_live_entries();++i)
    {
      assert(i < valid_vector.size());
      valid_vector[i] = (buffer_ptr->get_size_of_variable_length_field(i) > 0u);
    }
}
