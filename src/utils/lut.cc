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

#include "lut.h"

using namespace std;

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::LUTBase()
{
  m_num_input_vcfs = 0u;
  m_num_merged_fields = 0u;
  clear();
}

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::LUTBase(int64_t numInputGVCFs, int64_t numMergedFields)
{
  m_num_input_vcfs = numInputGVCFs;
  m_num_merged_fields = numMergedFields;
  clear();
  resize_inputs_2_merged_lut_if_needed(numInputGVCFs, numMergedFields);
  resize_merged_2_inputs_lut_if_needed(numInputGVCFs, numMergedFields);
}

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
void LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::clear()
{
  for(auto& vec : m_inputs_2_merged_lut)
    vec.clear();
  m_inputs_2_merged_lut.clear();
  for(auto& vec : m_merged_2_inputs_lut)
    vec.clear();
  m_merged_2_inputs_lut.clear();
}

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
void LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::reset_vector(vector<int64_t>& vec, int64_t from)
{
  for(auto i=from;i<static_cast<int64_t>(vec.size());++i)
    vec[i] = lut_missing_value;
}

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
void LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_and_reset_vector(vector<int64_t>& vec, int64_t new_size)
{
  auto old_size = vec.size();
  if(new_size > static_cast<int64_t>(old_size))
  {
    vec.resize(new_size);
    reset_vector(vec, old_size);
  }
}

template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
void LUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_and_reset_lut
(vector<vector<int64_t>>& lut, int64_t new_lut_size, int64_t new_vector_size, int64_t& numRowsVar, int64_t& numColsVar)
{
  auto old_lut_size = lut.size();
  if(new_lut_size > static_cast<int64_t>(old_lut_size))
  {
    lut.resize(new_lut_size);
    numRowsVar = new_lut_size;
  }
  auto old_vector_size = (lut.size() > 0u) ? lut[0].size() : 0u;
  //Begin resizing of vectors at start_idx
  auto start_idx = old_lut_size;
  if(new_vector_size > static_cast<int64_t>(old_vector_size))     //every vector needs to be resized
  {
    start_idx = 0u;
    numColsVar = new_vector_size;
  }
  else
    new_vector_size = old_vector_size;      //new vector size is smaller, don't bother reducing the size of existing rows
  for(auto i=start_idx;static_cast<int64_t>(i)<new_lut_size;++i)
    resize_and_reset_vector(lut[i], new_vector_size);
}
//explicit initialization to avoid link errors
template class LUTBase<true,true>;
template class LUTBase<true,false>;
template class LUTBase<false,true>;
template class LUTBase<false,false>;

//explicit initialization to avoid link errors
template class MergedAllelesIdxLUT<true,true>;

//explicit initialization to avoid link errors
template class GoldLUTTemplate<true, true>;

