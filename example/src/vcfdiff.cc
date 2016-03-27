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

#include <getopt.h>
#include <mpi.h>
#include "vcfdiff.h"
#include "vid_mapper.h"
#include "json_config.h"

double g_threshold = 1e-5; //fp comparison threshold
int64_t g_num_callsets = INT64_MAX;

#define VERIFY_OR_THROW(X) if(!(X)) throw VCFDiffException(#X);

VCFDiffFile::VCFDiffFile(const std::string& filename, const std::string& regions)
: m_samples_lut(1u, 1u), m_fields_lut(1u, 1u), m_contigs_lut(1u, 1u)
{
  m_filename = filename;
  m_regions = "";
  m_reader = bcf_sr_init();
  m_num_lines_read = 0ll;
  m_hdr = 0;
  m_line = 0;
  //Tmp buffer
  m_tmp_hts_string.m = 65536;    //64KB
  m_tmp_hts_string.s = (char*)malloc(m_tmp_hts_string.m*sizeof(char));
  VERIFY_OR_THROW(m_tmp_hts_string.s);
  m_tmp_hts_string.l = 0;
  //Read hdr
  auto fptr = bcf_open(filename.c_str(), "r");
  VERIFY_OR_THROW(fptr);
  m_hdr = bcf_hdr_read(fptr);
  bcf_close(fptr);
  //Contigs
  for(auto i=0;i<m_hdr->n[BCF_DT_CTG];++i)
    if(m_contigs.find(bcf_hdr_id2name(m_hdr, i)) == m_contigs.end())
        m_contigs.insert(bcf_hdr_id2name(m_hdr, i));
  m_contigs_lut.resize_luts_if_needed(1u, m_contigs.size());
  //Fields
  m_field_idx_to_known_field_enum.resize_luts_if_needed(m_hdr->n[BCF_DT_ID], GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0;i<m_hdr->n[BCF_DT_ID];++i)
  {
    if(m_fields.find(bcf_hdr_int2id(m_hdr, BCF_DT_ID, i)) == m_fields.end())
    {
      std::string field_name = std::move(std::string(bcf_hdr_int2id(m_hdr, BCF_DT_ID, i)));
      m_fields.insert(field_name);
      unsigned known_field_enum = 0;
      auto status = KnownFieldInfo::get_known_field_enum_for_name(field_name, known_field_enum);
      if(status)
        m_field_idx_to_known_field_enum.add_schema_idx_known_field_mapping(i, known_field_enum);
    }
  }
  m_fields_lut.resize_luts_if_needed(1u, m_fields.size());
  //Samples
  for(auto i=0;i<bcf_hdr_nsamples(m_hdr);++i)
    if(m_samples.find(bcf_hdr_int2id(m_hdr, BCF_DT_SAMPLE, i)) == m_samples.end())
      m_samples.insert(bcf_hdr_int2id(m_hdr, BCF_DT_SAMPLE, i));
  m_samples_lut.resize_luts_if_needed(1u, m_samples.size());
}

VCFDiffFile::~VCFDiffFile()
{
  m_filename.clear();
  m_regions.clear();
  m_line = 0;   //line is managed by reader
  if(m_reader)
    bcf_sr_destroy(m_reader);
  m_reader = 0;
  if(m_hdr)
    bcf_hdr_destroy(m_hdr);
  m_hdr = 0;
  m_contigs.clear();
  m_fields.clear();
  m_samples.clear();
  m_contigs_lut.clear();
  m_fields_lut.clear();
  m_samples_lut.clear();
  m_fields_in_gold_line.clear();
  m_fields_in_test_line.clear();
  m_gold_genotype_idx_to_test_idx.clear();
  if(m_tmp_hts_string.m && m_tmp_hts_string.s)
    free(m_tmp_hts_string.s);
  m_tmp_hts_string.s = 0;
  m_tmp_hts_string.m = 0;
}

void VCFDiffFile::setup_lut(const std::set<std::string>& gold_set, const std::set<std::string>& test_set,
    const int bcf_dt_type, GoldLUT& lut, const VCFDiffFile& gold)
{
  lut.resize_luts_if_needed(1u, std::max(gold_set.size(), test_set.size()));
  std::vector<std::string> common_elements(std::min(gold_set.size(), test_set.size()));
  //Find common fields
  auto iter = std::set_intersection(gold_set.begin(), gold_set.end(), test_set.begin(), test_set.end(),
      common_elements.begin());
  common_elements.resize(iter - common_elements.begin());
  for(const auto& x : common_elements)
  {
    auto idx = bcf_hdr_id2int(m_hdr, bcf_dt_type, x.c_str());
    VERIFY_OR_THROW(idx >= 0);
    auto gold_idx = bcf_hdr_id2int(gold.m_hdr, bcf_dt_type, x.c_str());
    VERIFY_OR_THROW(gold_idx >= 0);
    lut.add_input_merged_idx_pair(0, idx, gold_idx);
  }
}

void VCFDiffFile::setup_luts(const VCFDiffFile& gold, const bool use_callsets_file_for_samples)
{
  if(!use_callsets_file_for_samples && gold.m_samples != m_samples)
    throw VCFDiffException("ERROR: Sample names in the 2 file headers are different - cannot perform any further comparisons");
  setup_lut(gold.m_contigs, m_contigs, BCF_DT_CTG, m_contigs_lut, gold);
  setup_lut(gold.m_fields, m_fields, BCF_DT_ID, m_fields_lut, gold);
  if(!use_callsets_file_for_samples)
    setup_lut(gold.m_samples, m_samples, BCF_DT_SAMPLE, m_samples_lut, gold);
  m_fields_in_gold_line.resize(std::max(gold.m_fields.size(), m_fields.size()));
  m_fields_in_test_line.resize(std::max(gold.m_fields.size(), m_fields.size()));
  reset_field_to_line_idx_mapping();
}

void VCFDiffFile::set_regions_and_open_file()
{
  bcf_sr_set_regions(m_reader, m_regions.c_str(), 0);
  if(bcf_sr_add_reader(m_reader, m_filename.c_str()) != 1)
    throw VCFDiffException(std::string("Could not open file ")+m_filename+" or its index - VCF/BCF files must be block compressed and indexed");
  bcf_sr_next_line(m_reader);
  m_line = bcf_sr_get_line(m_reader, 0);
  if(m_line)
    ++m_num_lines_read;
}

void VCFDiffFile::read_and_advance()
{
  assert(m_line);
  bcf_sr_next_line(m_reader);
  m_line = bcf_sr_get_line(m_reader, 0);
  if(m_line)
    ++m_num_lines_read;
}

void VCFDiffFile::seek_and_read(const int rid, const int pos)
{
  assert(rid >= 0 && rid < m_hdr->n[BCF_DT_CTG]);
  auto contig = bcf_hdr_id2name(m_hdr, rid);
  bcf_sr_seek(m_reader, contig, pos);
  bcf_sr_next_line(m_reader);
  m_line = bcf_sr_get_line(m_reader, 0);
  if(m_line)
    ++m_num_lines_read;
}

void VCFDiffFile::print_line(std::ostream& fptr)
{
  m_tmp_hts_string.l = 0;
  vcf_format(m_hdr, m_line, &m_tmp_hts_string);
  fptr << m_tmp_hts_string.s << "\n";
}

inline int GET_NUM_FIELDS(const bcf1_t* line, const int bcf_field_type)
{
  if(bcf_field_type == BCF_HL_INFO)
    return line->n_info;
  else
  {
    assert(bcf_field_type == BCF_HL_FMT);
    return line->n_fmt;
  }
}

inline int GET_FIELD_IDX(const bcf1_t* line, const int bcf_field_type, const int idx)
{
  assert(idx >= 0);
  if(bcf_field_type == BCF_HL_INFO)
  {
    assert(idx < line->n_info);
    return line->d.info[idx].key;
  }
  else
  {
    assert(bcf_field_type == BCF_HL_FMT);
    assert(idx < line->n_fmt);
    return line->d.fmt[idx].id;
  }
}

inline int GET_FIELD_BT_TYPE(const bcf1_t* line, const int bcf_field_type, const int idx)
{
  assert(idx >= 0);
  if(bcf_field_type == BCF_HL_INFO)
  {
    assert(idx < line->n_info);
    return line->d.info[idx].type;
  }
  else
  {
    assert(bcf_field_type == BCF_HL_FMT);
    assert(idx < line->n_fmt);
    return line->d.fmt[idx].type;
  }
}

inline int GET_NUM_ELEMENTS(const bcf1_t* line, const int bcf_field_type, const int idx)
{
  assert(idx >= 0);
  if(bcf_field_type == BCF_HL_INFO)
  {
    assert(idx < line->n_info);
    return line->d.info[idx].len;
  }
  else
  {
    assert(bcf_field_type == BCF_HL_FMT);
    assert(idx < line->n_fmt);
    return line->d.fmt[idx].n;
  }
}

template<class T>
inline T GET_DATA_PTR(const bcf1_t* line, const int bcf_field_type, const int idx, const int sample_idx, const int num_elements_per_sample)
{
  assert(idx >= 0);
  if(bcf_field_type == BCF_HL_INFO)
  {
    assert(idx < line->n_info);
    return reinterpret_cast<T>(line->d.info[idx].vptr);
  }
  else
  {
    assert(bcf_field_type == BCF_HL_FMT);
    assert(idx < line->n_fmt);
    return reinterpret_cast<T>(line->d.fmt[idx].p) + sample_idx*num_elements_per_sample;
  }
}

bool bcf_is_empty_field(const bcf_hdr_t* hdr, const bcf1_t* line, const int bcf_field_type, const int idx)
{
  assert(idx >= 0);
  auto num_elements_per_sample = GET_NUM_ELEMENTS(line, bcf_field_type, idx);
  auto total_num_elements = (bcf_field_type == BCF_HL_FMT) ? num_elements_per_sample*bcf_hdr_nsamples(hdr) : num_elements_per_sample;
  switch(GET_FIELD_BT_TYPE(line, bcf_field_type, idx))
  {
    case BCF_BT_INT8:
      {
        auto ptr = GET_DATA_PTR<const int8_t*>(line, bcf_field_type, idx, 0, num_elements_per_sample);
        for(auto i=0;i<total_num_elements;++i)
          if(is_bcf_valid_value<int8_t>(ptr[i]))
            return false;
        break;
      }
    case BCF_BT_INT16:
      {
        auto ptr = GET_DATA_PTR<const int16_t*>(line, bcf_field_type, idx, 0, num_elements_per_sample);
        for(auto i=0;i<total_num_elements;++i)
          if(is_bcf_valid_value<int16_t>(ptr[i]))
            return false;
        break;
      }
    case BCF_BT_INT32:
      {
        auto ptr = GET_DATA_PTR<const int32_t*>(line, bcf_field_type, idx, 0, num_elements_per_sample);
        for(auto i=0;i<total_num_elements;++i)
          if(is_bcf_valid_value<int32_t>(ptr[i]))
            return false;
        break;
      }
    case BCF_BT_FLOAT:
      {
        auto ptr = GET_DATA_PTR<const float*>(line, bcf_field_type, idx, 0, num_elements_per_sample);
        for(auto i=0;i<total_num_elements;++i)
          if(is_bcf_valid_value<float>(ptr[i]))
            return false;
        break;
      }
    case BCF_BT_CHAR:
      {
        auto ptr = GET_DATA_PTR<const char*>(line, bcf_field_type, idx, 0, num_elements_per_sample);
        for(auto i=0;i<total_num_elements;++i)
          if(is_bcf_valid_value<char>(ptr[i]))
            return false;
        break;
      }
    default:
      throw VCFDiffException("Unknown BCF_BT_TYPE "+std::to_string(GET_FIELD_BT_TYPE(line, bcf_field_type, idx)));
      break;
  }
  return true;
}

template<class T1, class T2>
inline bool compare_unequal(const T1 a, const T2 b)
{
  return !(
      (a == b) ||
      (is_bcf_missing_value<T1>(a) && is_bcf_missing_value<T2>(b)) ||
      (is_bcf_vector_end_value<T1>(a) && is_bcf_vector_end_value<T2>(b))
      );
}

//Specialization for float - values should be close
template<>
inline bool compare_unequal(const float a, const float b)
{
  if((is_bcf_missing_value<float>(a) && is_bcf_missing_value<float>(b)) ||
      (is_bcf_vector_end_value<float>(a) && is_bcf_vector_end_value<float>(b)))
    return false;
  float diff = fabsf(a-b);
  float rel_diff = ((a != 0) ? fabsf(diff/a) : 0);
  return (diff > g_threshold && rel_diff > g_threshold);
}

template<class T1, class T2>
bool VCFDiffFile::compare_unequal_vector(const bcf_hdr_t* gold_hdr, const bcf1_t* gold_line, const int bcf_field_type,
    int gold_line_field_pos_idx, int test_line_field_pos_idx)
{
  auto num_samples = (bcf_field_type == BCF_HL_INFO) ? 1 : bcf_hdr_nsamples(gold_hdr);
  auto num_gold_elements = GET_NUM_ELEMENTS(gold_line, bcf_field_type, gold_line_field_pos_idx);
  auto num_test_elements = GET_NUM_ELEMENTS(m_line, bcf_field_type, test_line_field_pos_idx);
  auto min_num_per_sample = std::min(num_gold_elements, num_test_elements);
  auto max_num_per_sample = std::max(num_gold_elements, num_test_elements);
  //Field idxs
  auto gold_field_idx = GET_FIELD_IDX(gold_line, bcf_field_type, gold_line_field_pos_idx);
  auto test_field_idx = GET_FIELD_IDX(m_line, bcf_field_type, test_line_field_pos_idx);
  //Length descriptor
  auto known_field_enum = m_field_idx_to_known_field_enum.get_known_field_enum_for_schema_idx(test_field_idx);
  auto length_descriptor = BCF_VL_FIXED;
  if(m_field_idx_to_known_field_enum.is_defined_value(known_field_enum))
    length_descriptor = KnownFieldInfo::get_length_descriptor_for_known_field_enum(known_field_enum);
  else
  {
    auto gold_length_descriptor = bcf_hdr_id2length(gold_hdr, bcf_field_type, gold_field_idx);
    auto test_length_descriptor = bcf_hdr_id2length(m_hdr, bcf_field_type, test_field_idx);
    if(gold_length_descriptor != test_length_descriptor)
      return true;
    length_descriptor = gold_length_descriptor;
  }
  //Vector lengths must be identical for these types of fields
  //Also, if the alleles don't match, don't bother checking further
  if((length_descriptor == BCF_VL_G || length_descriptor == BCF_VL_R || length_descriptor == BCF_VL_A)
    && (num_test_elements != num_gold_elements || m_diff_alleles_flag))
    return true;
  for(auto j=0;j<num_samples;++j)
  {
    int lut_test_sample_idx = (bcf_field_type == BCF_HL_INFO) ? 0 : m_samples_lut.get_test_idx_for_gold(0, j);
    assert(!GoldLUT::is_missing_value(lut_test_sample_idx) && lut_test_sample_idx < bcf_hdr_nsamples(m_hdr));
    auto gold_ptr = GET_DATA_PTR<const T1*>(gold_line, bcf_field_type, gold_line_field_pos_idx, j, num_gold_elements);
    auto test_ptr = GET_DATA_PTR<const T2*>(m_line, bcf_field_type, test_line_field_pos_idx, lut_test_sample_idx, num_test_elements);
    for(auto k=0;k<min_num_per_sample;++k)
    {
      auto test_vector_idx = k;
      switch(length_descriptor)
      {
        case BCF_VL_A:
          test_vector_idx = m_alleles_lut.get_input_idx_for_merged(0, k+1) - 1;   //ALT only
          break;
        case BCF_VL_R:
          test_vector_idx = m_alleles_lut.get_input_idx_for_merged(0, k);
          break;
        case BCF_VL_G:
          assert(static_cast<size_t>(k) < m_gold_genotype_idx_to_test_idx.size());
          test_vector_idx = m_gold_genotype_idx_to_test_idx[k];
          break;
        default:
          break;        //test_vector_idx = k
      }
      assert(!CombineAllelesLUT::is_missing_value(test_vector_idx) && test_vector_idx < min_num_per_sample);
      if(compare_unequal<T1, T2>(gold_ptr[k], test_ptr[test_vector_idx]))
        return true;    //unequal
    }
    if(num_test_elements > num_gold_elements)
    {
      for(auto k=min_num_per_sample;k<max_num_per_sample;++k)
        if(is_bcf_valid_value<T2>(test_ptr[k]))
          return true;  //extra elements
    }
    else
    {
      for(auto k=min_num_per_sample;k<max_num_per_sample;++k)
        if(is_bcf_valid_value<T1>(gold_ptr[k]))
          return true;  //extra elements
    }
  }
  return false;
}

bool VCFDiffFile::compare_unequal_fields(const bcf_hdr_t* gold_hdr, const bcf1_t* gold_line, const int bcf_field_type, std::string& error_message)
{
  reset_field_to_line_idx_mapping();
  auto diff_fields = false;
  for(auto i=0;i<GET_NUM_FIELDS(gold_line, bcf_field_type);++i)
  {
    auto field_idx = GET_FIELD_IDX(gold_line, bcf_field_type, i);
    assert(field_idx >= 0 && static_cast<size_t>(field_idx) < m_fields_in_gold_line.size());
    m_fields_in_gold_line[field_idx] = i;
  }
  for(auto i=0;i<GET_NUM_FIELDS(m_line, bcf_field_type);++i)
  {
    auto field_idx = GET_FIELD_IDX(m_line, bcf_field_type, i);
    assert(field_idx >= 0 && static_cast<size_t>(field_idx) < m_fields_in_test_line.size());
    m_fields_in_test_line[field_idx] = i;
  }
  for(auto i=0;i<GET_NUM_FIELDS(gold_line, bcf_field_type);++i)
  {
    auto gold_field_idx = GET_FIELD_IDX(gold_line, bcf_field_type, i);
    auto lut_test_field_idx = m_fields_lut.get_test_idx_for_gold(0, gold_field_idx);
    //Skip comparing END values as the big outer loop handles that
    if(!GoldLUT::is_missing_value(lut_test_field_idx))
    {
      auto known_field_enum = m_field_idx_to_known_field_enum.get_known_field_enum_for_schema_idx(lut_test_field_idx);
      if(m_field_idx_to_known_field_enum.is_defined_value(known_field_enum) && known_field_enum == GVCF_END_IDX)
        continue;
    }
    assert(GoldLUT::is_missing_value(lut_test_field_idx) || static_cast<size_t>(lut_test_field_idx) < m_fields_in_test_line.size());
    //Missing field in test
    if(GoldLUT::is_missing_value(lut_test_field_idx) || m_fields_in_test_line[lut_test_field_idx] < 0)
    {
      if(!bcf_is_empty_field(gold_hdr, gold_line, bcf_field_type, i))
      {
        diff_fields = true;
        error_message +=  (std::string("ERROR: ") + (bcf_field_type == BCF_HL_INFO ? "INFO" : "FORMAT")
           + " field " + bcf_hdr_int2id(gold_hdr, BCF_DT_ID, gold_field_idx) + " missing in test line\n");
      }
    }
    else
    {
      auto test_line_field_pos_idx = m_fields_in_test_line[lut_test_field_idx];
      //Fields are of different type
      if(bcf_hdr_id2type(gold_hdr, bcf_field_type, gold_field_idx) != bcf_hdr_id2type(m_hdr, bcf_field_type, lut_test_field_idx))
      {
        diff_fields = true;
        error_message += (std::string("ERROR: ") + (bcf_field_type == BCF_HL_INFO ? "INFO" : "FORMAT")
          + " field " + bcf_hdr_int2id(gold_hdr, BCF_DT_ID, gold_field_idx) + " type is different in gold and test\n");
        continue;
      }
      //Boolean field, both lines have the field - nothing to do
      if(bcf_hdr_id2type(gold_hdr, bcf_field_type, gold_field_idx) == BCF_HT_FLAG)
        continue;
      //Big switch
      auto field_bt_type_combo = GET_FIELD_BT_TYPE(gold_line, bcf_field_type, i)  << 8 |
          GET_FIELD_BT_TYPE(m_line, bcf_field_type, test_line_field_pos_idx);
      auto field_vector_diff = false;
      switch(field_bt_type_combo)
      {
        case BCF_BT_INT8 << 8 | BCF_BT_INT8:
          field_vector_diff = compare_unequal_vector<int8_t, int8_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT8 << 8 | BCF_BT_INT16:
          field_vector_diff = compare_unequal_vector<int8_t, int16_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT8 << 8 | BCF_BT_INT32:
          field_vector_diff = compare_unequal_vector<int8_t, int32_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT16 << 8 | BCF_BT_INT8:
          field_vector_diff = compare_unequal_vector<int16_t, int8_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT16 << 8 | BCF_BT_INT16:
          field_vector_diff = compare_unequal_vector<int16_t, int16_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT16 << 8 | BCF_BT_INT32:
          field_vector_diff = compare_unequal_vector<int16_t, int32_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT32 << 8 | BCF_BT_INT8:
          field_vector_diff = compare_unequal_vector<int32_t, int8_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT32 << 8 | BCF_BT_INT16:
          field_vector_diff = compare_unequal_vector<int32_t, int16_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_INT32 << 8 | BCF_BT_INT32:
          field_vector_diff = compare_unequal_vector<int32_t, int32_t>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_FLOAT << 8 | BCF_BT_FLOAT:
          field_vector_diff = compare_unequal_vector<float, float>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        case BCF_BT_CHAR << 8 | BCF_BT_CHAR:
          field_vector_diff = compare_unequal_vector<char, char>(gold_hdr, gold_line, bcf_field_type, i, test_line_field_pos_idx);
          break;
        default:
          throw VCFDiffException("Unknown type comparison "+std::to_string(field_bt_type_combo));
          break;
      }
      diff_fields = diff_fields || field_vector_diff;
      if(field_vector_diff)
        error_message += (std::string("ERROR: Gold and test differ in the ")+(bcf_field_type == BCF_HL_INFO ? "INFO" : "FORMAT")
            +" field "+bcf_hdr_int2id(gold_hdr, BCF_DT_ID, gold_field_idx)+"\n");
    }
  }
  for(auto i=0;i<GET_NUM_FIELDS(m_line, bcf_field_type);++i)
  {
    auto test_field_idx = GET_FIELD_IDX(m_line, bcf_field_type, i);
    //Skip comparing END values as the big outer loop handles that
    auto known_field_enum = m_field_idx_to_known_field_enum.get_known_field_enum_for_schema_idx(test_field_idx);
    if(m_field_idx_to_known_field_enum.is_defined_value(known_field_enum) && known_field_enum == GVCF_END_IDX)
      continue;
    auto lut_gold_field_idx = m_fields_lut.get_gold_idx_for_test(0, test_field_idx);
    assert(GoldLUT::is_missing_value(lut_gold_field_idx) || static_cast<size_t>(lut_gold_field_idx) < m_fields_in_gold_line.size());
    //Missing field in gold
    if(GoldLUT::is_missing_value(lut_gold_field_idx) || m_fields_in_gold_line[lut_gold_field_idx] < 0)
    {
      if(!bcf_is_empty_field(m_hdr, m_line, bcf_field_type, i))
      {
        diff_fields = true;
        error_message += (std::string("ERROR: ") + (bcf_field_type == BCF_HL_INFO ? "INFO" : "FORMAT")
          + " field " + bcf_hdr_int2id(m_hdr, BCF_DT_ID, test_field_idx) + " added in test line\n");
      }
    }
  }
  return diff_fields;
}

void VCFDiffFile::compare_line(const bcf_hdr_t* gold_hdr, bcf1_t* gold_line)
{
  //Ignore chr,pos as it's already handled previously
  bcf_unpack(m_line, BCF_UN_ALL);
  bcf_unpack(gold_line, BCF_UN_ALL);
  std::string error_message = "";
  auto diff_line_flag = false;
  //ID field
  auto diff_ID_flag = (strcmp(m_line->d.id, gold_line->d.id) != 0);
  error_message += (diff_ID_flag ? "ID field different\n" : "");
  diff_line_flag = diff_ID_flag || diff_line_flag ;
  //REF + ALT
  m_diff_alleles_flag = (m_line->n_allele != gold_line->n_allele);
  if(!m_diff_alleles_flag)
  {
    m_alleles_lut.resize_luts_if_needed(1, m_line->n_allele);
    std::unordered_map<std::string, int> allele_to_gold_idx;
    //Ignore REF
    for(auto i=1;i<gold_line->n_allele;++i)
      allele_to_gold_idx[gold_line->d.allele[i]] = i;
    m_alleles_lut.add_input_merged_idx_pair(0, 0, 0);   //REF mapping
    //Ignore REF
    for(auto i=1;i<m_line->n_allele;++i)
    {
      auto iter = allele_to_gold_idx.find(m_line->d.allele[i]);
      if(iter == allele_to_gold_idx.end())
      {
        m_diff_alleles_flag = true;
        break;
      }
      m_alleles_lut.add_input_merged_idx_pair(0, i, (*iter).second);
    }
    //Setup genotypes map
    if(!m_diff_alleles_flag)
    {
      auto num_alleles = m_line->n_allele;
      auto num_gts = (num_alleles*(num_alleles+1))/2;
      m_gold_genotype_idx_to_test_idx.resize(num_gts);
      for(auto i=0;i<gold_line->n_allele;++i)
      {
        auto lut_test_allele_idx_i = m_alleles_lut.get_input_idx_for_merged(0, i);
        for(auto j=i;j<gold_line->n_allele;++j)
        {
          auto lut_test_allele_idx_j = m_alleles_lut.get_input_idx_for_merged(0, j);
          m_gold_genotype_idx_to_test_idx[bcf_alleles2gt(i, j)] = bcf_alleles2gt(lut_test_allele_idx_i, lut_test_allele_idx_j);
        }
      }
    }
  }
  error_message += (m_diff_alleles_flag ? "Allele list different\n" : "");
  diff_line_flag = m_diff_alleles_flag || diff_line_flag;
  //QUAL
  auto diff_QUAL_flag = compare_unequal<float>(m_line->qual, gold_line->qual);
  error_message += (diff_QUAL_flag ? "QUAL field different\n" : "");
  diff_line_flag = diff_QUAL_flag || diff_line_flag;
  //FILTER
  auto diff_FILTER_flag = (m_line->d.n_flt != gold_line->d.n_flt);
  if(!diff_FILTER_flag)
  {
    reset_field_to_line_idx_mapping();
    for(auto i=0;i<gold_line->d.n_flt;++i)
      m_fields_in_gold_line[gold_line->d.flt[i]] = i;
    for(auto i=0;i<m_line->d.n_flt;++i)
    {
      auto lut_gold_flt_idx = m_fields_lut.get_gold_idx_for_test(0, m_line->d.flt[i]);
      assert(GoldLUT::is_missing_value(lut_gold_flt_idx) ||
          static_cast<size_t>(lut_gold_flt_idx) < m_fields_in_gold_line.size());
      if(GoldLUT::is_missing_value(lut_gold_flt_idx) || m_fields_in_gold_line[lut_gold_flt_idx] < 0)
      {
        diff_FILTER_flag = true;
        break;
      }
    }
  }
  error_message += (diff_FILTER_flag ? "FILTER list different\n" : "");
  diff_line_flag = diff_FILTER_flag || diff_line_flag;
  //INFO fields
  auto  diff_INFO_fields = compare_unequal_fields(gold_hdr, gold_line, BCF_HL_INFO, error_message);
  diff_line_flag = diff_INFO_fields || diff_line_flag;
  //FORMAT fields
  auto diff_FORMAT_fields = compare_unequal_fields(gold_hdr, gold_line, BCF_HL_FMT, error_message);
  diff_line_flag = diff_FORMAT_fields || diff_line_flag;
  if(diff_line_flag)
  {
    std::cerr << "=====================================================================\n";
    m_tmp_hts_string.l = 0;
    vcf_format(gold_hdr, gold_line, &m_tmp_hts_string);
    std::cerr << m_tmp_hts_string.s << "\n";
    print_line();
    std::cerr << error_message;
    std::cerr << "=====================================================================\n";
  }
}

std::string VCFDiffFile::create_region(const std::string& regions,
    const std::unordered_map<std::string, std::pair<int64_t, int64_t>>& regions_contig_to_interval, const std::string& contig)
{
  if(regions.length())
  {
    auto iter = regions_contig_to_interval.find(contig);
    if(iter == regions_contig_to_interval.end())
      return "";
    auto pair = (*iter).second;
    if(pair.first == 1ll && pair.second == INT64_MAX)
      return (contig + ",");    //full contig
    if(pair.second == INT64_MAX)
    {
      auto contig_idx = bcf_hdr_name2id(m_hdr, contig.c_str());
      VERIFY_OR_THROW(contig_idx >= 0);
      auto contig_length = bcf_hdr_id2contig_length(m_hdr, contig_idx);
      return (contig + ":" + std::to_string(pair.first) + "-" + std::to_string(contig_length) + ",");
    }
    else
      return (contig + ":" + std::to_string(pair.first) + "-" + std::to_string(pair.second) + ",");
  }
  else
    return (contig + ",");
}

/*
   Set regions to traverse based on contigs in the header
*/
void set_regions(VCFDiffFile& gold, VCFDiffFile& test, const std::string& regions="")
{
  const auto& gold_contigs = gold.m_contigs;
  const auto& test_contigs = test.m_contigs;
  std::vector<std::string> common_contigs(std::min(gold_contigs.size(), test_contigs.size()));
  {
    //Find common contigs
    auto iter = std::set_intersection(gold_contigs.begin(), gold_contigs.end(), test_contigs.begin(), test_contigs.end(),
        common_contigs.begin());
    common_contigs.resize(iter-common_contigs.begin());
  }
  //Sort as strings
  std::sort(common_contigs.begin(), common_contigs.end());
  //Regions set - parse regions string
  std::unordered_map<std::string, std::pair<int64_t, int64_t>> regions_contig_to_interval;
  if(regions.length())
  {
    std::set<std::string> regions_contig_set;
    //regions of the form 'chr'|'chr:pos'|'chr:from-to'|'chr:from-[,...]
    //Same format as bcftools
    std::istringstream f1(regions);
    std::string contig_region;
    while(std::getline(f1, contig_region, ','))
    {
      std::string contig;
      int64_t begin = 1;
      int64_t end = INT64_MAX;
      std::istringstream f2(contig_region);
      std::string tmp_s;
      auto idx = 0u;
      while(std::getline(f2, tmp_s, ':'))
      {
        if(idx == 0u)
          contig = tmp_s;
        else
        {
          std::istringstream f3(tmp_s);
          std::string boundaries;
          auto idx2 = 0u;
          while(std::getline(f3, boundaries, '-'))
          {
            if(idx2 == 0u)
            {
              begin = strtoll(boundaries.c_str(), 0, 10);
              if(tmp_s.find("-") == std::string::npos)      //"-" not found
                end = begin;
            }
            else
              end = strtoll(boundaries.c_str(), 0, 10);
            ++idx2;
          }
        }
        ++idx;
      }
      if(regions_contig_set.find(contig) == regions_contig_set.end())
      {
        regions_contig_set.insert(contig);
        regions_contig_to_interval[contig] = std::pair<int64_t, int64_t>(begin, end);
      }
    }
    //Find common contigs among gold, test, regions
    auto copy_common = common_contigs;
    {
      auto iter = std::set_intersection(copy_common.begin(), copy_common.end(), regions_contig_set.begin(), regions_contig_set.end(),
          common_contigs.begin());
      common_contigs.resize(iter-common_contigs.begin());
    }
    std::sort(common_contigs.begin(), common_contigs.end());
  }
  //Contigs only in gold
  std::vector<std::string> only_gold(gold_contigs.size());
  {
    auto iter = std::set_difference(gold_contigs.begin(), gold_contigs.end(), common_contigs.begin(), common_contigs.end(),
        only_gold.begin());
    only_gold.resize(iter - only_gold.begin());
  }
  gold.m_regions="";
  for(auto& x : common_contigs)
    gold.m_regions += (gold.create_region(regions, regions_contig_to_interval, x));
  for(auto& x : only_gold)
    gold.m_regions += (gold.create_region(regions, regions_contig_to_interval, x));
  if(gold.m_regions.length())
    gold.m_regions.pop_back();  //remove trailing ,
  //Contigs only in test
  std::vector<std::string> only_test(test_contigs.size());
  {
    auto iter = std::set_difference(test_contigs.begin(), test_contigs.end(), common_contigs.begin(), common_contigs.end(),
        only_test.begin());
    only_test.resize(iter - only_test.begin());
  }
  test.m_regions="";
  for(auto& x : common_contigs)
    test.m_regions += (test.create_region(regions, regions_contig_to_interval, x));
  for(auto& x : only_test)
    test.m_regions += (test.create_region(regions, regions_contig_to_interval, x));
  if(test.m_regions.length())
    test.m_regions.pop_back();  //remove trailing ,
  gold.set_regions_and_open_file();
  test.set_regions_and_open_file();
}

void construct_regions_for_partitions(const std::string& loader_json_filename, VidMapper*& vid_mapper,
    JSONConfigBase& json_config_base, const int rank,
    VCFDiffFile& gold, VCFDiffFile& test, std::string& regions)
{
  regions = "";
  //Parse JSON
  std::ifstream ifs(loader_json_filename.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  rapidjson::Document json;
  json.Parse(str.c_str());
  VERIFY_OR_THROW(json.HasMember("vid_mapping_file") && json["vid_mapping_file"].IsString());
  VERIFY_OR_THROW(json.HasMember("column_partitions"));
  //callsets file
  std::string callset_mapping_file = (json.HasMember("callset_mapping_file") && json["callset_mapping_file"].IsString())
    ? json["callset_mapping_file"].GetString() : "";
  vid_mapper = static_cast<VidMapper*>(new FileBasedVidMapper(json["vid_mapping_file"].GetString(), callset_mapping_file));
  //Limit #callsets
  if(json.HasMember("limit_callset_row_idx") && json["limit_callset_row_idx"].IsInt64())
    g_num_callsets = json["limit_callset_row_idx"].GetInt64() + 1;
  else
    g_num_callsets = vid_mapper->get_num_callsets();
  //Get partitions
  json_config_base.read_from_file(loader_json_filename);
  auto column_interval = json_config_base.get_column_partition(rank);
  //Find all contig regions in the current partition
  auto current_position = column_interval.first;
  std::string contig_name;
  while(current_position <= column_interval.second)
  {
    int64_t contig_position = -1ll;
    auto status = vid_mapper->get_contig_location(current_position, contig_name, contig_position);
    if(status)
    {
      regions += contig_name;
      regions += (":" + std::to_string(contig_position+1)+"-");    //regions string is 1-based
      ContigInfo info;
      VERIFY_OR_THROW(vid_mapper->get_contig_info(contig_name, info));
      //partition ends here
      if(info.m_tiledb_column_offset + info.m_length > column_interval.second)
        regions += std::to_string(column_interval.second - info.m_tiledb_column_offset + 1); //regions string is 1-based
      regions += ",";
    }
    status = vid_mapper->get_next_contig_location(current_position, contig_name, current_position);
    if(!status) //no more contigs found
      break;
  }
  if(regions.length())   //delete last comma
    regions.pop_back();
}

enum ArgsIdxEnum
{
  ARGS_USE_CALLSETS_FILE_FOR_SAMPLE_IDX=1000
};

void setup_samples_lut(const std::string& test_to_gold_callset_map_file, VCFDiffFile& gold, VCFDiffFile& test, const VidMapper* vid_mapper)
{
  //Parse JSON
  std::ifstream ifs(test_to_gold_callset_map_file.c_str());
  VERIFY_OR_THROW(ifs.is_open());
  std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  rapidjson::Document json;
  json.Parse(str.c_str());
  VERIFY_OR_THROW(json.HasMember("test_to_gold_callset_map") && json["test_to_gold_callset_map"].IsObject());
  const auto& callset_map = json["test_to_gold_callset_map"];
  //if(vid_mapper->get_num_callsets() != callset_map.MemberCount())
    //throw VCFDiffException("ERROR: Mismatch in the #samples in the callsets JSON file and test-to-gold callset mapping file. Callsets JSON has "
        //+std::to_string(vid_mapper->get_num_callsets())+" samples while mapping file has "+std::to_string(callset_map.MemberCount())+
        //" - cannot perform any further comparisons");
  //Everyone should have the same number of samples
  if(g_num_callsets != bcf_hdr_nsamples(gold.m_hdr))
    throw VCFDiffException("ERROR: Mismatch in the #samples in the gold VCF and the callsets JSON file. Gold has "+
        std::to_string(bcf_hdr_nsamples(gold.m_hdr))+" samples while callsets JSON has "+std::to_string(g_num_callsets)+
        " - cannot perform any further comparisons");
  if(g_num_callsets != bcf_hdr_nsamples(test.m_hdr))
    throw VCFDiffException("ERROR: Mismatch in the #samples in the test VCF and the callsets JSON file. Test has "+
        std::to_string(bcf_hdr_nsamples(test.m_hdr))+" samples while callsets JSON has "+std::to_string(g_num_callsets)+
        " - cannot perform any further comparisons");
  auto num_samples_found = 0u;
  for(auto b=callset_map.MemberBegin(), e=callset_map.MemberEnd();b!=e;++b)
  {
    VERIFY_OR_THROW((*b).value.IsString());
    auto test_sample_name = (*b).name.GetString();
    auto gold_sample_name = (*b).value.GetString();
    auto gold_sample_idx = bcf_hdr_id2int(gold.m_hdr, BCF_DT_SAMPLE, gold_sample_name);
    auto test_sample_idx = bcf_hdr_id2int(test.m_hdr, BCF_DT_SAMPLE, test_sample_name);
    if(gold_sample_idx >= 0 && test_sample_idx >= 0)
    {
      test.m_samples_lut.add_input_merged_idx_pair(0, test_sample_idx, gold_sample_idx);
      ++num_samples_found;
    }
  }
  VERIFY_OR_THROW(num_samples_found == g_num_callsets && "Test-to-gold callset mapping file does not have mapping for all samples");
}

int main(int argc, char** argv)
{
#ifdef HTSDIR
  //Initialize MPI environment
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get my world rank
  int my_world_mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
  static struct option long_options[] =
  {
    {"threshold",1,0,'t'},
    {"regions",1,0,'r'},
    {"loader-config",1,0,'l'},
    {"process-rank",1,0,'p'},
    {"test_to_gold_callset_map_file",1,0, ARGS_USE_CALLSETS_FILE_FOR_SAMPLE_IDX}
  };
  int c = -1;
  std::string regions = "";
  std::string loader_json_filename = "";
  std::string test_to_gold_callset_map_file = "";
  while((c=getopt_long(argc, argv, "t:r:l:p:", long_options, NULL)) >= 0)
  {
    switch(c)
    {
      case 't':
        g_threshold = strtod(optarg, 0);
        break;
      case 'r':
        regions = std::move(std::string(optarg));
        break;
      case 'l':
        loader_json_filename = std::move(std::string(optarg));
        break;
      case 'p':
        my_world_mpi_rank = strtoll(optarg, 0, 10);
        break;
      case ARGS_USE_CALLSETS_FILE_FOR_SAMPLE_IDX:
        test_to_gold_callset_map_file = std::move(std::string(optarg));
        break;
      default:
        throw VCFDiffException(std::string("Unknown argument: ")+argv[optind-1]);
        break;
    }
  }
  if(optind+2 > argc)
  {
    std::cerr << "Needs 2 VCF files as input <gold> <test>\n";
    exit(-1);
  }
  VCFDiffFile gold(argv[optind]);
  VCFDiffFile test(argv[optind+1]);
  auto use_loader_json_file = (loader_json_filename.length() && regions.length() == 0u);
  //Setup luts
  test.setup_luts(gold, use_loader_json_file && (test_to_gold_callset_map_file.length() > 0u));
  //Loader input json - compare partitions created by vcf2tiledb
  VidMapper* vid_mapper = 0;
  JSONConfigBase json_config_base;
  if(use_loader_json_file)
  {
    construct_regions_for_partitions(loader_json_filename, vid_mapper, json_config_base, my_world_mpi_rank,
        gold, test, regions);
    if(test_to_gold_callset_map_file.length())
      setup_samples_lut(test_to_gold_callset_map_file, gold, test, vid_mapper);
  }
  //Regions
  set_regions(gold, test, regions);
  bool have_data = test.m_line && gold.m_line;
  //Control variables
  int curr_gold_contig_idx_in_vid = -1;
  ContigInfo info;
  auto both_must_have_valid_lines = false;
  while(have_data)
  {
    auto lut_gold_contig_idx = test.m_contigs_lut.get_gold_idx_for_test(0, test.m_line->rid);
    auto lut_test_contig_idx = test.m_contigs_lut.get_test_idx_for_gold(0, gold.m_line->rid);
    both_must_have_valid_lines = false;
    //Different contig
    while(have_data && gold.m_line->rid != lut_gold_contig_idx)
    {
      //Gold has reached a contig not in test
      if(GoldLUT::is_missing_value(lut_test_contig_idx))
      {
        std::cerr << "ERROR: Gold vcf has extra line(s) - printing the first extra line\n";
        gold.print_line();
        break;    //because common contigs are handled first
      }
      //test has reached a contig not in gold
      if(GoldLUT::is_missing_value(lut_gold_contig_idx))
      {
        std::cerr << "ERROR: Test vcf has extra line(s) - printing the first extra line\n";
        test.print_line();
        break;    //because common contigs are handled first
      }
      //Both have entries in common_contigs, gold is at "lower" contig
      if(std::string(bcf_hdr_id2name(gold.m_hdr, gold.m_line->rid)) < std::string(bcf_hdr_id2name(gold.m_hdr, lut_gold_contig_idx)))
      {
        std::cerr << "ERROR: Gold vcf has extra line(s) - printing the first extra line\n";
        gold.print_line();
        gold.seek_and_read(lut_gold_contig_idx, test.m_line->pos);  //gold seeks to test's position
        lut_test_contig_idx = gold.m_line ? test.m_contigs_lut.get_test_idx_for_gold(0, gold.m_line->rid) : lut_missing_value;
      }
      else      //test seeks to gold's position
      {
        std::cerr << "ERROR: Test vcf has extra line(s) - printing the first extra line\n";
        test.print_line();
        test.seek_and_read(lut_test_contig_idx, gold.m_line->pos);
        lut_gold_contig_idx = test.m_line ? test.m_contigs_lut.get_gold_idx_for_test(0, test.m_line->rid) : lut_missing_value; 
      }
      have_data = test.m_line && gold.m_line;
    }
    //Same contig now
    if(have_data && gold.m_line->rid == lut_gold_contig_idx)
    {
      auto do_compare = true;
      //Set END points
      gold.m_line->m_end_point = -1;
      test.m_line->m_end_point = -1; 
      //Unpacks INFO fields
      bcf_unpack(gold.m_line, BCF_UN_INFO);
      bcf_unpack(test.m_line, BCF_UN_INFO);
      bcf_set_end_point_from_info(gold.m_hdr, gold.m_line);
      bcf_set_end_point_from_info(test.m_hdr, test.m_line);
      //Variant types
      bcf_get_variant_types(gold.m_line);
      bcf_get_variant_types(test.m_line);
      //NON_REF blocks
      auto gold_is_NON_REF_block = (gold.m_line->d.var_type == VCF_NON_REF);
      auto test_is_NON_REF_block = (test.m_line->d.var_type == VCF_NON_REF);
      auto partial_overlap = ( (gold.m_line->pos <= test.m_line->pos && gold.m_line->m_end_point >= test.m_line->pos)
              || (test.m_line->pos <= gold.m_line->pos && test.m_line->m_end_point >= gold.m_line->pos) );
      auto full_overlap = (gold.m_line->pos == test.m_line->pos && gold.m_line->m_end_point == test.m_line->m_end_point);
      //Different positions
      if(!full_overlap)
      {
        do_compare = false;
        auto before_begin = false;
        auto after_end = false;
        //Check if variant is before or after the specified interval
        if(vid_mapper)
        {
          if(curr_gold_contig_idx_in_vid != gold.m_line->rid)
          {
            auto status = vid_mapper->get_contig_info(bcf_hdr_id2name(gold.m_hdr, gold.m_line->rid), info);
            curr_gold_contig_idx_in_vid = status ? gold.m_line->rid : -1;
          }
          if(curr_gold_contig_idx_in_vid == gold.m_line->rid)
          {
            //Variant starts before column partition
            if((info.m_tiledb_column_offset + static_cast<int64_t>(gold.m_line->pos)) <
                json_config_base.get_column_partition(my_world_mpi_rank).first)
              before_begin = true;
            //Variant past the column interval
            if((info.m_tiledb_column_offset + static_cast<int64_t>(gold.m_line->m_end_point)) >
                json_config_base.get_column_partition(my_world_mpi_rank).second)
              after_end = true;
          }
        }
        //Overlap and if this is a NON_REF interval, print warning and do comparison
        if(gold_is_NON_REF_block && test_is_NON_REF_block && partial_overlap)
        {
          if(!before_begin && !after_end)
          {
            std::cout << "WARNING: Gold and test REF blocks overlap, but do not match exactly at position "<<
              bcf_hdr_id2name(test.m_hdr, test.m_line->rid) << ","<< (std::min(test.m_line->pos, gold.m_line->pos)+1)<<"\n";
            gold.print_line(std::cout);
            test.print_line(std::cout);
          }
          do_compare = true;
        }
        else
        {
          if(!before_begin)
          {
            std::cerr << "ERROR: Lines with different positions found - resetting file ptr to next match position\n";
            gold.print_line();
            test.print_line();
          }
          //No overlap
          if(!partial_overlap)
          {
            if(test.m_line->pos > gold.m_line->pos)
              gold.seek_and_read(gold.m_line->rid, test.m_line->pos);
            else
              test.seek_and_read(test.m_line->rid, gold.m_line->pos);
          }
          else          //partial overlap
          {
            if(test.m_line->m_end_point > gold.m_line->m_end_point)
            {
              //test contains deletion and gold contains spanning deletion, advance test
              if(test.m_line->rlen > 1 && (test.m_line->d.var_type|VCF_INDEL)
                  && (gold.m_line->d.var_type|VCF_SPANNING_DELETION) && test.m_line->pos <= gold.m_line->pos)
              {
                test.read_and_advance();
                both_must_have_valid_lines = before_begin;      //if !before_begin, then error message is already printed
              }
              else
                gold.read_and_advance();
            }
            else
              if(gold.m_line->m_end_point > test.m_line->m_end_point)
              {
                //gold contains deletion and test contains spanning deletion, advance gold
                if(gold.m_line->rlen > 1 && (gold.m_line->d.var_type|VCF_INDEL)
                    && (test.m_line->d.var_type|VCF_SPANNING_DELETION) && gold.m_line->pos <= test.m_line->pos)
                {
                  gold.read_and_advance();
                  both_must_have_valid_lines = before_begin;      //if !before_begin, then error message is already printed
                }
                else
                  test.read_and_advance();
              }
              else      //equal end points - advance both
              {
                gold.read_and_advance();
                test.read_and_advance();
                both_must_have_valid_lines = true;
              }
          }
        }
      }
      if(do_compare)
      {
        test.compare_line(gold.m_hdr, gold.m_line);
        if(full_overlap)
        {
          gold.read_and_advance();
          test.read_and_advance();
          both_must_have_valid_lines = true;
        }
        else    //must be partial_overlap
        {
          assert(partial_overlap);
          if(test.m_line->m_end_point > gold.m_line->m_end_point)
            gold.read_and_advance();
          else
            if(gold.m_line->m_end_point > test.m_line->m_end_point)
              test.read_and_advance();
            else      //equal end points - advance both
            {
              gold.read_and_advance();
              test.read_and_advance();
              both_must_have_valid_lines = true;
            }
        }
      }
      have_data = test.m_line && gold.m_line;
    }
    else
      break;    //either no more data or moved to contigs not present in the other file
  }
  have_data = test.m_line && gold.m_line;
  //Print error for the case where one of them has run out of lines, but the other has lines
  if(!have_data)
  {
    for(auto i=0u;i<2u;++i)
    {
      auto& diff_ref = (i == 0u) ? gold : test;
      auto name = (i == 0u) ? "Gold" : "Test";
      if(diff_ref.m_line)
      {
        if(both_must_have_valid_lines)
        {
          std::cerr << "ERROR: "<< name << " vcf has extra line(s) - printing the first extra line\n";
          diff_ref.print_line();
        }
        else
        {
          diff_ref.read_and_advance();
          if(diff_ref.m_line)
          {
            std::cerr << "ERROR: "<< name << " vcf has extra line(s) - printing the first extra line\n";
            diff_ref.print_line();
          }
        }
      }
    }
  }
  if(vid_mapper)
    delete vid_mapper;
  MPI_Finalize();
#else //ifdef HTSDIR
  std::cerr << "vcf diff needs htslib - recompile with HTSDIR set\n";
#endif  //ifdef HTSDIR
  return 0;
}
