/**
 * @file   query_processor.cc
 * @author Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * @section LICENSE
 *
 * The MIT License
 * 
 * @copyright Copyright (c) 2014 Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * This file implements the QueryProcessor class.
 */
  
#include "query_processor.h"
#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>

#include<algorithm>
#include<unordered_map>
#include "lut.h"

#ifndef HTSDIR
#include "vcf.h"
#else
#include "htslib/vcf.h"
#endif

/* ------------------- Genotyping functions ------------------ */

/******************************************************
********************* CONSTRUCTORS ********************
******************************************************/

QueryProcessor::GTColumn::GTColumn(int64_t col, uint64_t row_num) {
  ALT_.resize(row_num);
  col_ = col;
  REF_.resize(row_num);
  PL_.resize(row_num);
}

void QueryProcessor::GTColumn::reset()
{
  for(auto i=0ull;i<REF_.size();++i)
  {
    REF_[i].clear();
    ALT_[i].clear();
    PL_[i].clear();
  }
}

/******************************************************
*********************** OPERATIONS ********************
******************************************************/

std::string g_non_reference_allele = "<NON_REF>";
#define CHECK_MISSING_SAMPLE_GIVEN_REF(REF) (((REF).size() == 0) || ((REF)[0] == '$'))
#define CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(REF) ((REF).size() == 0 || ((REF) == ""))
#define CHECK_IN_THE_MIDDLE_REF(REF) ((REF)[0] == 'N')
#define IS_NON_REF_ALLELE(allele) ((allele)[0] == '&')
// TODO: Implement function for deriving the final ALT and PL values
/*
 * @brief - get the longest reference allele among all variants at this position and store its value in merged_reference_allele
 * For example, if we have the reference alleles T (SNP) and TG (deletion) in two GVCFs at the same location, the reference allele
 * in the merged GVCF should be TG. This modifies the alt alleles (as explained below)
 * Does a sanity check - ref alleles should be a prefix of the merged
 * @param reference_vector - vector of REF strings
 * @param merged_reference_allele - string to store the merged reference
 */
void merge_reference_allele(const std::vector<std::string>& reference_vector, std::string& merged_reference_allele)
{
  auto* longer_ref = &merged_reference_allele;
  auto merged_ref_length = merged_reference_allele.length();
  for(const auto& curr_ref : reference_vector)
  {
    if(CHECK_MISSING_SAMPLE_GIVEN_REF(curr_ref))
      continue;
    auto curr_ref_length = curr_ref.length();
    auto* shorter_ref = &curr_ref;
    const auto is_curr_ref_longer = (curr_ref_length > merged_ref_length);
    if(is_curr_ref_longer)
    {
      longer_ref = const_cast<std::string*>(&curr_ref);
      shorter_ref = &merged_reference_allele;
    }
#ifdef DEBUG
    //sanity check only - the shorter ref must be a prefix of the longer ref (since they begin at the same location)
    if(longer_ref->find(*shorter_ref) != 0)
    {
      throw std::invalid_argument(std::string{"When combining variants at a given position, the shorter reference allele should be a prefix of the longer reference allele: \'"} + *shorter_ref + " , " + *longer_ref);
    }
#endif
    if(is_curr_ref_longer)
    {
      if(curr_ref_length >= merged_reference_allele.capacity())
	merged_reference_allele.reserve(2*curr_ref_length+1);	//why 2, why not?
      if(merged_ref_length > 0 && CHECK_IN_THE_MIDDLE_REF(merged_reference_allele))
        merged_reference_allele = curr_ref;
      else      //append remaining chars to merged reference
        merged_reference_allele.append(curr_ref, merged_ref_length, curr_ref_length - merged_ref_length);
      merged_ref_length = curr_ref_length;
      longer_ref = &merged_reference_allele;
    }
    else
      if(CHECK_IN_THE_MIDDLE_REF(merged_reference_allele) && !CHECK_IN_THE_MIDDLE_REF(curr_ref))
        merged_reference_allele = curr_ref;
  }
}

/*
 * @brief - remap alt alleles for each input to the longer merged reference
 * For example, if we have the reference alleles T (SNP) and TG (deletion) in two GVCFs at the same location, the reference allele
 * in the merged GVCF should be TG. If the alt alleles were G and T respectively, then the alt alleles in the merged variant
 * become GG,T
 * @param reference_vector - vector of REF strings 
 * @param alt_alleles - vector of vector of alt allele strings
 * @param merged_reference_allele - the merged reference produced by merge_reference_allele
 * @param alleles_LUT -  LUT containing mapping of allele idxs between merged variant and input variants
 * @return vector of merged alt allele strings
 */
const std::vector<std::string> merge_alt_alleles(const std::vector<std::string>& reference_vector, 
    const std::vector<std::vector<std::string>> alt_alleles,
    const std::string& merged_reference_allele,
    CombineAllelesLUT& alleles_LUT) {
  // marking non_reference_allele as already seen will ensure it's not included in the middle
  auto seen_alleles = std::unordered_map<std::string, int>{{g_non_reference_allele,-1}};
  auto merged_alt_alleles = std::vector<std::string>{};
  auto merged_reference_length = merged_reference_allele.length();
  //invalidate all existing mappings in the LUT
  alleles_LUT.reset_luts();
  //vector to store idx mappings for NON_REF allele, update LUT at end as #ALT alleles are not known till end
  //Set everything to -1 (invalid mapping)
  auto input_non_reference_allele_idx = std::vector<int>(reference_vector.size(), -1);
  auto merged_allele_idx = 1u;	//why 1, ref is index 0, alt begins at 1
  auto input_sample_idx = 0u;
  for (const auto& curr_reference : reference_vector)
  {
    if(!CHECK_MISSING_SAMPLE_GIVEN_REF(curr_reference))
    {
      const auto& curr_reference_length = curr_reference.length();
      const auto& curr_allele_vector = alt_alleles[input_sample_idx];
      auto is_suffix_needed = false;
      auto suffix_length = 0u;
      if(curr_reference_length < merged_reference_length)
      {
        is_suffix_needed = true;
        suffix_length = merged_reference_length - curr_reference_length;
      }
      //mapping for reference allele 0 -> 0
      alleles_LUT.add_input_merged_idx_pair(input_sample_idx, 0, 0);
      auto input_allele_idx = 1u;	//why 1, ref is index 0, alt begins at 1
      //copy of allele if needed
      std::string copy_allele;
      for (const auto& allele : curr_allele_vector)
      {
        if(IS_NON_REF_ALLELE(allele))
          input_non_reference_allele_idx[input_sample_idx] = input_allele_idx;
        else
        {
          auto* allele_ptr = &allele;
          if(is_suffix_needed)
          {
            copy_allele = allele;
            copy_allele.append(merged_reference_allele, curr_reference_length, suffix_length);
            allele_ptr = &copy_allele; //allele_ptr now points to a copy of var.alt()[l] (+suffix), hence it's safe to move later
          }
          const auto& iter_pos = seen_alleles.find(*allele_ptr);
          if (iter_pos == seen_alleles.end()) { //allele seen for the first time
            seen_alleles[*allele_ptr] = merged_allele_idx;
            //always check whether LUT is big enough for alleles_LUT (since the #alleles in the merged variant is unknown)
            //Most of the time this function will return quickly (just a if condition check)
            alleles_LUT.resize_luts_if_needed(merged_allele_idx + 1); 
            alleles_LUT.add_input_merged_idx_pair(input_sample_idx, input_allele_idx, merged_allele_idx);
            if(is_suffix_needed)
              merged_alt_alleles.push_back(std::move(*allele_ptr)); //allele_ptr points to a copy - use move
            else
              merged_alt_alleles.push_back(*allele_ptr);	//allele_ptr points to curr_allele_vector[input_allele_idx] - copy to vector
            ++merged_allele_idx;
          }
          else
            alleles_LUT.add_input_merged_idx_pair(input_sample_idx, input_allele_idx, (*iter_pos).second);
        }
        ++input_allele_idx;
      }
    }
    ++input_sample_idx;
  }
  // always want non_reference_allele to be last
  merged_alt_alleles.push_back(g_non_reference_allele);
  auto non_reference_allele_idx = merged_alt_alleles.size(); //why not -1, include reference allele also
  //always check whether LUT is big enough for alleles_LUT (since the #alleles in the merged variant is unknown)
  alleles_LUT.resize_luts_if_needed(non_reference_allele_idx + 1); 
  //Add mappings for non_ref allele
  input_sample_idx = 0u;
  for (const auto& curr_reference : reference_vector)
  {
    if(!CHECK_MISSING_SAMPLE_GIVEN_REF(curr_reference))
    {
      if(input_non_reference_allele_idx[input_sample_idx] >= 0)
        alleles_LUT.add_input_merged_idx_pair(input_sample_idx, input_non_reference_allele_idx[input_sample_idx], non_reference_allele_idx);
    }
    ++input_sample_idx;
  }
  return merged_alt_alleles;
}

/*
  Copied from defunct gamgee library
  @param pls - vector of PL values for a given sample as stored in TileDB
  @param input_sample_idx 
  @param alleles_LUT LUT mapping alleles from input to merged alleles list
  @param num_merged_alleles  
  @param remapped_pls - matrix of PLs, each row corresponds to a single genotype, each column corresponds to a sample
  @num_missing_samples - keeps track of how many samples had unknown values for a given genotype
 */
void remap_pl(const std::vector<int32_t>& pls, const uint32_t input_sample_idx,
    const CombineAllelesLUT& alleles_LUT, const unsigned num_merged_alleles,
    std::vector<std::vector<int>>& remapped_pls,
    std::vector<uint64_t>& num_missing_samples) {
  //index of NON_REF in merged variant
  const auto merged_non_reference_allele_idx = static_cast<int>(num_merged_alleles-1);
  //index of NON_REF in input sample
  const auto input_non_reference_allele_idx = alleles_LUT.get_input_idx_for_merged(input_sample_idx, merged_non_reference_allele_idx);

  for (auto allele_j = 0u; allele_j < num_merged_alleles; ++allele_j) {
    auto input_j_allele = alleles_LUT.get_input_idx_for_merged(input_sample_idx, allele_j);
    if (CombineAllelesLUT::is_missing_value(input_j_allele))	//no mapping found for current allele in input gvcf
    {
      if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
      {
	//fill in missing values for all genotypes with allele_j as one component
	for(auto allele_k = allele_j; allele_k < num_merged_alleles;++allele_k)
        {
          auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
	  remapped_pls[gt_idx][input_sample_idx] = bcf_int32_missing;
          ++(num_missing_samples[gt_idx]);
        }
	continue;	//skip to next value of allele_j
      }
      else //input contains NON_REF allele, use its idx
	input_j_allele = input_non_reference_allele_idx;
    }
    for (auto allele_k = allele_j; allele_k < num_merged_alleles; ++allele_k) {
      auto gt_idx = bcf_alleles2gt(allele_j, allele_k);
      auto input_k_allele = alleles_LUT.get_input_idx_for_merged(input_sample_idx, allele_k);
      if (CombineAllelesLUT::is_missing_value(input_k_allele))	//no mapping found for current allele in input gvcf
      {
	if(CombineAllelesLUT::is_missing_value(input_non_reference_allele_idx))	//input did not have NON_REF allele
	{
	  remapped_pls[gt_idx][input_sample_idx] = bcf_int32_missing; //put missing value
          ++(num_missing_samples[gt_idx]);
	  continue;	//skip to next value of allele_k
	}
	else //input has NON_REF, use its idx
	  input_k_allele = input_non_reference_allele_idx;
      }
      remapped_pls[gt_idx][input_sample_idx] = pls[bcf_alleles2gt(input_j_allele, input_k_allele)];
    }
  }
}

// TODO: Implement the genotyping function
void do_dummy_genotyping(const QueryProcessor::GTColumn* gt_column, std::ostream& output)
{
  std::string merged_reference_allele;
  merged_reference_allele.reserve(10);
  merge_reference_allele(gt_column->REF_, merged_reference_allele);

  //initialize to number of samples
  CombineAllelesLUT alleles_LUT { static_cast<unsigned>(gt_column->REF_.size()) };
  auto& merged_alt_alleles = merge_alt_alleles(gt_column->REF_, gt_column->ALT_, merged_reference_allele, alleles_LUT);

  //Allocate space for remapped PL
  auto num_samples = gt_column->REF_.size();
  auto num_merged_alleles = merged_alt_alleles.size() + 1u;
  auto num_gts = (num_merged_alleles*(num_merged_alleles+1))/2;
  std::vector<std::vector<int>> remapped_PLs = std::vector<std::vector<int>>(num_gts, std::vector<int>(num_samples, bcf_int32_missing));
  std::vector<uint64_t> num_missing_samples = std::vector<uint64_t>(num_gts, 0ull);

  //Remap PL
  auto input_sample_idx = 0u;
  for(const auto& input_pl_vector : gt_column->PL_)
  {
    if(!CHECK_MISSING_SAMPLE_GIVEN_REF(gt_column->REF_[input_sample_idx]))
      remap_pl(input_pl_vector, input_sample_idx,
          alleles_LUT, num_merged_alleles,
          remapped_PLs,  num_missing_samples);
    else
      for(auto i=0u;i<num_gts;++i)
        ++(num_missing_samples[i]);
    ++input_sample_idx;
  }
  //Compute medians
  std::vector<int> median_vector;
  median_vector.resize(num_gts);
  for(auto i=0u;i<num_gts;++i)
    if(num_missing_samples[i] == num_samples)
      median_vector[i] = bcf_int32_missing;
    else
    {
      auto& curr_PL_vector = remapped_PLs[i];
      auto dec_order_median_idx = (num_samples - num_missing_samples[i])/2;
      //auto inc_order_median_idx = num_missing_samples[i] + (num_samples - num_missing_samples[i])/2;
      std::nth_element(curr_PL_vector.begin(), curr_PL_vector.begin() + dec_order_median_idx, curr_PL_vector.end(), std::greater<int>());
      //std::nth_element(curr_PL_vector.begin(), curr_PL_vector.begin() + inc_order_median_idx, curr_PL_vector.end());
      median_vector[i] = curr_PL_vector[dec_order_median_idx];
      //median_vector[i] = curr_PL_vector[inc_order_median_idx];
      assert(median_vector[i] != bcf_int32_missing);
    }
  output << gt_column->col_ << ",";
  output << merged_reference_allele;
  for(const auto& curr_alt_allele : merged_alt_alleles)
    output << "," << curr_alt_allele;
  for(auto value : median_vector)
    output << "," << value;
  output << "\n";
  return;
}

/* ----------------- QueryProcessor functions ---------------- */

/******************************************************
********************* CONSTRUCTORS ********************
******************************************************/

QueryProcessor::QueryProcessor(const std::string& workspace, 
                               StorageManager& storage_manager) 
    : storage_manager_(storage_manager) {
  set_workspace(workspace);
  create_workspace(); 
}

void QueryProcessor::export_to_CSV(const StorageManager::ArrayDescriptor* ad,
                                   const std::string& filename) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::string& array_name = array_schema.array_name();
  const unsigned int attribute_num = array_schema.attribute_num();
  const unsigned int dim_num = array_schema.dim_num();
  
  // Prepare CSV file
  CSVFile csv_file(filename, CSVFile::WRITE);

  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its = 
      new StorageManager::const_iterator[attribute_num+1];
  StorageManager::const_iterator tile_it_end;
  initialize_tile_its(ad, tile_its, tile_it_end);

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Iterate over all tiles
  while(tile_its[attribute_num] != tile_it_end) {
    // Iterate over all cells of each tile
    initialize_cell_its(tile_its, attribute_num, cell_its, cell_it_end);

    while(cell_its[attribute_num] != cell_it_end) { 
      csv_file << cell_to_csv_line(cell_its, attribute_num);
      advance_cell_its(attribute_num, cell_its);
    }
 
    advance_tile_its(attribute_num, tile_its);
  }

  // Clean up 
  delete [] tile_its;
  delete [] cell_its;
}

void QueryProcessor::handle_gvcf_ranges(gVCF_PQ& end_pq, std::vector<PQStruct>& PQ_end_vec, QueryProcessor::GTColumn* gt_column,
      std::unordered_map<uint64_t, GTTileIteratorsTracker>& tile_idx_2_iters, std::ostream& output_stream,
    int64_t current_start_position, int64_t next_start_position, bool is_last_call)
{
  while(!end_pq.empty() && (current_start_position < next_start_position || is_last_call))
  {
    int64_t top_end_pq = end_pq.top()->m_end_point;
    int64_t min_end_point = (is_last_call || (top_end_pq < (next_start_position - 1))) ? top_end_pq : (next_start_position-1);
    //Prepare gt_column
    gt_column->col_ = current_start_position;
    gt_column->reset();
    for(auto i=0ull;i<PQ_end_vec.size();++i)
    {
      auto& curr_struct = PQ_end_vec[i];
      if(curr_struct.m_needs_to_be_processed)
      {
	auto find_iter = tile_idx_2_iters.find(curr_struct.m_tile_idx);
	assert(find_iter != tile_idx_2_iters.end());
        gt_fill_row<StorageManager::const_iterator>(gt_column, i, curr_struct.m_array_column, curr_struct.m_cell_pos, 
	    &((*find_iter).second.m_iter_vector[0]));
      }
    }
    do_dummy_genotyping(gt_column, output_stream);
    //The following intervals have been completely processed
    while(!end_pq.empty() && end_pq.top()->m_end_point == min_end_point)
    {
      auto top_element = end_pq.top();
      top_element->m_needs_to_be_processed = false;
      auto find_iter = tile_idx_2_iters.find(top_element->m_tile_idx);
      assert(find_iter != tile_idx_2_iters.end());
      GTTileIteratorsTracker& current_iterator_tracker = find_iter->second;
      --(current_iterator_tracker.m_reference_counter);
      end_pq.pop();
    }
    current_start_position = min_end_point + 1;   //next start position, after the end
  }
}

void QueryProcessor::scan_and_operate(const StorageManager::ArrayDescriptor* ad, std::ostream& output_stream)
{
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
  unsigned int attribute_num = array_schema.attribute_num();
  
  unsigned num_queried_attributes = 7;
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  std::unordered_map<uint64_t, GTTileIteratorsTracker> tile_idx_2_iters;
  // END
  tile_its[0] = storage_manager_.begin(ad, 0);
  // REF
  tile_its[1] = storage_manager_.begin(ad, 1);
  // ALT
  tile_its[2] = storage_manager_.begin(ad, 2);
  // PL
  tile_its[3] = storage_manager_.begin(ad, 20);
  // NULL
  tile_its[4] = storage_manager_.begin(ad, 21);
  // OFFSETS
  tile_its[5] = storage_manager_.begin(ad, 22);
  // coordinates
  tile_its[6] = storage_manager_.begin(ad, attribute_num);
  
  StorageManager::const_iterator tile_it_end = storage_manager_.end(ad, attribute_num);
  
  QueryProcessor::GTColumn* gt_column = new GTColumn(0, total_num_samples);

  //Priority queue of END positions
  gVCF_PQ end_pq;
  //Vector of PQStruct - pre-allocate to eliminate allocations inside the while loop
  //Elements of the priority queue end_pq are pointers to elements of this vector
  auto PQ_end_vec = std::vector<PQStruct>(total_num_samples, PQStruct{});
  for(auto i=0ull;i<PQ_end_vec.size();++i)
    PQ_end_vec[i].m_sample_idx = i;
  //Current gVCF position being operated on
  int64_t current_start_position = -1ll;
  //Get first valid position in the array
  if(tile_its[num_queried_attributes-1] != tile_it_end)
  {
    Tile::const_iterator cell_it = (*(tile_its[num_queried_attributes-1])).begin();
    Tile::const_iterator cell_it_end = (*(tile_its[num_queried_attributes-1])).end();
    if(cell_it != cell_it_end)
    {
      std::vector<int64_t> current_coord = *cell_it;
      current_start_position = current_coord[1];
    }
  }
  int64_t next_start_position = -1ll;
  uint64_t tile_idx = 0ull;
  for(;tile_its[num_queried_attributes-1] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))
  {
    //Setup map for tile id to iterator
    auto insert_iter = tile_idx_2_iters.emplace(std::pair<uint64_t, GTTileIteratorsTracker>(tile_idx, GTTileIteratorsTracker(num_queried_attributes)));
    GTTileIteratorsTracker& current_iterator_tracker = insert_iter.first->second;
    for(auto i=0u;i<num_queried_attributes;++i)
      current_iterator_tracker.m_iter_vector[i] = tile_its[i];
    // Initialize cell iterators for the coordinates
    Tile::const_iterator cell_it = (*tile_its[num_queried_attributes-1]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[num_queried_attributes-1]).end();
    for(;cell_it != cell_it_end;++cell_it) {
      std::vector<int64_t> next_coord = *cell_it;
      if(next_coord[1] != current_start_position) //have found cell with next gVCF position, handle accumulated values
      {
        next_start_position = next_coord[1];
        assert(next_coord[1] > current_start_position);
        handle_gvcf_ranges(end_pq, PQ_end_vec, gt_column, tile_idx_2_iters, output_stream,
            current_start_position, next_start_position, false);
        assert(end_pq.empty() || end_pq.top()->m_end_point >= next_start_position);  //invariant
        current_start_position = next_start_position;
      }
      //Accumulate cells with position == current_start_position
      uint64_t sample_idx = next_coord[0];
      auto& curr_struct = PQ_end_vec[sample_idx];
      //Store array column idx corresponding to this cell
      curr_struct.m_array_column = current_start_position;
      //Get END corresponding to this cell
      curr_struct.m_end_point = static_cast<const AttributeTile<int64_t>& >(*tile_its[0]).cell(cell_it.pos());
      assert(curr_struct.m_end_point >= current_start_position);
      //Store tile idx
      curr_struct.m_tile_idx = tile_idx;
      ++(current_iterator_tracker.m_reference_counter);
      //Store position of cell wrt tile
      curr_struct.m_cell_pos = cell_it.pos();
      assert(cell_it.pos() < (*tile_its[0]).cell_num());
      curr_struct.m_needs_to_be_processed = true;
      end_pq.push(&(PQ_end_vec[sample_idx]));
      assert(end_pq.size() <= total_num_samples);
    }
    for(auto map_iter = tile_idx_2_iters.begin(), map_end = tile_idx_2_iters.end();map_iter != map_end;)
    {
      if(map_iter->second.m_reference_counter == 0ull)
      {
	auto tmp_iter = map_iter;
	map_iter++;
	tile_idx_2_iters.erase(tmp_iter);
      }
      else
	map_iter++;
    }
    ++tile_idx;
  }
  handle_gvcf_ranges(end_pq, PQ_end_vec, gt_column, tile_idx_2_iters, output_stream,
      current_start_position, 0, true);
  delete[] tile_its;
  delete gt_column;
}

QueryProcessor::GTColumn* QueryProcessor::gt_get_column(
    const StorageManager::ArrayDescriptor* ad, uint64_t col, QueryProcessor::GTProfileStats* stats) const {
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t row_num = dim_domains[0].second - dim_domains[0].first + 1;
  unsigned int attribute_num = array_schema.attribute_num();

  // Check that column falls into the domain of the second dimension
  assert(col >= dim_domains[1].first && col <= dim_domains[1].second);

  // Indicates how many rows have been filled.
  uint64_t filled_rows = 0;

  // Initialize reverse tile iterators for 
  // END, REF, ALT, PL, NULL, OFFSETS, coordinates
  // The reverse tile iterator will start with the tiles
  // of the various attributes that have the largest
  // id that either intersect with col, or precede col.
  StorageManager::const_reverse_iterator* tile_its;
  StorageManager::const_reverse_iterator tile_it_end;  
  unsigned int gt_attribute_num = 
      gt_initialize_tile_its(ad, tile_its, tile_it_end, col);

  // Create and initialize GenotypingColumn members
  QueryProcessor::GTColumn* gt_column = new GTColumn(col, row_num);

  // Create cell iterators
  Tile::const_reverse_iterator cell_it, cell_it_end;
#ifdef DO_PROFILING
  uint64_t num_cells_touched = 0;
  uint64_t num_tiles_touched = 0;
  bool first_sample = true;
#endif
  // Fill the genotyping column
  while(tile_its[gt_attribute_num] != tile_it_end && filled_rows < row_num) {
    // Initialize cell iterators for the coordinates
    cell_it = (*tile_its[gt_attribute_num]).rbegin();
    cell_it_end = (*tile_its[gt_attribute_num]).rend();
    while(cell_it != cell_it_end && filled_rows < row_num) {
      std::vector<int64_t> next_coord = *cell_it;
#ifdef DO_PROFILING
      ++num_cells_touched;
#endif
      // If next cell is not on the right of col, and corresponds to 
      // uninvestigated row
      if(next_coord[1] <= col && CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(gt_column->REF_[next_coord[0]])) {
        gt_fill_row<StorageManager::const_reverse_iterator>(gt_column, next_coord[0], next_coord[1], cell_it.pos(), tile_its);
        ++filled_rows;
#ifdef DO_PROFILING
	if(first_sample)
	{
	  stats->m_sum_num_cells_first_sample += num_cells_touched;
	  stats->m_sum_sq_num_cells_first_sample += (num_cells_touched*num_cells_touched);
	  first_sample = false;
	}
	else
	{
	  ++(stats->m_num_samples);
	  stats->m_sum_num_cells_touched += num_cells_touched;
	  stats->m_sum_sq_num_cells_touched += (num_cells_touched*num_cells_touched);
	}
	num_cells_touched = 0;
#endif
      }
      ++cell_it;
    }
    advance_tile_its(gt_attribute_num, tile_its);
#ifdef DO_PROFILING
    ++num_tiles_touched;
#endif
  }

  //No need for this assertion
  //assert(filled_rows == row_num);

  delete [] tile_its;

#ifdef DO_PROFILING
  if(num_cells_touched > 0)	//last iteration till invalid
  {
    stats->m_sum_num_cells_last_iter += num_cells_touched;
    stats->m_sum_sq_num_cells_last_iter += (num_cells_touched*num_cells_touched);
    num_cells_touched = 0;
  }
  stats->m_sum_num_tiles_touched += num_tiles_touched;
  stats->m_sum_sq_num_tiles_touched += (num_tiles_touched*num_tiles_touched);
#endif

  return gt_column;
}

void QueryProcessor::join(const StorageManager::ArrayDescriptor* ad_A, 
                          const StorageManager::ArrayDescriptor* ad_B,
                          const std::string& result_array_name) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_info()->array_schema_;
  const ArraySchema& array_schema_B = ad_B->array_info()->array_schema_;

  std::pair<bool,std::string> can_join = 
      ArraySchema::join_compatible(array_schema_A, array_schema_B);

  if(!can_join.first)
    throw QueryProcessorException(std::string("[QueryProcessor] Input arrays "
                                  " are not join-compatible.") + 
                                  can_join.second);

  ArraySchema array_schema_C = ArraySchema::create_join_result_schema(
                                   array_schema_A, 
                                   array_schema_B, 
                                   result_array_name);
  
  if(array_schema_A.has_regular_tiles())
    join_regular(ad_A, ad_B, array_schema_C);
  else 
    join_irregular(ad_A, ad_B, array_schema_C);
} 

void QueryProcessor::subarray(const StorageManager::ArrayDescriptor* ad,
                              const Tile::Range& range,
                              const std::string& result_array_name) const { 
  if(ad->array_schema().has_regular_tiles())
    subarray_regular(ad, range, result_array_name);
  else 
    subarray_irregular(ad, range, result_array_name);
}

/******************************************************
******************* PRIVATE METHODS *******************
******************************************************/

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_iterator* cell_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
      ++cell_its[i];
}

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_iterator* cell_its,
                                      int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    cell_its[i] += step;
}

inline
void QueryProcessor::advance_cell_its(
    unsigned int attribute_num, Tile::const_reverse_iterator* cell_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
      ++cell_its[i];
}

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_reverse_iterator* cell_its,
                                      int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    cell_its[i] += step;
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_iterator* tile_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    ++tile_its[i];
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_iterator* tile_its, 
    int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    tile_its[i] += step;
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_reverse_iterator* tile_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    ++tile_its[i];
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_reverse_iterator* tile_its, 
    int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    tile_its[i] += step;
}

inline
void QueryProcessor::append_cell(const Tile::const_iterator* cell_its,
                                 Tile** tiles,
                                 unsigned int attribute_num) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    *tiles[i] << cell_its[i]; 
}

inline
void QueryProcessor::append_cell(const Tile::const_iterator* cell_its_A,
                                 const Tile::const_iterator* cell_its_B,
                                 Tile** tiles_C,
                                 unsigned int attribute_num_A,
                                 unsigned int attribute_num_B) const {
  for(unsigned int i=0; i<attribute_num_A; i++) 
    *tiles_C[i] << cell_its_A[i]; 
  for(unsigned int i=0; i<=attribute_num_B; i++)
    *tiles_C[attribute_num_A+i] << cell_its_B[i]; 
}

inline
CSVLine QueryProcessor::cell_to_csv_line(const Tile::const_iterator* cell_its,
                                         unsigned int attribute_num) const {
  CSVLine csv_line;

  // Append coordinates first
  cell_its[attribute_num] >> csv_line;
  // Append attribute values next
  for(unsigned int i=0; i<attribute_num; i++)
    cell_its[i] >> csv_line;

  return csv_line;
}

void QueryProcessor::create_workspace() const {
  struct stat st;
  stat(workspace_.c_str(), &st);

  // If the workspace does not exist, create it
  if(!S_ISDIR(st.st_mode)) { 
    int dir_flag = mkdir(workspace_.c_str(), S_IRWXU);
    assert(dir_flag == 0);
  }
}

inline
void QueryProcessor::get_tiles(
    const StorageManager::ArrayDescriptor* ad, 
    uint64_t tile_id, const Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Get attribute tiles
  for(unsigned int i=0; i<=attribute_num; i++) 
   tiles[i] = storage_manager_.get_tile(ad, i, tile_id);
}

bool QueryProcessor::path_exists(const std::string& path) const {
  struct stat st;
  stat(path.c_str(), &st);
  return S_ISDIR(st.st_mode);
}

template<class ITER>
void QueryProcessor::gt_fill_row(
    GTColumn* gt_column, int64_t row, int64_t column, int64_t pos,
    const ITER* tile_its) const {
  // First check if the row is NULL
  int64_t END_v = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[0]).cell(pos);
  if(END_v < gt_column->col_) {
    gt_column->REF_[row] = "$";
    return;
  }

  // Retrieve the offsets
  const AttributeTile<int64_t>& OFFSETS_tile = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[5]);
  int64_t REF_offset = OFFSETS_tile.cell(pos*5);
  int64_t ALT_offset = OFFSETS_tile.cell(pos*5+1);
  int64_t PL_offset = OFFSETS_tile.cell(pos*5+4);

  // Retrieve the NULL bitmap
  const AttributeTile<int>& NULL_tile = 
      static_cast<const AttributeTile<int>& >(*tile_its[4]);
  int NULL_bitmap = NULL_tile.cell(pos);

  char c;
  int i;

  // Fill the REF
  //If the queried column is identical to the cell's column, then the REF value stored is correct
  //Else, the cell stores an interval and the REF value is set to "N" which means could be anything
  if(column == gt_column->col_)
  {
    const AttributeTile<char>& REF_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[1]);
    std::string REF_s = "";
    i = 0;
    while((c = REF_tile.cell(REF_offset+i)) != '\0') { 
      REF_s.push_back(c);
      ++i;
    }
    gt_column->REF_[row] = REF_s;
  }
  else
    gt_column->REF_[row] = "N";

  // Fill the ALT values
  const AttributeTile<char>& ALT_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[2]);
  i = 0;
  std::string ALT_s = "";
  while((c = ALT_tile.cell(ALT_offset+i)) != '&') {
    if(c == '\0') {
      gt_column->ALT_[row].push_back(ALT_s);
      ALT_s = "";
    } else {
      ALT_s.push_back(c);
    }
    i++;
  }
  assert(ALT_s == "");
  gt_column->ALT_[row].push_back("&");

  // Fill the PL values
  if((NULL_bitmap & 1) == 0) { // If the PL values are
    const AttributeTile<int>& PL_tile = 
        static_cast<const AttributeTile<int>& >(*tile_its[3]);
    int ALT_num = gt_column->ALT_[row].size(); 
    int PL_num = (ALT_num+1)*(ALT_num+2)/2;
    for(int i=0; i<PL_num; i++) 
      gt_column->PL_[row].push_back(PL_tile.cell(PL_offset+i));
  }
}

inline
unsigned int QueryProcessor::gt_initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    StorageManager::const_reverse_iterator*& tile_its, 
    StorageManager::const_reverse_iterator& tile_it_end,
    uint64_t col) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Create reverse iterators
  tile_its = new StorageManager::const_reverse_iterator[7];
  // Find the rank of the tile the left sweep starts from.
  uint64_t start_rank = storage_manager_.get_left_sweep_start_rank(ad, col);

  // END
  tile_its[0] = storage_manager_.rbegin(ad, 0, start_rank);
  // REF
  tile_its[1] = storage_manager_.rbegin(ad, 1, start_rank);
  // ALT
  tile_its[2] = storage_manager_.rbegin(ad, 2, start_rank);
  // PL
  tile_its[3] = storage_manager_.rbegin(ad, 20, start_rank);
  // NULL
  tile_its[4] = storage_manager_.rbegin(ad, 21, start_rank);
  // OFFSETS
  tile_its[5] = storage_manager_.rbegin(ad, 22, start_rank);
  // coordinates
  tile_its[6] = storage_manager_.rbegin(ad, attribute_num, start_rank);
  tile_it_end = storage_manager_.rend(ad, attribute_num);

  // The number of attributes is 6, and the coordinates is the extra one
  return 6;
}

inline
void QueryProcessor::initialize_cell_its(
    const Tile** tiles, unsigned int attribute_num,
    Tile::const_iterator* cell_its, Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = tiles[i]->begin();
  cell_it_end = tiles[attribute_num]->end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_iterator* tile_its, unsigned int attribute_num,
    Tile::const_iterator* cell_its, Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).begin();
  cell_it_end = (*tile_its[attribute_num]).end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_reverse_iterator* tile_its, 
    unsigned int attribute_num,
    Tile::const_iterator* cell_its, 
    Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).begin();
  cell_it_end = (*tile_its[attribute_num]).end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_reverse_iterator* tile_its, 
    unsigned int attribute_num,
    Tile::const_reverse_iterator* cell_its, 
    Tile::const_reverse_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).rbegin();
  cell_it_end = (*tile_its[attribute_num]).rend();
}

inline
void QueryProcessor::initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    StorageManager::const_iterator* tile_its, 
    StorageManager::const_iterator& tile_it_end) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  for(unsigned int i=0; i<=attribute_num; i++)
    tile_its[i] = storage_manager_.begin(ad, i);
  tile_it_end = storage_manager_.end(ad, attribute_num);
}

void QueryProcessor::join_irregular(const StorageManager::ArrayDescriptor* ad_A, 
                                    const StorageManager::ArrayDescriptor* ad_B,
                                    const ArraySchema& array_schema_C) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_schema();
  const ArraySchema& array_schema_B = ad_B->array_schema();
  unsigned int attribute_num_A = array_schema_A.attribute_num();
  unsigned int attribute_num_B = array_schema_B.attribute_num();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Prepare result array
  const StorageManager::ArrayDescriptor* ad_C = 
      storage_manager_.open_array(array_schema_C);

  // Create tiles 
  const Tile** tiles_A = new const Tile*[attribute_num_A+1];
  const Tile** tiles_B = new const Tile*[attribute_num_B+1];
  Tile** tiles_C = new Tile*[attribute_num_C+1];
  
  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its_A = 
      new StorageManager::const_iterator[attribute_num_A+1];
  StorageManager::const_iterator *tile_its_B = 
      new StorageManager::const_iterator[attribute_num_B+1];
  StorageManager::const_iterator tile_it_end_A;
  StorageManager::const_iterator tile_it_end_B;
  initialize_tile_its(ad_A, tile_its_A, tile_it_end_A);
  initialize_tile_its(ad_B, tile_its_B, tile_it_end_B);
  
  // Create cell iterators
  Tile::const_iterator* cell_its_A = 
      new Tile::const_iterator[attribute_num_A+1];
  Tile::const_iterator* cell_its_B = 
      new Tile::const_iterator[attribute_num_B+1];
  Tile::const_iterator cell_it_end_A, cell_it_end_B;

  // Auxiliary variables storing the number of skipped tiles when joining.
  // It is used to advance only the coordinates iterator when a tile is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a tile joins.
  int64_t skipped_tiles_A = 0;
  int64_t skipped_tiles_B = 0;

  // Initialize tiles with id 0 for C (result array)
  new_tiles(array_schema_C, 0, tiles_C); 

  // Join algorithm
  while(tile_its_A[attribute_num_A] != tile_it_end_A &&
        tile_its_B[attribute_num_B] != tile_it_end_B) {
    // Potential join result generation
    if(may_join(tile_its_A[attribute_num_A], tile_its_B[attribute_num_B])) {
      // Update iterators in A
      if(skipped_tiles_A) {
        advance_tile_its(attribute_num_A, tile_its_A, skipped_tiles_A);
        skipped_tiles_A = 0;
        initialize_cell_its(tile_its_A, attribute_num_A, 
                            cell_its_A, cell_it_end_A);
      }
      // Update iterators in B
      if(skipped_tiles_B) {
        advance_tile_its(attribute_num_B, tile_its_B, skipped_tiles_B);
        skipped_tiles_B = 0;
        initialize_cell_its(tile_its_B, attribute_num_B, 
                            cell_its_B, cell_it_end_B);
      }
      // Join the tiles
      join_tiles_irregular(attribute_num_A, cell_its_A, cell_it_end_A, 
                           attribute_num_B, cell_its_B, cell_it_end_B,
                           ad_C, tiles_C);
    }

    // Check which tile precedes the other in the global order
    // Note that operator '<', when the operands are from different
    // arrays, returns true if the first tile precedes the second
    // in the global order by checking their bounding coordinates.
    if(tile_its_A[attribute_num_A] < tile_its_B[attribute_num_B]) {
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
    }
    else {
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    }
  }
  
  // Send the lastly created tiles to storage manager
  store_tiles(ad_C, tiles_C);

  // Close result array
  storage_manager_.close_array(ad_C);

  // Clean up
  delete [] tiles_A;
  delete [] tiles_B;
  delete [] tiles_C;
  delete [] tile_its_A;
  delete [] tile_its_B;
  delete [] cell_its_A;
  delete [] cell_its_B;
}

void QueryProcessor::join_regular(const StorageManager::ArrayDescriptor* ad_A, 
                                  const StorageManager::ArrayDescriptor* ad_B,
                                  const ArraySchema& array_schema_C) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_schema();
  const ArraySchema& array_schema_B = ad_B->array_schema();
  unsigned int attribute_num_A = array_schema_A.attribute_num();
  unsigned int attribute_num_B = array_schema_B.attribute_num();
  unsigned int attribute_num_C = array_schema_C.attribute_num();
  uint64_t tile_id_A, tile_id_B;

  // Prepare result array
  const StorageManager::ArrayDescriptor* ad_C = 
      storage_manager_.open_array(array_schema_C);

  // Create tiles 
  const Tile** tiles_A = new const Tile*[attribute_num_A+1];
  const Tile** tiles_B = new const Tile*[attribute_num_B+1];
  Tile** tiles_C = new Tile*[attribute_num_C+1];
  
  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its_A = 
      new StorageManager::const_iterator[attribute_num_A+1];
  StorageManager::const_iterator *tile_its_B = 
      new StorageManager::const_iterator[attribute_num_B+1];
  StorageManager::const_iterator tile_it_end_A;
  StorageManager::const_iterator tile_it_end_B;
  initialize_tile_its(ad_A, tile_its_A, tile_it_end_A);
  initialize_tile_its(ad_B, tile_its_B, tile_it_end_B);
  
  // Create cell iterators
  Tile::const_iterator* cell_its_A = 
      new Tile::const_iterator[attribute_num_A+1];
  Tile::const_iterator* cell_its_B = 
      new Tile::const_iterator[attribute_num_B+1];
  Tile::const_iterator cell_it_end_A, cell_it_end_B;

  // Auxiliary variables storing the number of skipped tiles when joining.
  // It is used to advance only the coordinates iterator when a tile is
  // finished/skipped, and then efficiently advance the attribute iterators only
  // when a tile joins.
  int64_t skipped_tiles_A = 0;
  int64_t skipped_tiles_B = 0;

  // Join algorithm
  while(tile_its_A[attribute_num_A] != tile_it_end_A &&
        tile_its_B[attribute_num_B] != tile_it_end_B) {
    tile_id_A = tile_its_A[attribute_num_A].tile_id();
    tile_id_B = tile_its_B[attribute_num_B].tile_id();

    // Potential join result generation
    if(tile_id_A == tile_id_B) {
      // Update iterators in A
      if(skipped_tiles_A) {
        advance_tile_its(attribute_num_A, tile_its_A, skipped_tiles_A);
        skipped_tiles_A = 0;
        initialize_cell_its(tile_its_A, attribute_num_A, 
                            cell_its_A, cell_it_end_A);
      }
      // Update iterators in B
      if(skipped_tiles_B) {
        advance_tile_its(attribute_num_B, tile_its_B, skipped_tiles_B);
        skipped_tiles_B = 0;
        initialize_cell_its(tile_its_B, attribute_num_B, 
                            cell_its_B, cell_it_end_B);
      }

      // Initialize tiles for C (result array)
      new_tiles(array_schema_C, tile_id_A, tiles_C);
      // Join the tiles
      join_tiles_regular(attribute_num_A, cell_its_A, cell_it_end_A, 
                         attribute_num_B, cell_its_B, cell_it_end_B,
                         ad_C, tiles_C);
      // Send the created tiles to storage manager
      store_tiles(ad_C, tiles_C);
    }

    // Tile precedence in the case of regular tiles is simply determined
    // by the order of the tile ids.
    if(tile_id_A < tile_id_B) {
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
    } else if(tile_id_A > tile_id_B) {
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    } else { // tile_id_A == tile_id_B, advance both
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    }
  }
  
  // Close result array
  storage_manager_.close_array(ad_C);

  // Clean up
  delete [] tiles_A;
  delete [] tiles_B;
  delete [] tiles_C;
  delete [] tile_its_A;
  delete [] tile_its_B;
  delete [] cell_its_A;
  delete [] cell_its_B;
}

void QueryProcessor::join_tiles_irregular(
    unsigned int attribute_num_A, Tile::const_iterator* cell_its_A,
    Tile::const_iterator& cell_it_end_A, 
    unsigned int attribute_num_B, Tile::const_iterator* cell_its_B,
    Tile::const_iterator& cell_it_end_B,
    const StorageManager::ArrayDescriptor* ad_C, Tile** tiles_C) const {
  // For easy reference
  const ArraySchema& array_schema_C = ad_C->array_schema();
  uint64_t capacity = array_schema_C.capacity();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Auxiliary variables storing the number of skipped cells when joining.
  // It is used to advance only the coordinates iterator when a cell is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a cell joins.
  int64_t skipped_cells_A = 0;
  int64_t skipped_cells_B = 0;

  while(cell_its_A[attribute_num_A] != cell_it_end_A &&
        cell_its_B[attribute_num_B] != cell_it_end_B) {
    // If the coordinates are equal
    // Note that operator '==', when the operands correspond to different
    // tiles, returns true if the cell values pointed by the iterators
    // are equal.
    if(cell_its_A[attribute_num_A] == cell_its_B[attribute_num_B]) {      
      advance_cell_its(attribute_num_A, cell_its_A, skipped_cells_A);
      advance_cell_its(attribute_num_B, cell_its_B, skipped_cells_B);
      skipped_cells_A = 0;
      skipped_cells_B = 0;
      if(tiles_C[attribute_num_C]->cell_num() == capacity) {
        uint64_t new_tile_id = tiles_C[attribute_num_C]->tile_id() + 1;
        store_tiles(ad_C, tiles_C);
        new_tiles(array_schema_C, new_tile_id, tiles_C); 
      }
      append_cell(cell_its_A, cell_its_B, tiles_C, 
                  attribute_num_A, attribute_num_B);
      advance_cell_its(attribute_num_A, cell_its_A);
      advance_cell_its(attribute_num_B, cell_its_B);
    // Otherwise check which cell iterator to advance
    } else {
      if(array_schema_C.precedes(cell_its_A[attribute_num_A],
                                 cell_its_B[attribute_num_B])) {
        ++cell_its_A[attribute_num_A];
        ++skipped_cells_A;
      } else {
        ++cell_its_B[attribute_num_B];
        ++skipped_cells_B;
      }
    }
  }
}

void QueryProcessor::join_tiles_regular(
    unsigned int attribute_num_A, Tile::const_iterator* cell_its_A,
    Tile::const_iterator& cell_it_end_A, 
    unsigned int attribute_num_B, Tile::const_iterator* cell_its_B,
    Tile::const_iterator& cell_it_end_B,
    const StorageManager::ArrayDescriptor* ad_C, Tile** tiles_C) const {
  // For easy reference
  const ArraySchema& array_schema_C = ad_C->array_schema();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Auxiliary variables storing the number of skipped cells when joining.
  // It is used to advance only the coordinates iterator when a cell is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a cell joins.
  int64_t skipped_cells_A = 0;
  int64_t skipped_cells_B = 0;

  while(cell_its_A[attribute_num_A] != cell_it_end_A &&
        cell_its_B[attribute_num_B] != cell_it_end_B) {
    // If the coordinates are equal
    // Note that operator '==', when the operands correspond to different
    // tiles, returns true if the cell values pointed by the iterators
    // are equal.
    if(cell_its_A[attribute_num_A] == cell_its_B[attribute_num_B]) {      
      advance_cell_its(attribute_num_A, cell_its_A, skipped_cells_A);
      advance_cell_its(attribute_num_B, cell_its_B, skipped_cells_B);
      skipped_cells_A = 0;
      skipped_cells_B = 0;
      append_cell(cell_its_A, cell_its_B, tiles_C, 
                  attribute_num_A, attribute_num_B);
      advance_cell_its(attribute_num_A, cell_its_A);
      advance_cell_its(attribute_num_B, cell_its_B);
    // Otherwise check which cell iterator to advance
    } else {
      if(array_schema_C.precedes(cell_its_A[attribute_num_A],
                                 cell_its_B[attribute_num_B])) {
        ++cell_its_A[attribute_num_A];
        ++skipped_cells_A;
      } else {
        ++cell_its_B[attribute_num_B];
        ++skipped_cells_B;
      }
    }
  }
}

bool QueryProcessor::may_join(
    const StorageManager::const_iterator& it_A,
    const StorageManager::const_iterator& it_B) const {
  // For easy reference
  const ArraySchema& array_schema_A = it_A.array_schema();
  const MBR& mbr_A = it_A.mbr();
  const MBR& mbr_B = it_B.mbr();

  // Check if the tile MBRs overlap
  if(array_schema_A.has_irregular_tiles() && !overlap(mbr_A, mbr_B))
    return false;

  // For easy reference 
  BoundingCoordinatesPair bounding_coordinates_A = it_A.bounding_coordinates();
  BoundingCoordinatesPair bounding_coordinates_B = it_B.bounding_coordinates();

  // Check if the cell id ranges (along the global order) intersect
  if(!overlap(bounding_coordinates_A, bounding_coordinates_B, array_schema_A))
    return false;

  return true;
}

inline
void QueryProcessor::new_tiles(const ArraySchema& array_schema, 
                               uint64_t tile_id, Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = array_schema.attribute_num();
  uint64_t capacity = array_schema.capacity();

  for(unsigned int i=0; i<=attribute_num; i++)
    tiles[i] = storage_manager_.new_tile(array_schema, i, tile_id, capacity);
}

bool QueryProcessor::overlap(const MBR& mbr_A, const MBR& mbr_B) const {
  assert(mbr_A.size() == mbr_B.size());
  assert(mbr_A.size() % 2 == 0);

  // For easy rederence
  unsigned int dim_num = mbr_A.size() / 2;

  for(unsigned int i=0; i<dim_num; i++) 
    if(mbr_A[2*i+1] < mbr_B[2*i] || mbr_A[2*i] > mbr_B[2*i+1])
      return false;

  return true;
}

bool QueryProcessor::overlap(
    const BoundingCoordinatesPair& bounding_coordinates_A,
    const BoundingCoordinatesPair& bounding_coordinates_B, 
    const ArraySchema& array_schema) const {
  if(array_schema.precedes(bounding_coordinates_A.second, 
                           bounding_coordinates_B.first) ||
     array_schema.succeeds(bounding_coordinates_A.first, 
                           bounding_coordinates_B.second))
    return false;
  else
    return true;
}

inline
void QueryProcessor::set_workspace(const std::string& path) {
  workspace_ = path;
  
  // Replace '~' with the absolute path
  if(workspace_[0] == '~') {
    workspace_ = std::string(getenv("HOME")) +
                 workspace_.substr(1, workspace_.size()-1);
  }

  // Check if the input path is an existing directory 
  assert(path_exists(workspace_));
 
  workspace_ += "/QueryProcessor";
}

inline
void QueryProcessor::store_tiles(const StorageManager::ArrayDescriptor* ad,
                                 Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Append attribute tiles
  for(unsigned int i=0; i<=attribute_num; i++)
    storage_manager_.append_tile(tiles[i], ad, i);
} 

void QueryProcessor::subarray_irregular(
    const StorageManager::ArrayDescriptor* ad,
    const Tile::Range& range, const std::string& result_array_name) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  unsigned int attribute_num = array_schema.attribute_num();
  uint64_t capacity = array_schema.capacity();

  // Create tiles
  const Tile** tiles = new const Tile*[attribute_num+1];
  Tile** result_tiles = new Tile*[attribute_num+1];

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Prepare result array
  ArraySchema result_array_schema = array_schema.clone(result_array_name);
  const StorageManager::ArrayDescriptor* result_ad = 
      storage_manager_.open_array(result_array_schema);
  
  // Get the tile ids that overlap with the range
  std::vector<std::pair<uint64_t, bool> > overlapping_tile_ids;
  storage_manager_.get_overlapping_tile_ids(ad, range, &overlapping_tile_ids);

  // Initialize tile iterators
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it = 
      overlapping_tile_ids.begin();
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it_end =
      overlapping_tile_ids.end();
      
  // Create result tiles and load input array tiles 
  uint64_t tile_id = 0;
  new_tiles(result_array_schema, tile_id, result_tiles); 

  // Auxiliary variable storing the number of skipped cells when investigating
  // a tile partially overlapping the range. It is used to advance only the
  // coordinates iterator when a cell is not in range, and then efficiently
  // advance the attribute iterators only when a cell falls in the range.
  int64_t skipped;

  // Iterate over all tiles
  for(; tile_id_it != tile_id_it_end; ++tile_id_it) {
    get_tiles(ad, tile_id_it->first, tiles);
    initialize_cell_its(tiles, attribute_num, cell_its, cell_it_end); 
    skipped = 0;

    if(tile_id_it->second) { // Full overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(result_tiles[attribute_num]->cell_num() == capacity) {
          store_tiles(result_ad, result_tiles);
          new_tiles(result_array_schema, ++tile_id, result_tiles); 
        }
        append_cell(cell_its, result_tiles, attribute_num);
        advance_cell_its(attribute_num, cell_its);
      }
    } else { // Partial overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(cell_its[attribute_num].cell_inside_range(range)) {
          if(result_tiles[attribute_num]->cell_num() == capacity) {
            store_tiles(result_ad, result_tiles);
            new_tiles(result_array_schema, ++tile_id, result_tiles); 
          }
          advance_cell_its(attribute_num, cell_its, skipped);
          skipped = 0;
          append_cell(cell_its, result_tiles, attribute_num);
          advance_cell_its(attribute_num, cell_its);
        } else { // Advance only the coordinates cell iterator
          skipped++;
          ++cell_its[attribute_num];
        }
      }
    }
  } 

  // Send the lastly created tiles to storage manager
  store_tiles(result_ad, result_tiles);
  
  // Close result array
  storage_manager_.close_array(result_ad);

  // Clean up 
  delete [] tiles;
  delete [] result_tiles;
  delete [] cell_its;
}

void QueryProcessor::subarray_regular(
    const StorageManager::ArrayDescriptor* ad,
    const Tile::Range& range, const std::string& result_array_name) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  unsigned int attribute_num = array_schema.attribute_num();

  // Create tiles 
  const Tile** tiles = new const Tile*[attribute_num+1];
  Tile** result_tiles = new Tile*[attribute_num+1];

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Prepare result array
  ArraySchema result_array_schema = array_schema.clone(result_array_name);
  const StorageManager::ArrayDescriptor* result_ad = 
      storage_manager_.open_array(result_array_schema);
    
  // Get the tile ids that overlap with the range
  std::vector<std::pair<uint64_t, bool> > overlapping_tile_ids;
  storage_manager_.get_overlapping_tile_ids(ad, range, &overlapping_tile_ids);

  // Initialize tile iterators
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it = 
      overlapping_tile_ids.begin();
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it_end =
      overlapping_tile_ids.end();


  // Auxiliary variable storing the number of skipped cells when investigating
  // a tile partially overlapping the range. It is used to advance only the
  // coordinates iterator when a cell is not in range, and then efficiently
  // advance the attribute iterators only when a cell falls in the range.
  int64_t skipped;

  // Iterate over all overlapping tiles
  for(; tile_id_it != tile_id_it_end; ++tile_id_it) {
    // Create result tiles and load input array tiles 
    new_tiles(result_array_schema, tile_id_it->first, result_tiles); 
    get_tiles(ad, tile_id_it->first, tiles); 
    initialize_cell_its(tiles, attribute_num, cell_its, cell_it_end); 
    skipped = 0;
 
    if(tile_id_it->second) { // Full overlap
      while(cell_its[attribute_num] != cell_it_end) {
        append_cell(cell_its, result_tiles, attribute_num);
        advance_cell_its(attribute_num, cell_its);
      }
    } else { // Partial overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(cell_its[attribute_num].cell_inside_range(range)) {
          advance_cell_its(attribute_num, cell_its, skipped);
          skipped = 0;
          append_cell(cell_its, result_tiles, attribute_num);
          advance_cell_its(attribute_num, cell_its);
        } else { // Advance only the coordinates cell iterator
          ++skipped;
          ++cell_its[attribute_num];
        }
      }
    }
      
    // Send new tiles to storage manager
    store_tiles(result_ad, result_tiles);
  } 
  
  // Close result array
  storage_manager_.close_array(result_ad);

  // Clean up 
  delete [] tiles;
  delete [] result_tiles;
  delete [] cell_its;
}

