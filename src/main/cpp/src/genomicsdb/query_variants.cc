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

#include "gt_common.h"
#include "query_variants.h"
#include "timer.h"

using namespace std;

#if 0
//Utility functions
//Search for cell with given search value
/*
 * Given a column idx and the co-ordinate tile, returns the iterator in the tile pointing
 * to the first cell with the next valid column (>column)
 */
Tile::const_iterator get_first_cell_after(const Tile& tile, uint64_t column)
{
  uint64_t search_value = column+1;
  //Invariant: low always points to cell with value < search_value, high to cell with value>=search_value
  //Hence, low and high initially "point" outside the tile
  int64_t low = -1ll;
  int64_t high = tile.cell_num();
  int64_t mid = 0;
  while(high > low + 1)
  {
    mid = (high + low)/2;
    auto mid_value = static_cast<const CoordinateTile<int64_t>&>(tile).cell(mid)[1];
    if(mid_value >= search_value)
      high = mid;
    else
      low = mid;
#ifdef DEBUG
    //Check invariant
    if(high < tile.cell_num())    //not end of tile
    {
      auto value = static_cast<const CoordinateTile<int64_t>&>(tile).cell(high)[1];
      assert(value >= search_value);
      if(high > 0)
      {
        //previous cell should be less than search position
        auto value = static_cast<const CoordinateTile<int64_t>&>(tile).cell(high-1)[1];
        assert(value < search_value);
      }
    }
#endif
  }
  return Tile::const_iterator(&tile, high);
}
#endif

//Profiling functions
GTProfileStats::GTProfileStats()
{
  m_num_queries = 0;
  m_stats_sum_vector.resize(GTProfileStats::GT_NUM_STATS);
  m_stats_sum_sq_vector.resize(GTProfileStats::GT_NUM_STATS);
  m_stats_tmp_count_vector.resize(GTProfileStats::GT_NUM_STATS);
  for(auto i=0u;i<m_stats_sum_sq_vector.size();++i)
    m_stats_sum_sq_vector[i] = m_stats_sum_vector[i] = m_stats_tmp_count_vector[i] = 0;
  m_stats_name_vector = std::vector<std::string>{
    "GT_NUM_CELLS",   //total #cells traversed in the query, coord cells are accessed for every time
      "GT_NUM_CELLS_IN_LEFT_SWEEP",//total #cells traversed in the query left sweep, coord cells are accessed for every time
      "GT_NUM_VALID_CELLS_IN_QUERY",//#valid cells actually returned in query 
      "GT_NUM_ATTR_CELLS_ACCESSED",//#attribute cells accessed in the query
      "GT_NUM_PQ_FLUSHES_DUE_TO_OVERLAPPING_CELLS",//#times PQ gets flushed due to overlapping cells in the input
      "GT_NUM_OPERATOR_INVOCATIONS" //#times operator gets invoked
  };
}

void GTProfileStats::print_stats(std::ostream& fptr) const
{
  fptr << "stat_name,sum,sum_sq,mean,std-dev\n";
  for(auto i=0u;i<GTProfileStats::GT_NUM_STATS;++i)
  {
    fptr << m_stats_name_vector[i];
    fptr << "," << m_stats_sum_vector[i] << "," << std::setprecision(6) << m_stats_sum_sq_vector[i];
    if(m_num_queries == 0)
      fptr << ",*,*";
    else
    {
      double mean = ((double)m_stats_sum_vector[i])/m_num_queries;
      double std_dev = sqrt(abs((m_stats_sum_sq_vector[i]/m_num_queries) - (mean*mean)));
      fptr << "," << std::setprecision(6) << mean << "," << std_dev;
    }
    fptr << "\n";
  }
  m_interval_sweep_timer.print("Sweep at query begin position", std::cerr);
  m_operator_timer.print("Operator time", std::cerr);
  m_genomicsdb_cell_fill_timer.print("GenomicsDB cell fill timer", std::cerr);
}

//Static members
bool VariantQueryProcessor::m_are_static_members_initialized = false;
unordered_map<type_index, shared_ptr<VariantFieldCreatorBase>> VariantQueryProcessor::m_type_index_to_creator;

//Initialize static members function
void VariantQueryProcessor::initialize_static_members()
{
  VariantQueryProcessor::m_type_index_to_creator.clear();
  //Map type_index to creator functions
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(bool))] =
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<uint8_t, unsigned>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(int8_t))] =
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<int8_t, int>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(uint8_t))] =
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<uint8_t, unsigned>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(int))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<int>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(unsigned))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<unsigned>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(int64_t))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<int64_t>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(uint64_t))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<uint64_t>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(float))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<float>>());
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(double))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldPrimitiveVectorData<double>>());
  //Char becomes string instead of vector<char>
  VariantQueryProcessor::m_type_index_to_creator[std::type_index(typeid(char))] = 
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldData<std::string>>()); 
  //Set initialized flag
  VariantQueryProcessor::m_are_static_members_initialized = true;
}

void VariantQueryProcessor::initialize_known(const VariantArraySchema& schema)
{
  //Initialize schema idx <--> known_field_enum mapping
  m_schema_idx_to_known_variant_field_enum_LUT.resize_luts_if_needed(schema.attribute_num(), GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0ull;i<schema.attribute_num();++i)
  {
    auto iter = g_known_variant_field_name_to_enum.find(schema.attribute_name(i));
    if(iter != g_known_variant_field_name_to_enum.end())
      m_schema_idx_to_known_variant_field_enum_LUT.add_schema_idx_known_field_mapping(i, (*iter).second);
  }
}

void VariantQueryProcessor::initialize_v0(const VariantArraySchema& schema)
{
  m_GT_schema_version = GT_SCHEMA_V0;
}

//Added AF,AN and AC fields
void VariantQueryProcessor::initialize_v1(const VariantArraySchema& schema)
{
  //Check if any attributes in V2 schema
  const auto v1_fields = std::unordered_set<std::string>{ "AF", "AN", "AC" };
  for(auto i=0ull;i<schema.attribute_num();++i)
    if(v1_fields.find(schema.attribute_name(i)) != v1_fields.end())
    {
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V1;
      break;
    }
}

//Added GT and PS fields
void VariantQueryProcessor::initialize_v2(const VariantArraySchema& schema)
{
  //Check if any attributes in V2 schema
  const auto v2_fields = std::unordered_set<std::string>{ "GT", "PS" };
  for(auto i=0ull;i<schema.attribute_num();++i)
    if(v2_fields.find(schema.attribute_name(i)) != v2_fields.end())
    {
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V2;
      break;
    }
}

void VariantQueryProcessor::initialize_version(const VariantArraySchema& schema)
{
  initialize_known(schema);
  //Initialize to v0 schema by default
  initialize_v0(schema);
  initialize_v1(schema);
  initialize_v2(schema);
}

void VariantQueryProcessor::register_field_creators(const VariantArraySchema& schema, const VidMapper& vid_mapper)
{
  m_field_factory.resize(schema.attribute_num());
  for(auto i=0ull;i<schema.attribute_num();++i)
  {
    type_index t = schema.type(i);
    auto iter = VariantQueryProcessor::m_type_index_to_creator.find(t);
    if(iter == VariantQueryProcessor::m_type_index_to_creator.end())
      throw UnknownAttributeTypeException("Unknown type of schema attribute "+std::string(t.name()));
    //For known fields, check for special creators
    unsigned enumIdx = m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(i);
    if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(enumIdx) && KnownFieldInfo::requires_special_creator(enumIdx))
      m_field_factory.Register(i, KnownFieldInfo::get_field_creator(enumIdx));
    else
      m_field_factory.Register(i, (*iter).second);
    //TileDB does not have a way to distinguish between char, string and int8_t fields
    //Hence, a series of possibly messy checks here
    const auto& field_name = schema.attribute_name(i);
    const auto* vid_field_info = vid_mapper.get_field_info(field_name);
    if(vid_field_info)
      if(vid_field_info->m_bcf_ht_type == BCF_HT_FLAG)
      {
        auto iter = VariantQueryProcessor::m_type_index_to_creator.find(std::type_index(typeid(int8_t)));
        assert(iter != VariantQueryProcessor::m_type_index_to_creator.end());
        m_field_factory.Register(i, (*iter).second);
      }
  }
}

VariantQueryProcessor::VariantQueryProcessor(VariantStorageManager* storage_manager, const std::string& array_name,
    const VidMapper& vid_mapper)
{
  //initialize static members
  if(!VariantQueryProcessor::m_are_static_members_initialized)
    VariantQueryProcessor::initialize_static_members();
  clear();
  m_storage_manager = storage_manager;
  m_ad = storage_manager->open_array(array_name, &vid_mapper, "r");
  if(m_ad < 0)
    throw VariantQueryProcessorException("Could not open array "+array_name+" at workspace: "+storage_manager->get_workspace());
  m_array_schema = new VariantArraySchema();
  auto status = storage_manager->get_array_schema(m_ad, m_array_schema);
  assert(status == TILEDB_OK);
  m_vid_mapper = new VidMapper(vid_mapper);
  initialize();
}

VariantQueryProcessor::VariantQueryProcessor(const VariantArraySchema& array_schema, const VidMapper& vid_mapper)
{
  //initialize static members
  if(!VariantQueryProcessor::m_are_static_members_initialized)
    VariantQueryProcessor::initialize_static_members();
  clear();
  m_storage_manager = 0;
  m_array_schema = new VariantArraySchema(array_schema);
  m_vid_mapper = new VidMapper(vid_mapper);
  initialize();
}

void VariantQueryProcessor::initialize()
{
  assert(m_array_schema);
  //Initialize versioning information
  initialize_version(*m_array_schema); 
  //Register creators in factory
  register_field_creators(*m_array_schema, *m_vid_mapper);
}

void VariantQueryProcessor::obtain_TileDB_attribute_idxs(const VariantArraySchema& schema, VariantQueryConfig& queryConfig) const
{
  if(queryConfig.get_num_queried_attributes() == 0u)  //add all attributes
  {
    for(auto i=0ull;i<schema.attribute_num();++i)
      queryConfig.add_attribute_to_query(schema.attribute_name(i), i);
  }
  else
    for(auto i=0ull;i<schema.attribute_num();++i)
    {
      const auto& name = schema.attribute_name(i);
      unsigned query_idx = 0u;
      if(queryConfig.get_query_idx_for_name(name, query_idx))
        queryConfig.set_schema_idx_for_query_idx(query_idx, i);
    }
  for(auto i=0u;i<queryConfig.get_num_queried_attributes();++i)
    if(!queryConfig.is_schema_idx_defined_for_query_idx(i))
      throw UnknownQueryAttributeException("Invalid query attribute : "+queryConfig.get_query_attribute_name(i));
}

void VariantQueryProcessor::handle_gvcf_ranges(VariantCallEndPQ& end_pq,
    const VariantQueryConfig& query_config, Variant& variant,
    SingleVariantOperatorBase& variant_operator,
    int64_t& current_start_position, int64_t next_start_position, bool is_last_call, uint64_t& num_calls_with_deletions,
    GTProfileStats* stats_ptr) const
{
#ifdef DO_PROFILING
  assert(stats_ptr);
#endif
  while(!end_pq.empty() && (current_start_position < next_start_position || is_last_call) && !(variant_operator.overflow()))
  {
    int64_t top_end_pq = end_pq.top()->get_column_end();
    int64_t min_end_point = (is_last_call || (top_end_pq < (next_start_position - 1))) ? top_end_pq : (next_start_position-1);
    //If deletions, single position stepping
    min_end_point = num_calls_with_deletions ? current_start_position : min_end_point;
    //Prepare variant for aligned column interval
    variant.set_column_interval(current_start_position, min_end_point);
#ifdef DO_PROFILING
    stats_ptr->m_operator_timer.start();
    stats_ptr->update_stat(GTProfileStats::GT_NUM_OPERATOR_INVOCATIONS, 1u);
#endif
    variant_operator.operate(variant, query_config);
#ifdef DO_PROFILING
    stats_ptr->m_operator_timer.stop();
#endif
    //The following intervals have been completely processed
    while(!end_pq.empty() && static_cast<int64_t>(end_pq.top()->get_column_end()) == min_end_point)
    {
      auto top_element = end_pq.top();
      if(top_element->contains_deletion())
        --num_calls_with_deletions;
      top_element->mark_valid(false);
      end_pq.pop();
    }
    current_start_position = min_end_point + 1;   //next start position, after the end
  }
}

void VariantQueryProcessor::scan_and_operate(
    const int ad,
    const VariantQueryConfig& query_config,
    SingleVariantOperatorBase& variant_operator, unsigned column_interval_idx, bool handle_spanning_deletions,
    VariantQueryProcessorScanState* scan_state) const
{
  GTProfileStats* stats_ptr = 0;
#ifdef DO_PROFILING
  GTProfileStats stats;
  stats_ptr = &stats;
#endif
  assert(query_config.is_bookkeeping_done());
  //Priority queue of VariantCalls ordered by END positions
  VariantCallEndPQ local_end_pq;
  VariantCallEndPQ& end_pq = scan_state ? scan_state->get_end_pq() : local_end_pq;
  //Rank of tile from which scan should start
  int64_t start_column = 0;
  //Current gVCF position being operated on
  int64_t current_start_position = -1ll;
  //Variant object
  Variant local_variant;
  Variant& variant = scan_state ? scan_state->get_variant() : local_variant;
  variant.set_query_config(&query_config);
  variant.resize_based_on_query();
  //Number of calls with deletions
  uint64_t num_calls_with_deletions = scan_state ? scan_state->get_num_calls_with_deletions() : 0ull;
  //Used when deletions have to be treated as intervals and the PQ needs to be emptied
  std::vector<VariantCall*> tmp_pq_buffer(query_config.get_num_rows_to_query());
  //Forward iterator
  VariantArrayCellIterator* forward_iter = 0;
  if(scan_state && scan_state->m_iter && scan_state->m_current_start_position >= 0) //resuming a previous scan
  {
    current_start_position = scan_state->m_current_start_position;
    forward_iter = scan_state->m_iter;
#ifdef DO_PROFILING
    stats_ptr = &(scan_state->m_stats);
#endif
  }
  else //new scan
  {
    //Scan only queried interval, not whole array
    if(query_config.get_num_column_intervals() > 0u)
    {
      //If the queried interval is [100:200], then the first part of this function gets a Variant object
      //containing Calls that intersect with 100. Some of them could start before 100 and extend beyond. This
      //information is recorded in the Call
      //This part of the code accumulates such Calls, sets the current_start_position to query column interval begin
      //and lets the code in the for loop nest (forward scan) handle calling handle_gvcf_ranges()
      gt_get_column(ad, query_config, column_interval_idx, variant, stats_ptr);
      //Insert valid calls produced by gt_get_column into the priority queue
      for(Variant::valid_calls_iterator iter=variant.begin();iter != variant.end();++iter)
      {
        auto& curr_call = *iter;
        end_pq.push(&curr_call);
        if(handle_spanning_deletions && curr_call.contains_deletion())
          ++num_calls_with_deletions;
        assert(end_pq.size() <= query_config.get_num_rows_to_query());
      }
      //Valid calls were found, start position == query colum interval begin
      if(end_pq.size() > 0)
        current_start_position = query_config.get_column_begin(column_interval_idx);
      //All cells with column == query_column (or intersecting with query_column) will be handled
      //by gt_get_column(). Hence, must start from next column
      start_column = query_config.get_column_begin(column_interval_idx) + 1;
    }
    //Initialize forward scan iterators
    gt_initialize_forward_iter(ad, query_config, start_column, forward_iter);
  }
  //If uninitialized, store first column idx of forward scan in current_start_position
  if(current_start_position < 0 && !(forward_iter->end()))
  {
    auto& cell = **forward_iter;
    //Coordinates are at the start of the cell
    current_start_position = cell.get_begin_column();
  }
  //Set current column for variant (end is un-important as Calls are used to track end of intervals)
  variant.set_column_interval(current_start_position, current_start_position);
  //Next co-ordinate to consider
  int64_t next_start_position = -1ll;
  auto end_loop = false;
  for(;!(forward_iter->end()) && !end_loop && (scan_state == 0 || !(variant_operator.overflow()));++(*forward_iter))
  {
    auto& cell = **forward_iter;
#ifdef DO_PROFILING
    stats_ptr->update_stat(GTProfileStats::GT_NUM_CELLS, 1u);
    stats_ptr->update_stat(GTProfileStats::GT_NUM_ATTR_CELLS_ACCESSED, query_config.get_num_queried_attributes());
#endif
#ifdef DUPLICATE_CELL_AT_END
    //Ignore cell copies at END positions
    auto cell_column_value = cell.get_begin_column();
    auto END_v = *(cell.get_field_ptr_for_query_idx<int64_t>(query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX)));
    if(cell_column_value > END_v)
      continue;
#endif
    end_loop = scan_handle_cell(query_config, column_interval_idx, variant, variant_operator, cell,
        end_pq, tmp_pq_buffer, current_start_position, next_start_position, num_calls_with_deletions, handle_spanning_deletions, stats_ptr);
    //Do not increment the iterator if buffer overflows in the operator
    if(scan_state && variant_operator.overflow())
      break;
  }
  //Loop is over - no more data available from TileDB array
  if(end_loop || forward_iter->end())
  {
    auto is_last_call = false;
    if(query_config.get_num_column_intervals() > 0u)
    {
      next_start_position = query_config.get_column_end(column_interval_idx); //terminate at queried end
      if(next_start_position != INT64_MAX) //avoid wraparound
        ++next_start_position;
    }
    else
    {
      next_start_position = 0; //don't bother with next_start_position, forward_iter->end() must be true
      is_last_call = true;
    }
    //handle last interval
    handle_gvcf_ranges(end_pq, query_config, variant, variant_operator, current_start_position, next_start_position,
        is_last_call, num_calls_with_deletions, stats_ptr);
    if(!variant_operator.overflow())
      delete forward_iter;
#ifdef DO_PROFILING
    stats_ptr->print_stats(std::cerr);
#endif
    if(scan_state)
    {
      if(variant_operator.overflow()) //buffer full
        scan_state->set_scan_state(forward_iter, current_start_position, num_calls_with_deletions);
      else //totally done
      {
        //Invalidate iterator
        scan_state->invalidate();
        scan_state->m_done = true;
      }
    }
  }
  else  //more data available in TileDB array, but buffer is full in operator
  {
    assert(scan_state && variant_operator.overflow());
    scan_state->set_scan_state(forward_iter, current_start_position, num_calls_with_deletions);
  }
}

bool VariantQueryProcessor::scan_handle_cell(const VariantQueryConfig& query_config, unsigned column_interval_idx,
    Variant& variant, SingleVariantOperatorBase& variant_operator,
    const BufferVariantCell& cell,
    VariantCallEndPQ& end_pq, std::vector<VariantCall*>& tmp_pq_buffer,
    int64_t& current_start_position, int64_t& next_start_position,
    uint64_t& num_calls_with_deletions, bool handle_spanning_deletions,
    GTProfileStats* stats_ptr) const
{
  //If only interval requested and end of interval crossed, then done
  if(query_config.get_num_column_intervals() > 0u &&
      cell.get_begin_column() > static_cast<int64_t>(query_config.get_column_end(column_interval_idx)))
    return true;
  if(cell.get_begin_column() != current_start_position) //have found cell with next gVCF position, handle accumulated values
  {
    next_start_position = cell.get_begin_column();
    assert(cell.get_begin_column() > current_start_position);
    handle_gvcf_ranges(end_pq, query_config, variant, variant_operator, current_start_position,
        next_start_position, false, num_calls_with_deletions, stats_ptr);
    assert(end_pq.empty() || static_cast<int64_t>(end_pq.top()->get_column_end()) >= next_start_position || variant_operator.overflow());  //invariant
    //Buffer overflow, don't process anymore
    if(variant_operator.overflow())
      return false;
    //Set new start for next interval
    current_start_position = next_start_position;
    variant.set_column_interval(current_start_position, current_start_position);
    //Do not reset variant as some of the Calls that are long intervals might still be valid 
  }
  //Accumulate cells with position == current_start_position
  //Include only if row is part of query
  if(query_config.is_queried_array_row_idx(cell.get_row()))
  {
    auto& curr_call = variant.get_call(query_config.get_query_row_idx_for_array_row_idx(cell.get_row()));
    //Overlapping intervals for current call - spans across next position
    //Have to ignore rest of this interval - overwrite with the new info from the cell
    if(curr_call.is_valid() && static_cast<int64_t>(curr_call.get_column_end()) >= cell.get_begin_column())
    {
      //Have to cycle through priority queue and remove this call
      auto found_curr_call = false;
      auto num_entries_in_tmp_pq_buffer = 0ull;
      while(!end_pq.empty() && !found_curr_call)
      {
        auto top_call = end_pq.top();
        if(top_call == &curr_call)
          found_curr_call = true;
        else
        {
          assert(num_entries_in_tmp_pq_buffer < query_config.get_num_rows_to_query());
          tmp_pq_buffer[num_entries_in_tmp_pq_buffer++] = top_call;
        }
        end_pq.pop();
      }
      assert(found_curr_call);
      for(auto i=0ull;i<num_entries_in_tmp_pq_buffer;++i)
        end_pq.push(tmp_pq_buffer[i]);
      //Can handle overlapping deletions and reference blocks - if something else, throw error
      if(!curr_call.contains_deletion() && !curr_call.is_reference_block())
	throw VariantQueryProcessorException("Unhandled overlapping variants at columns "+std::to_string(curr_call.get_column_begin())+" and "
	      + std::to_string(cell.get_begin_column())+" for row "+std::to_string(cell.get_row()));
      if(curr_call.contains_deletion())
      {
	//Reduce #calls with deletions 
	assert(num_calls_with_deletions > 0u);
	--num_calls_with_deletions;
      }
    }
    curr_call.reset_for_new_interval();
    gt_fill_row(variant, cell.get_row(), cell.get_begin_column(), query_config, cell, stats_ptr);
    //When cells are duplicated at the END, then the VariantCall object need not be valid
    if(curr_call.is_valid())
    {
      end_pq.push(&curr_call);
      if(handle_spanning_deletions && curr_call.contains_deletion())
        ++num_calls_with_deletions;
      assert(end_pq.size() <= query_config.get_num_rows_to_query());
    }
  }
  return false;
}

void VariantQueryProcessor::iterate_over_cells(
    const int ad,
    const VariantQueryConfig& query_config, 
    SingleCellOperatorBase& variant_operator,
    const bool use_common_array_object) const
{
  assert(query_config.is_bookkeeping_done());
  //Initialize forward scan iterators
  SingleCellTileDBIterator* columnar_forward_iter = get_storage_manager()->begin_columnar_iterator(ad, query_config,
      use_common_array_object);
  for(;!(columnar_forward_iter->end());++(*columnar_forward_iter))
  {
    auto& cell = **columnar_forward_iter;
    auto coords = cell.get_coordinates();
    if(query_config.is_queried_array_row_idx(coords[0]))       //If row is part of query, process cell
      variant_operator.operate_on_columnar_cell(cell, query_config, get_array_schema());
  }
  variant_operator.finalize();
  delete columnar_forward_iter;
}

void VariantQueryProcessor::do_query_bookkeeping(const VariantArraySchema& array_schema,
    VariantQueryConfig& query_config, const VidMapper& vid_mapper, const bool alleles_required) const
{
  obtain_TileDB_attribute_idxs(array_schema, query_config);
  //Add END as a query attribute by default
  unsigned END_schema_idx = 
          m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_END_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(END_schema_idx));
  query_config.add_attribute_to_query("END", END_schema_idx);
  //Check if REF, ALT needs to be added as part of queried attributes
  unsigned ALT_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_ALT_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(ALT_schema_idx));
  unsigned REF_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_REF_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(REF_schema_idx));
  unsigned GT_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_GT_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(GVCF_GT_IDX));
  auto added_ALT_REF = false;
  auto added_GT = false;
  //Required by the caller
  if(alleles_required)
  {
    query_config.add_attribute_to_query("ALT", ALT_schema_idx);
    query_config.add_attribute_to_query("REF", REF_schema_idx);
    added_ALT_REF = true;
  }
  //As attributes are added to the query, the #queried attributes increases
  //So, do not store this number into a scalar before the loop
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    const auto& field_name = query_config.get_query_attribute_name(i);
    const auto* vid_field_info = vid_mapper.get_field_info(field_name);
    auto length_descriptor = FieldLengthDescriptor();
    if(vid_field_info)
    {
      length_descriptor = vid_field_info->m_length_descriptor;
      query_config.set_query_attribute_info_parameters(i,
          array_schema.type(schema_idx), length_descriptor,
          vid_field_info->m_VCF_field_combine_operation);
    }
    else //No information in vid file, see if something can be gleaned from known fields
    {
      auto known_field_enum = 0u;
      if(KnownFieldInfo::get_known_field_enum_for_name(field_name, known_field_enum))
      {
        length_descriptor.set_length_descriptor(0,
            KnownFieldInfo::get_length_descriptor_for_known_field_enum(known_field_enum));
        query_config.set_query_attribute_info_parameters(i,
            array_schema.type(schema_idx), length_descriptor,
            KnownFieldInfo::get_VCF_field_combine_operation_for_known_field_enum(known_field_enum)
            );
      }
    }
    //Does the length of the field depend on the number of alleles? If yes, add ALT and REF as query fields
    if(!added_ALT_REF && length_descriptor.is_length_allele_dependent())
    {
      query_config.add_attribute_to_query("ALT", ALT_schema_idx);
      query_config.add_attribute_to_query("REF", REF_schema_idx);
      added_ALT_REF = true;
    }
    //Does the length of the field depend on the ploidy? If yes, add GT field
    if(!added_GT && length_descriptor.is_length_genotype_dependent())
    {
      query_config.add_attribute_to_query("GT", GT_schema_idx);
      added_GT = true;
    }
  }
  //Re-order query fields so that special fields are first
  query_config.reorder_query_fields();
  //Set known field enum within query
  query_config.resize_LUT(GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    unsigned schema_idx = query_config.get_schema_idx_for_query_idx(i);
    unsigned known_variant_field_enum = 
      m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(schema_idx);
    if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(known_variant_field_enum))
    {
      query_config.add_query_idx_known_field_enum_mapping(i, known_variant_field_enum);
      assert(g_known_variant_field_names[known_variant_field_enum] == query_config.get_query_attribute_name(i));
    }
  }
  //Either not added or both REF and ALT must be part of query
  assert(!added_ALT_REF
      || (query_config.is_defined_query_idx_for_known_field_enum(GVCF_REF_IDX)
        && query_config.is_defined_query_idx_for_known_field_enum(GVCF_ALT_IDX)));
  assert(!added_GT
      || query_config.is_defined_query_idx_for_known_field_enum(GVCF_GT_IDX));
  //Set number of rows in the array
  auto& dim_domains = array_schema.dim_domains();
  uint64_t row_num = m_storage_manager ? m_storage_manager->get_num_valid_rows_in_array(m_ad) :   //may read from array metadata
    (dim_domains[0].second - dim_domains[0].first + 1);
  query_config.set_num_rows_in_array(row_num, static_cast<int64_t>(dim_domains[0].first));
  query_config.setup_array_row_idx_to_query_row_idx_map();
  //Bounds checking for query
  for(auto i=0u;i<query_config.get_num_column_intervals();++i)
  {
    const auto& col_range = query_config.get_column_interval(i);
    if(col_range.first < dim_domains[1].first || col_range.first > dim_domains[1].second
        || col_range.second < dim_domains[1].first || col_range.second > dim_domains[1].second)
      throw OutOfBoundsQueryException("Query interval "+std::to_string(i)+" : "+std::to_string(col_range.first)+", "
          +std::to_string(col_range.second)+" is out of bounds");
  }
  //If specific rows requested
  if(!query_config.query_all_rows())
    for(uint64_t i=0ull;i<query_config.get_num_rows_to_query();++i)
    {
      auto row_idx = query_config.get_rows_to_query()[i];
      if(row_idx < dim_domains[0].first || row_idx > dim_domains[0].second)
        throw OutOfBoundsQueryException("Queried row index "+std::to_string(row_idx)+" is out of bounds");
    }
  //Done with bookkeeping
  query_config.set_done_bookkeeping(true);
}

void VariantQueryProcessor::gt_get_column_interval(
    const int ad,
    const VariantQueryConfig& query_config, unsigned column_interval_idx,
    vector<Variant>& variants, GA4GHPagingInfo* paging_info, GTProfileStats* stats_ptr) const {
#ifdef DO_PROFILING
  assert(stats_ptr);
#endif
  if(paging_info)
    paging_info->init_page_query();
  uint64_t start_variant_idx = variants.size();
  //Will be used later in the function to produce Variants with one CallSet
  VariantQueryConfig subset_query_config(query_config);
  vector<int64_t> subset_rows = vector<int64_t>(1u, query_config.get_smallest_row_idx_in_array());
  subset_query_config.update_rows_to_query(subset_rows);  //only 1 row, row 0
  //Structure that helps merge multiple Calls into a single variant if the GA4GH specific merging
  //conditions are satisfied
  GA4GHCallInfoToVariantIdx call_info_2_variant;
  //If the queried interval is [100:200], then the first part of the function gets a Variant object
  //containing Calls that intersect with 100. Some of them could start before 100 and extend beyond. This
  //information is recorded in the Call
  //With paging, if this is a continuation of the query and the continuation should occur <= the same column 
  //as the start of the query, the left sweep operation should be repeated. However, if the continuation is beyond
  //the start of the query, skip the left sweep operation (as all variants accessed by the left sweep would have
  //been returned in a previous page)
  if(paging_info == 0 || paging_info->get_last_column() <= query_config.get_column_begin(column_interval_idx))
  {
    Variant interval_begin_variant(&query_config);
    interval_begin_variant.resize_based_on_query();
    //If cells are duplicated, no claim can be made  on the order in which cells are traversed since
    //the order of END cells has no bearing on the order of begin values
    //If not duplicated, row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
    //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
    vector<uint64_t> query_row_idx_in_order = vector<uint64_t>(query_config.get_num_rows_to_query(), UNDEFINED_NUM_ROWS_VALUE);
#if VERBOSE>0
    std::cerr << "[query_variants:gt_get_column_interval] Getting " << query_config.get_num_rows_to_query() << " rows" << std::endl;
#endif
    gt_get_column(ad, query_config, column_interval_idx, interval_begin_variant, stats_ptr,
#ifdef DUPLICATE_CELL_AT_END
        0
#else
        &query_row_idx_in_order
#endif
        );
    //This interval contains many Calls, likely un-aligned (no common start/end). Split this variant 
    //into multiple Variants, each containing calls that are satisfy the GA4GH properties for merging calls
    interval_begin_variant.move_calls_to_separate_variants(query_config, variants, query_row_idx_in_order,
        call_info_2_variant, paging_info);
  }
  //If this is not a single position query and paging limit is not hit, need to fetch more cells
  if((query_config.get_column_end(column_interval_idx) >
        query_config.get_column_begin(column_interval_idx)) && 
      !(paging_info && paging_info->is_page_limit_hit(variants.size())))
  {
    //Get the iterator to the first cell that has column > query_column_start. All cells with column == query_column_start
    //(or intersecting) would have been handled by gt_get_column(). Hence, must start from next column
    uint64_t start_column_forward_sweep = query_config.get_column_interval(column_interval_idx).first+1u;
    //If paging, continue at the last column that was handled in the previous page
    start_column_forward_sweep = paging_info ? std::max<uint64_t>(paging_info->get_last_column(), start_column_forward_sweep) 
      : start_column_forward_sweep;
    VariantArrayCellIterator* forward_iter = 0;
    gt_initialize_forward_iter(ad, query_config, query_config.get_column_interval(column_interval_idx).first+1, forward_iter);
    //Used to store single call variants  - one variant per cell
    //Multiple variants could be merged later on
    Variant tmp_variant(&subset_query_config);
    tmp_variant.resize_based_on_query();
#if VERBOSE>0
    std::cerr << "[query_variants:gt_get_column_interval] Fetching columns from " << query_config.get_column_begin(column_interval_idx) + 1;
    std::cerr << " to " << query_config.get_column_end(column_interval_idx) << std::endl;
#endif
    //Used for paging
    auto last_column_idx = start_column_forward_sweep;
    //Num handled variants - if beginning from same column as end of last page, get from paging info
    //else, reset to 0
    auto num_last_column_variants_handled_after_curr_page = paging_info ? 
      (paging_info->get_last_column() == last_column_idx) ? paging_info->get_num_handled_variants_in_last_column() : 0u
      : 0u;
    num_last_column_variants_handled_after_curr_page = 0u;
    auto curr_column_idx = start_column_forward_sweep;
    auto curr_row_idx = 0ull;
    bool stop_inserting_new_variants = false;
    for(;!(forward_iter->end());++(*forward_iter))
    {
#ifdef DO_PROFILING
      stats_ptr->update_stat(GTProfileStats::GT_NUM_CELLS, 1u);
      //FIXME: in the current implementation, every *iter accesses every attribute
      stats_ptr->update_stat(GTProfileStats::GT_NUM_ATTR_CELLS_ACCESSED, query_config.get_num_queried_attributes());
#endif
      auto& cell = **forward_iter;
      curr_row_idx = cell.get_row();
      curr_column_idx = cell.get_begin_column();
      if(curr_column_idx > query_config.get_column_end(column_interval_idx))    //Genomic interval begins after end of query interval
        break;
      //Check variants handled in previous page
      //TODO: should never enter this if statement, I think
      if(paging_info && paging_info->handled_previously(curr_row_idx, curr_column_idx))
        continue;
      if(query_config.is_queried_array_row_idx(curr_row_idx))       //If row is part of query, process cell
      {
        //Create Variant with single Call (subset_query_config contains one row)
        subset_rows[0] = curr_row_idx;
        subset_query_config.update_rows_to_query(subset_rows);
        tmp_variant.resize_based_on_query();
        assert(tmp_variant.get_num_calls() == 1u);      //exactly 1 call
        tmp_variant.reset_for_new_interval();
        tmp_variant.get_call(0u).set_row_idx(curr_row_idx); //set row idx
        tmp_variant.set_column_interval(curr_column_idx, curr_column_idx);
        gt_fill_row(tmp_variant, curr_row_idx, curr_column_idx, subset_query_config, cell, stats_ptr);
        assert(tmp_variant.get_num_calls() == 1u);      //exactly 1 call
        //When cells are duplicated at the END, then the VariantCall object need not be valid
        if(tmp_variant.get_call(0).is_valid())
        {
          //Move call to variants vector, creating new Variant if necessary
          auto newly_inserted = move_call_to_variant_vector(subset_query_config, tmp_variant.get_call(0), variants, call_info_2_variant,
              stop_inserting_new_variants);
          //Check if page limit hit
          PAGE_END_CHECK_LOGIC
        }
      }
      last_column_idx = curr_column_idx;
    }
#if VERBOSE>0
    std::cerr << "[query_variants:gt_get_column_interval] Fetching columns complete " << std::endl;
#endif
    delete forward_iter;
  }
  if(paging_info)
  {
    //Exited loop without hitting page limit - only reason, end of query
    if(!(paging_info->is_page_limit_hit(variants.size())))
      paging_info->set_query_completed(variants);
    //Remove variants handled in the previous page(s)
    paging_info->shift_left_variants(variants);
  }
#if VERBOSE>0
  std::cerr << "[query_variants:gt_get_column_interval] re-arrangement of variants " << std::endl;
#endif
  GA4GHOperator variant_operator(query_config);
  for(auto i=start_variant_idx;i<variants.size();++i)
    if(variants[i].get_num_calls() > 1u) //possible re-arrangement of PL/AD/GT fields needed
    {
#ifdef DO_PROFILING
      stats_ptr->m_operator_timer.start();
#endif
      variant_operator.operate(variants[i], query_config);
      variant_operator.copy_back_remapped_fields(variants[i]); //copy back fields that have been remapped
#ifdef DO_PROFILING
      stats_ptr->m_operator_timer.stop();
#endif
    }
  if(paging_info)
    paging_info->serialize_page_end(m_array_schema->array_name());
#if VERBOSE>0
  std::cerr << "[query_variants:gt_get_column_interval] query complete " << std::endl;
#endif
}

void VariantQueryProcessor::gt_get_column(
    const int ad,
    const VariantQueryConfig& query_config, unsigned column_interval_idx,
    Variant& variant, GTProfileStats* stats_ptr, std::vector<uint64_t>* query_row_idx_in_order) const {
#ifdef DO_PROFILING
  assert(stats_ptr);
  stats_ptr->m_interval_sweep_timer.start();
#endif
  //New interval starts
  variant.reset_for_new_interval();

  assert(query_config.get_num_column_intervals() > 0u && column_interval_idx < query_config.get_num_column_intervals());
  assert(query_config.is_bookkeeping_done() && "Before calling gt_get_column(), do_query_bookkeeping() function must be called");
  assert(query_config.get_first_normal_field_query_idx() >= 1u); //END is required by default and must be the first attributes (idx 0)
  
  variant.set_column_interval(query_config.get_column_interval(column_interval_idx).first,
          query_config.get_column_interval(column_interval_idx).second);
  //TODO: Still single position query
  uint64_t col = query_config.get_column_interval(column_interval_idx).first;
#if VERBOSE>0
  std::cerr << "[query_variants:gt_get_column] Fetching column : " << col << std::endl;
#endif
#ifdef DUPLICATE_CELL_AT_END
  //If cells are duplicated at the end, we only need a forward iterator starting at col
  //i.e. start at the smallest cell with co-ordinate >= col
  VariantArrayCellIterator* cell_iter = 0;
  gt_initialize_forward_iter(ad, query_config, query_config.get_column_interval(column_interval_idx).first, cell_iter);
#endif //ifdef DUPLICATE_CELL_AT_END
  // Indicates how many rows have been filled.
  uint64_t filled_rows = 0;
  uint64_t num_valid_rows = 0;
  // Fill the genotyping column
  while(!(cell_iter->end()) && filled_rows < query_config.get_num_rows_to_query()) {
#ifdef DO_PROFILING
    stats_ptr->update_stat(GTProfileStats::GT_NUM_CELLS, 1u);
    stats_ptr->update_stat(GTProfileStats::GT_NUM_CELLS_IN_LEFT_SWEEP, 1u);
    //FIXME: in the current implementation, every *iter accesses every attribute
    stats_ptr->update_stat(GTProfileStats::GT_NUM_ATTR_CELLS_ACCESSED, query_config.get_num_queried_attributes());
#endif
    auto& cell = **cell_iter;
#ifdef DUPLICATE_CELL_AT_END
    // If next cell is not on the left of col, and
    // The rowIdx is being queried and
    // The row/call is uninitialized (uninvestigated) in the Variant
    if(cell.get_begin_column() >= static_cast<int64_t>(col) && query_config.is_queried_array_row_idx(cell.get_row()))
#else
    // If next cell is not on the right of col, and
    // The rowIdx is being queried and
    // The row/call is uninitialized (uninvestigated) in the Variant
    if(cell.get_begin_column() <= static_cast<int64_t>(col) && query_config.is_queried_array_row_idx(cell.get_row()))
#endif
    {
      auto curr_query_row_idx = query_config.get_query_row_idx_for_array_row_idx(cell.get_row());
      auto& curr_call = variant.get_call(curr_query_row_idx);
      if(!(curr_call.is_initialized()))
      {
        gt_fill_row(variant, cell.get_row(), cell.get_begin_column(), query_config, cell, stats_ptr
#ifdef DUPLICATE_CELL_AT_END
            , true
#endif
            );
        ++filled_rows;
#ifndef DUPLICATE_CELL_AT_END
        //If cells are NOT duplicated, then the order in which cells are traversed is reverse of column-major order
        //If cells are duplicated, no claim can be made since the order of END cells has no bearing on the order of
        //begin values
        if(curr_call.is_valid() && query_row_idx_in_order)
        {
          assert(num_valid_rows < query_row_idx_in_order->size());
          (*query_row_idx_in_order)[num_valid_rows++] = curr_query_row_idx;
        }
#endif
      }
    }
    ++(*cell_iter);
  }
  //Free memory 
  //if(cell.cell())
  //free(const_cast<void*>(cell.cell()));
  delete cell_iter;

#ifndef DUPLICATE_CELL_AT_END
  if(query_row_idx_in_order)
    query_row_idx_in_order->resize(num_valid_rows);
#endif

#ifdef DO_PROFILING
  stats_ptr->m_interval_sweep_timer.stop();
#endif
}

void VariantQueryProcessor::fill_field_prep(std::unique_ptr<VariantFieldBase>& field_ptr,
    const VariantQueryConfig& query_config, const unsigned query_idx) const
{
  auto schema_idx = query_config.get_schema_idx_for_query_idx(query_idx);
  if(field_ptr.get() == nullptr)       //Allocate only if null
    field_ptr = std::move(m_field_factory.Create(schema_idx, m_array_schema->is_variable_length_field(schema_idx)));
  field_ptr->set_valid(true);  //mark as valid
}

void VariantQueryProcessor::fill_field(std::unique_ptr<VariantFieldBase>& field_ptr,
    const BufferVariantCell::FieldsIter& attr_iter,
    const VariantQueryConfig& query_config, const unsigned query_idx
    ) const
{
  fill_field_prep(field_ptr, query_config, query_idx);
  //This function might mark the field as invalid - some fields are  determined to be invalid only
  //after accessing the data and comparing to NULL_* values
  field_ptr->copy_data_from_tile(attr_iter);
}

void VariantQueryProcessor::binary_deserialize(Variant& variant, const VariantQueryConfig& query_config,
    const vector<uint8_t>& buffer, uint64_t& offset) const
{
  assert(offset < buffer.size());
  //deserialize header
  variant.binary_deserialize_header(buffer, offset, query_config.get_num_queried_attributes());
  //VariantCall info
  for(auto i=0ull;i<variant.get_num_calls();++i)
  {
    auto& curr_call = variant.get_call(i);
    curr_call.binary_deserialize_header(buffer, offset);
    //Fields
    assert(query_config.get_num_queried_attributes() == curr_call.get_num_fields());
    for(auto j=0u;j<curr_call.get_num_fields();++j)
    {
      //check if field is valid
      auto is_valid_field = *(reinterpret_cast<const bool*>(&(buffer[offset])));
      offset += sizeof(bool);
      if(is_valid_field)
      {
        auto& field_ptr = curr_call.get_field(j); 
        fill_field_prep(field_ptr, query_config, j);
        auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(j);
        field_ptr->binary_deserialize(reinterpret_cast<const char*>(&(buffer[0])), offset,
            !length_descriptor.is_fixed_length_field(),
            length_descriptor.is_fixed_length_field() ? length_descriptor.get_num_elements() : 0u);
      }
    }
  }
  //Common fields in the Variant object
  for(auto i=0u;i<variant.get_num_common_fields();++i)
  {
    //Flag representing whether common field is not null
    auto is_valid_field = *(reinterpret_cast<const bool*>(&(buffer[offset])));
    offset += sizeof(bool);
    //Query idx for common field
    auto query_idx = *(reinterpret_cast<const unsigned*>(&(buffer[offset])));
    offset += sizeof(unsigned);
    variant.set_query_idx_for_common_field(i, query_idx);
    if(is_valid_field)
    {
      std::unique_ptr<VariantFieldBase>& field_ptr = variant.get_common_field(i); 
      fill_field_prep(field_ptr, query_config, query_idx);
      auto length_descriptor = query_config.get_length_descriptor_for_query_attribute_idx(query_idx);
      field_ptr->binary_deserialize(reinterpret_cast<const char*>(&(buffer[0])), offset,
          !length_descriptor.is_fixed_length_field(),
          length_descriptor.is_fixed_length_field() ? length_descriptor.get_num_elements() : 0u);
    }
  }
}

void VariantQueryProcessor::gt_fill_row(
    Variant& variant, int64_t row, int64_t column,
    const VariantQueryConfig& query_config,
    const BufferVariantCell& cell, GTProfileStats* stats_ptr
#ifdef DUPLICATE_CELL_AT_END
    , bool traverse_end_copies
#endif
    ) const {
#ifdef DO_PROFILING
  assert(stats_ptr);
  stats_ptr->m_genomicsdb_cell_fill_timer.start();
#endif
#if VERBOSE>1
  std::cerr << "[query_variants:gt_fill_row] Fill Row " << row << " column " << column << std::endl;
#endif
  //Current row should be part of query
  assert(query_config.is_queried_array_row_idx(row));
  VariantCall& curr_call = variant.get_call(query_config.get_query_row_idx_for_array_row_idx(row));
  //Curr call will be initialized, one way or the other
  curr_call.mark_initialized(true);
  curr_call.set_contains_deletion(false);
  curr_call.set_is_reference_block(false);
  //Column values
  auto query_column_value = static_cast<int64_t>(variant.get_column_begin());
  auto cell_begin_value = column;
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto* END_ptr = cell.get_field_ptr_for_query_idx<int64_t>(query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto END_v = *END_ptr;
  // First check if the row contains valid data, i.e., check whether the interval intersects with the current queried interval
#ifdef DUPLICATE_CELL_AT_END
  //Counterpart of reverse iterator when traverse_end_copies == true
  //If this is a begin cell and it begins after the query interval, mark invalid
  //If this is an END cell, then its begin value is guaranteed to be before the query_column_value. This
  //is because if the begin is AFTER query_column_value, then the begin cell would have been traversed first
  //by gt_get_column and the curr_call would already be marked initialized and invalid.
  //
  //When traverse_end_copies == false, we are doing a forward traversal and any END copies of cells
  //should be returned as invalid (and must be handled by the caller)
  if((traverse_end_copies && cell_begin_value <= END_v && cell_begin_value > query_column_value)
      || (!traverse_end_copies && cell_begin_value > END_v))
#else
  //When no duplicate cells are present and the reverse iterator is used, if the cell ends before the query column,
  //mark the cell as initialized but invalid as the row has no valid data for this query position
  if(END_v < query_column_value)
#endif
  {
    curr_call.mark_valid(false);
    return;
  }
  curr_call.mark_valid(true);   //contains valid data for this query
#ifdef DO_PROFILING
  stats_ptr->update_stat(GTProfileStats::GT_NUM_VALID_CELLS_IN_QUERY, 1u);
#endif
#ifdef DUPLICATE_CELL_AT_END
  if(column > END_v)
    std::swap(column, END_v);
#endif
  //Set begin,end of the Call - NOTE: need not be same as Variant's begin,end
  curr_call.set_column_interval(column, END_v);
  //END should be the first queried attribute - see reorder_query_fields()
  assert(query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX) < 1u);
  assert(query_config.get_first_normal_field_query_idx() >= 1u);
  assert(query_config.get_first_normal_field_query_idx() != UNDEFINED_ATTRIBUTE_IDX_VALUE);
  //Variables to store special fields
  //Num alternate alleles
  unsigned num_ALT_alleles = 0u;
  //ploidy
  unsigned ploidy = 0u;
  //Iterate over attributes
  auto attr_iter = cell.begin();
  ++attr_iter;  //skip the END field
  //First, load special fields up to and including ALT
  for(auto i=1u;i<query_config.get_first_normal_field_query_idx();++i,++attr_iter)
  {
    //Read from Tile
    fill_field(curr_call.get_field(i), attr_iter,
        query_config, i
        );
  }
  //Initialize ALT field, if needed
  const auto* ALT_field_ptr = get_known_field_if_queried<VariantFieldALTData, true>(curr_call, query_config, GVCF_ALT_IDX); 
  if(ALT_field_ptr && ALT_field_ptr->is_valid())
    num_ALT_alleles = ALT_field_ptr->get().size();   //ALT field data is vector<string>
  //Go over all normal query fields and fetch data
  for(auto i=query_config.get_first_normal_field_query_idx();i<query_config.get_num_queried_attributes();++i, ++attr_iter)
  {
    //Read from Tile
    fill_field(curr_call.get_field(i), attr_iter,
        query_config, i
        );     
  }
  //Initialize REF field, if queried
  const auto* REF_field_ptr = get_known_field_if_queried<VariantFieldString, true>(curr_call, query_config, GVCF_REF_IDX); 
  //Check for deletion
  if(REF_field_ptr && REF_field_ptr->is_valid() && ALT_field_ptr && ALT_field_ptr->is_valid())
  {
    auto has_deletion = VariantUtils::contains_deletion(REF_field_ptr->get(), ALT_field_ptr->get());
    curr_call.set_contains_deletion(has_deletion);
    curr_call.set_is_reference_block(VariantUtils::is_reference_block(REF_field_ptr->get(), ALT_field_ptr->get()));
  }
#ifdef DO_PROFILING
  stats_ptr->m_genomicsdb_cell_fill_timer.stop();
#endif
}

inline
unsigned int VariantQueryProcessor::gt_initialize_forward_iter(
    const int ad,
    const VariantQueryConfig& query_config, const int64_t column,
    VariantArrayCellIterator*& forward_iter) const {
  assert(query_config.is_bookkeeping_done());
  //Num attributes in query
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  //Assign forward iterator
  vector<int64_t> query_range = { query_config.get_smallest_row_idx_in_array(),
    static_cast<int64_t>(query_config.get_num_rows_in_array()+query_config.get_smallest_row_idx_in_array()-1),
    column, INT64_MAX };
  forward_iter = get_storage_manager()->begin(ad, &(query_range[0]), query_config.get_query_attributes_schema_idxs());
  return num_queried_attributes - 1;
}

void VariantQueryProcessor::clear()
{
  m_schema_idx_to_known_variant_field_enum_LUT.reset_luts();
  m_field_factory.clear();
}

