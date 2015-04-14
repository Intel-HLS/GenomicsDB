#include "gt_common.h"
#include "query_variants.h"
#include "variant_operations.h"

using namespace std;

//Static members
bool VariantQueryProcessor::m_are_static_members_initialized = false;
vector<string> VariantQueryProcessor::m_known_variant_field_names = vector<string>{
    "END",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "BASEQRANKSUM",
    "CLIPPINGRANKSUM",
    "MQRANKSUM",
    "READPOSRANKSUM",
    "DP",
    "MQ",
    "MQ0",
    "DP_FMT",
    "MIN_DP",
    "GQ",
    "SB_1",
    "SB_2",
    "SB_3",
    "SB_4",
    "AD",
    "PL",
    "AF",
    "AN",
    "AC",
    "NULL",
    "OFFSETS",
    "COORDINATES"
};
unordered_map<string, unsigned> VariantQueryProcessor::m_known_variant_field_name_to_enum;
//Initialize static members function
void VariantQueryProcessor::initialize_static_members()
{
  for(auto i=0u;i<VariantQueryProcessor::m_known_variant_field_names.size();++i)
    VariantQueryProcessor::m_known_variant_field_name_to_enum[VariantQueryProcessor::m_known_variant_field_names[i]] = i;
  VariantQueryProcessor::m_are_static_members_initialized = true;
}

void VariantQueryProcessor::initialize_known(const StorageManager::ArrayDescriptor* ad)
{
  //Initialize NULL bit idx for all known variant fields
  m_NULL_bit_enum_idx_vec[GVCF_AC_IDX] =                        GVCF_AC_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_AN_IDX] =                        GVCF_AN_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_AF_IDX] =                        GVCF_AF_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_PL_IDX] =                        GVCF_PL_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_AD_IDX] =                        GVCF_AD_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_SB_4_IDX] =                      GVCF_SB_4_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_SB_3_IDX] =                      GVCF_SB_3_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_SB_2_IDX] =                      GVCF_SB_2_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_SB_1_IDX] =                      GVCF_SB_1_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_GQ_IDX] =                        GVCF_GQ_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_MIN_DP_IDX] =                    GVCF_MIN_DP_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_DP_FMT_IDX] =                    GVCF_DP_FMT_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_MQ0_IDX] =                       GVCF_MQ0_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_MQ_IDX] =                        GVCF_MQ_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_DP_IDX] =                        GVCF_DP_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_READPOSRANKSUM_IDX] =            GVCF_READPOSRANKSUM_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_MQRANKSUM_IDX] =                 GVCF_MQRANKSUM_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_CLIPPINGRANKSUM_IDX] =           GVCF_CLIPPINGRANKSUM_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_BASEQRANKSUM_IDX] =              GVCF_BASEQRANKSUM_NULL_BITIDX;
  m_NULL_bit_enum_idx_vec[GVCF_QUAL_IDX] =                      GVCF_QUAL_NULL_BITIDX;
  //Initialize schema idx <--> known_field_enum mapping
  //+1 for the coordinates
  const auto& schema = ad->array_schema();
  m_schema_idx_to_known_variant_field_enum_LUT.resize_luts_if_needed(schema.attribute_num()+1, GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0u;i<schema.attribute_num();++i)
  {
    auto iter = VariantQueryProcessor::m_known_variant_field_name_to_enum.find(schema.attribute_name(i));
    if(iter != VariantQueryProcessor::m_known_variant_field_name_to_enum.end())
      m_schema_idx_to_known_variant_field_enum_LUT.add_schema_idx_known_field_mapping(i, (*iter).second);
  }
  m_schema_idx_to_known_variant_field_enum_LUT.add_schema_idx_known_field_mapping(schema.attribute_num(), GVCF_COORDINATES_IDX); 
  //Initialize which fields use OFFSETS
  m_known_field_enum_uses_offset[GVCF_REF_IDX] = true;
  m_known_field_enum_uses_offset[GVCF_ALT_IDX] = true;
  m_known_field_enum_uses_offset[GVCF_FILTER_IDX] = true;
  m_known_field_enum_uses_offset[GVCF_PL_IDX] = true;
  m_known_field_enum_uses_offset[GVCF_AD_IDX] = true;
}

void VariantQueryProcessor::initialize_v0(const StorageManager::ArrayDescriptor* ad)
{
  //invalidate null bit idx for AF, AC, AN
  for(auto i=GVCF_PL_IDX+1;i<GVCF_NUM_KNOWN_FIELDS;++i)
    invalidate_NULL_bitidx(i);
  //subtract GVCF_PL_NULL_BITIDX from other defined NULL bitmap values 
  for(auto i=0u;i<=GVCF_PL_IDX;++i)
  {
    if(is_NULL_bitidx_defined_for_known_field_enum(i))
    {
      assert(m_NULL_bit_enum_idx_vec[i] >= GVCF_PL_NULL_BITIDX);
      m_NULL_bit_enum_idx_vec[i] -= GVCF_PL_NULL_BITIDX;
    }
  }
  m_GT_schema_version = GT_SCHEMA_V0;
}

void VariantQueryProcessor::initialize_v1(const StorageManager::ArrayDescriptor* ad)
{
  const auto& schema = ad->array_schema();
  //Check if any attributes in V2 schema
  const auto v2_fields = std::unordered_set<std::string>{ "AF", "AN", "AC" };
  for(auto i=0u;i<schema.attribute_num();++i)
    if(v2_fields.find(schema.attribute_name(i)) != v2_fields.end())
    {
      unsigned diff = (GVCF_PL_NULL_BITIDX - GVCF_AC_NULL_BITIDX);    //v2 only knows upto field AC
      //Reverse of v1, increment NULL bitix of defined fields
      for(auto i=0u;i<=GVCF_PL_IDX;++i)
      {
        if((is_NULL_bitidx_defined_for_known_field_enum(i)))
          m_NULL_bit_enum_idx_vec[i] +=  diff;
      }
      diff = (GVCF_AF_NULL_BITIDX - GVCF_AC_NULL_BITIDX);
      //validate null bit idx for AF, AC, AN
      for(unsigned i=GVCF_AF_IDX;i<=GVCF_AC_IDX;++i)
        m_NULL_bit_enum_idx_vec[i] =  diff - (i - GVCF_AF_IDX);
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V1;
      break;
    }
}

void VariantQueryProcessor::initialize_version(const StorageManager::ArrayDescriptor* ad)
{
  initialize_known(ad);
  //Initialize to v0 schema by default
  initialize_v0(ad);
  initialize_v1(ad);
}

VariantQueryProcessor::VariantQueryProcessor(const std::string& workspace, StorageManager& storage_manager,
    const StorageManager::ArrayDescriptor* ad)
: QueryProcessor(workspace, storage_manager)
{
  if(!VariantQueryProcessor::m_are_static_members_initialized)
    VariantQueryProcessor::initialize_static_members();
  clear();
  //Invalidate m_NULL_bit_enum_idx_vec 
  m_NULL_bit_enum_idx_vec.resize(GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0u;i<m_NULL_bit_enum_idx_vec.size();++i)
    invalidate_NULL_bitidx(i);
  //set all fields NOT to use offset
  m_known_field_enum_uses_offset.resize(GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0u;i<m_known_field_enum_uses_offset.size();++i)
    m_known_field_enum_uses_offset[i] = false;
  //Initialize versioning information
  initialize_version(ad);
}

void VariantQueryProcessor::handle_gvcf_ranges(VariantIntervalPQ& end_pq, std::vector<PQStruct>& PQ_end_vec,
    const VariantQueryConfig& query_config, GTColumn* gt_column,
      std::unordered_map<uint64_t, GTTileIteratorsTracker>& tile_idx_2_iters, std::ostream& output_stream,
    int64_t current_start_position, int64_t next_start_position, bool is_last_call) const
{
  uint64_t num_deref_tile_iters = 0;
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
        gt_fill_row<StorageManager::const_iterator>(gt_column, i, curr_struct.m_array_column, curr_struct.m_cell_pos, query_config, 
	    &((*find_iter).second.m_iter_vector[0]), &num_deref_tile_iters);
      }
    }
    VariantOperations::do_dummy_genotyping(gt_column, output_stream);
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

void VariantQueryProcessor::iterate_over_all_tiles(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config) const
{
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
  
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  for(auto i=0u;i<num_queried_attributes;++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    tile_its[i] = get_storage_manager().begin(ad, schema_idx);
  }
  //Get tile end iter for COORDS
  StorageManager::const_iterator tile_it_end = get_storage_manager().end(ad,
      m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  //Get query idx for COORDS
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  
  GTColumn* gt_column = new GTColumn(0, total_num_samples);
  uint64_t num_deref_tile_iters = 0; 
  for(;tile_its[COORDS_query_idx] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))
  {
    Tile::const_iterator cell_it = (*tile_its[COORDS_query_idx]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[COORDS_query_idx]).end();
    std::vector<int64_t> next_coord = *cell_it;
    gt_column->reset();
    gt_column->col_ = next_coord[1];
    gt_fill_row<StorageManager::const_iterator>(gt_column, next_coord[0], next_coord[1], cell_it.pos(), query_config, tile_its,
	&num_deref_tile_iters);
  }
}

void VariantQueryProcessor::scan_and_operate(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config,
    std::ostream& output_stream) const
{
  assert(query_config.is_bookkeeping_done());
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
  
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  std::unordered_map<uint64_t, GTTileIteratorsTracker> tile_idx_2_iters;
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    tile_its[i] = get_storage_manager().begin(ad, schema_idx);
  }
  //Get tile end iter for COORDS
  StorageManager::const_iterator tile_it_end = get_storage_manager().end(ad,
      m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  //Get query idx for COORDS and END
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  unsigned END_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX);

  GTColumn* gt_column = new GTColumn(0, total_num_samples);

  //Priority queue of END positions
  VariantIntervalPQ end_pq;
  //Vector of PQStruct - pre-allocate to eliminate allocations inside the while loop
  //Elements of the priority queue end_pq are pointers to elements of this vector
  auto PQ_end_vec = std::vector<PQStruct>(total_num_samples, PQStruct{});
  for(auto i=0ull;i<PQ_end_vec.size();++i)
    PQ_end_vec[i].m_sample_idx = i;
  //Current gVCF position being operated on
  int64_t current_start_position = -1ll;
  //Get first valid position in the array
  if(tile_its[COORDS_query_idx] != tile_it_end)
  {
    Tile::const_iterator cell_it = (*(tile_its[COORDS_query_idx])).begin();
    Tile::const_iterator cell_it_end = (*(tile_its[COORDS_query_idx])).end();
    if(cell_it != cell_it_end)
    {
      std::vector<int64_t> current_coord = *cell_it;
      current_start_position = current_coord[1];
    }
  }
  int64_t next_start_position = -1ll;
  uint64_t tile_idx = 0ull;
  for(;tile_its[COORDS_query_idx] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))  //why -1, advance uses i<=N for loop
  {
    //Setup map for tile id to iterator
    //FIXME: possible use of boost::pool for reducing overhead of small allocs
    auto insert_iter = tile_idx_2_iters.emplace(std::pair<uint64_t, GTTileIteratorsTracker>(tile_idx, GTTileIteratorsTracker(num_queried_attributes)));
    GTTileIteratorsTracker& current_iterator_tracker = insert_iter.first->second;
    for(auto i=0u;i<num_queried_attributes;++i)
      current_iterator_tracker.m_iter_vector[i] = tile_its[i];
    // Initialize cell iterators for the coordinates
    Tile::const_iterator cell_it = (*tile_its[COORDS_query_idx]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[COORDS_query_idx]).end();
    for(;cell_it != cell_it_end;++cell_it) {
      std::vector<int64_t> next_coord = *cell_it;
      if(next_coord[1] != current_start_position) //have found cell with next gVCF position, handle accumulated values
      {
        next_start_position = next_coord[1];
        assert(next_coord[1] > current_start_position);
        handle_gvcf_ranges(end_pq, PQ_end_vec, query_config, gt_column, tile_idx_2_iters, output_stream,
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
      curr_struct.m_end_point = static_cast<const AttributeTile<int64_t>& >(*tile_its[END_query_idx]).cell(cell_it.pos());
      assert(curr_struct.m_end_point >= current_start_position);
      //Store tile idx
      curr_struct.m_tile_idx = tile_idx;
      ++(current_iterator_tracker.m_reference_counter);
      //Store position of cell wrt tile
      curr_struct.m_cell_pos = cell_it.pos();
      assert(cell_it.pos() < (*tile_its[COORDS_query_idx]).cell_num());
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
  handle_gvcf_ranges(end_pq, PQ_end_vec, query_config, gt_column, tile_idx_2_iters, output_stream,
      current_start_position, 0, true);
  delete[] tile_its;
  delete gt_column;
}

void VariantQueryProcessor::do_query_bookkeeping(const StorageManager::ArrayDescriptor* array_descriptor,
    VariantQueryConfig& query_config)
{
  QueryProcessor::obtain_TileDB_attribute_idxs(array_descriptor, query_config);
  //Add COORDS as a query attribute by default
  unsigned COORDS_schema_idx = 
          m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(COORDS_schema_idx));
  query_config.add_attribute_to_query("COORDINATES", COORDS_schema_idx);
  //Add END as a query attribute by default
  unsigned END_schema_idx = 
          m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_END_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(END_schema_idx));
  query_config.add_attribute_to_query("END", END_schema_idx);
  //Check if OFFSETS and NULL need to be added as part of queried attributes
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  unsigned OFFSETS_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_OFFSETS_IDX);
  unsigned NULL_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_NULL_IDX);
  for(auto i=0u;i<num_queried_attributes;++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    unsigned schema_idx =  query_config.get_schema_idx_for_query_idx(i);
    unsigned known_variant_field_enum = 
      m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(schema_idx);
    //known field
    if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(known_variant_field_enum))
    {
      //Does the field require OFFSETS?
      if(uses_OFFSETS_field(known_variant_field_enum))
      {
        assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(OFFSETS_schema_idx));
        query_config.add_attribute_to_query("OFFSETS", OFFSETS_schema_idx);
      }
      //Does the field require NULL?
      if(is_NULL_bitidx_defined_for_known_field_enum(known_variant_field_enum))
      {
        assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(NULL_schema_idx));
        query_config.add_attribute_to_query("NULL", NULL_schema_idx);
      }
    }
  }
  //Set enum within query
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
      assert(VariantQueryProcessor::m_known_variant_field_names[known_variant_field_enum] == query_config.get_query_attribute_name(i));
    }
  }
  query_config.set_done_bookkeeping(true);
}

GTColumn* VariantQueryProcessor::gt_get_column(
    const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config,
    GTProfileStats* stats) const {
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t row_num = dim_domains[0].second - dim_domains[0].first + 1;
  unsigned int attribute_num = array_schema.attribute_num();

  assert(query_config.get_num_column_ranges() > 0);
  assert(query_config.is_bookkeeping_done() && "Before calling gt_get_column(), do_query_bookkeeping() function must be called");
  unsigned column_range_idx = 0u;
  uint64_t col = query_config.get_column_range(column_range_idx).first;

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
      gt_initialize_tile_its(ad, query_config, column_range_idx, tile_its, tile_it_end);
  //Get query idx for COORDS
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);

  // Create and initialize GenotypingColumn members
  GTColumn* gt_column = new GTColumn(col, row_num);

  // Create cell iterators
  Tile::const_reverse_iterator cell_it, cell_it_end;
  uint64_t num_deref_tile_iters = 0;
#ifdef DO_PROFILING
  uint64_t num_cells_touched = 0;
  uint64_t num_tiles_touched = 0;
  bool first_sample = true;
#endif
  // Fill the genotyping column
  while(tile_its[COORDS_query_idx] != tile_it_end && filled_rows < row_num) {
    // Initialize cell iterators for the coordinates
    cell_it = (*tile_its[COORDS_query_idx]).rbegin();
    cell_it_end = (*tile_its[COORDS_query_idx]).rend();
#ifdef DO_PROFILING
    num_deref_tile_iters += 2;	//why 2, cell_it and cell_it_end
    ++num_tiles_touched;
#endif
    while(cell_it != cell_it_end && filled_rows < row_num) {
      std::vector<int64_t> next_coord = *cell_it;
#ifdef DO_PROFILING
      ++num_cells_touched;
#endif
      // If next cell is not on the right of col, and corresponds to 
      // uninvestigated row
      if(next_coord[1] <= col && CHECK_UNINITIALIZED_SAMPLE_GIVEN_REF(gt_column->REF_[next_coord[0]])) {
        gt_fill_row<StorageManager::const_reverse_iterator>(gt_column, next_coord[0], next_coord[1], cell_it.pos(), query_config,
            tile_its, &num_deref_tile_iters);
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
  stats->m_sum_num_deref_tile_iters += num_deref_tile_iters;
  stats->m_sum_sq_num_deref_tile_iters += (num_deref_tile_iters*num_deref_tile_iters);
  stats->m_sum_num_tiles_touched += num_tiles_touched;
  stats->m_sum_sq_num_tiles_touched += (num_tiles_touched*num_tiles_touched);
#endif

  return gt_column;
}

template<class ITER>
void VariantQueryProcessor::gt_fill_row(
    GTColumn* gt_column, int64_t row, int64_t column, int64_t pos,
    const VariantQueryConfig& query_config,
    const ITER* tile_its, uint64_t* num_deref_tile_iters) const {
  // First check if the row is NULL
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  int64_t END_v = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX)]).cell(pos);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  if(END_v < gt_column->col_) {
    gt_column->REF_[row] = "$";
    return;
  }

  // Retrieve the offsets
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_OFFSETS_IDX));
  const AttributeTile<int64_t>& OFFSETS_tile = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_OFFSETS_IDX)]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  int64_t REF_offset = OFFSETS_tile.cell(pos*5);
  int64_t ALT_offset = OFFSETS_tile.cell(pos*5+1);
  int64_t PL_offset = OFFSETS_tile.cell(pos*5+4);

  // Retrieve the NULL bitmap
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_NULL_IDX));
  const AttributeTile<int>& NULL_tile = 
      static_cast<const AttributeTile<int>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_NULL_IDX)]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  int NULL_bitmap = NULL_tile.cell(pos);

  char c;
  int i;

  // Fill the REF
  //If the queried column is identical to the cell's column, then the REF value stored is correct
  //Else, the cell stores an interval and the REF value is set to "N" which means could be anything
  if(column == gt_column->col_)
  {
    assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_REF_IDX));
    const AttributeTile<char>& REF_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_REF_IDX)]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
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
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_ALT_IDX));
  const AttributeTile<char>& ALT_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_ALT_IDX)]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
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
  // NULL IDX is the last IDX after all the attributes so shifting from that offset
  assert(is_NULL_bitidx_defined_for_known_field_enum(GVCF_PL_IDX));
  unsigned NULL_PL_bitidx = get_NULL_bitidx_for_known_field_enum(GVCF_PL_IDX);
  if(((NULL_bitmap >> NULL_PL_bitidx)  & 1) == 0) { // If the PL values are
    assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_PL_IDX));
    const AttributeTile<int>& PL_tile = 
        static_cast<const AttributeTile<int>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_PL_IDX)]);
#ifdef DO_PROFILING
    ++(*num_deref_tile_iters);
#endif
    int ALT_num = gt_column->ALT_[row].size(); 
    int PL_num = (ALT_num+1)*(ALT_num+2)/2;
    for(int i=0; i<PL_num; i++) 
      gt_column->PL_[row].push_back(PL_tile.cell(PL_offset+i));
  }

  if(m_GT_schema_version >= GT_SCHEMA_V1)
  {
    // Fill AF, AN, and AC
    if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_AF_IDX))
      fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, query_config.get_query_idx_for_known_field_enum(GVCF_AF_IDX),
          &(gt_column->AF_[row]));
    if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_AN_IDX))
      fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, query_config.get_query_idx_for_known_field_enum(GVCF_AN_IDX),
          &(gt_column->AN_[row]));
    if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_AC_IDX))
      fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, query_config.get_query_idx_for_known_field_enum(GVCF_AC_IDX),
          &(gt_column->AC_[row]));
  }
}

inline
unsigned int VariantQueryProcessor::gt_initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    const VariantQueryConfig& query_config, const unsigned col_range_idx,
    StorageManager::const_reverse_iterator*& tile_its, 
    StorageManager::const_reverse_iterator& tile_it_end) const {
  //Num attributes in query
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  // Create reverse iterators
  tile_its = new StorageManager::const_reverse_iterator[num_queried_attributes];
  // Find the rank of the tile the left sweep starts from.
  auto start_rank = get_storage_manager().get_left_sweep_start_rank(ad, query_config.get_column_range(col_range_idx).first);
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    tile_its[i] = get_storage_manager().rbegin(ad, schema_idx, start_rank);
  }
  tile_it_end = get_storage_manager().rend(ad, m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  return num_queried_attributes - 1;
}

// Helper function to fill the attribute given the tile pointer,
// position, and index
template<class ITER>
void VariantQueryProcessor::fill_cell_attribute(const int64_t& pos, const ITER* tile_its, 
        uint64_t* num_deref_tile_iters,
        const unsigned IDX, int *p_int_v) const {
    // Retrieve the tile corresponding to the IDX 
    const AttributeTile<int>& m_attribute_tile =
        static_cast<const AttributeTile<int>& >(*tile_its[IDX]);
#ifdef DO_PROFILING
    ++(*num_deref_tile_iters);
#endif
    *p_int_v = m_attribute_tile.cell(pos);
}

// Override of the function above for float type
template<class ITER>
void VariantQueryProcessor::fill_cell_attribute(int64_t pos, const ITER* tile_its, 
        uint64_t* num_deref_tile_iters,
        unsigned IDX, float *p_float_v) const {
    // Retrieve the tile corresponding to the IDX 
    const AttributeTile<float>& m_tile = 
        static_cast<const AttributeTile<float>& >(*tile_its[IDX]);
#ifdef DO_PROFILING
    ++(*num_deref_tile_iters);
#endif
    *p_float_v = m_tile.cell(pos);
}

void VariantQueryProcessor::clear()
{
  m_NULL_bit_enum_idx_vec.clear();
  m_known_field_enum_uses_offset.clear();
  m_schema_idx_to_known_variant_field_enum_LUT.reset_luts();
}
