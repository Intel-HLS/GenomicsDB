#include "query_variants.h"
#include "variant_operations.h"

void VariantQueryProcessor::initialize_version(const StorageManager::ArrayDescriptor* ad)
{
  const auto& schema = ad->array_schema();
  //Initialize to V1 schema by default
  m_GT_schema_version = GT_SCHEMA_V1;
  m_PL_NULL_bitidx = 0u; 
  //Check if any attributes in V2 schema
  const auto v2_fields = std::unordered_set<std::string>{ "AN", "AF", "AC" };
  for(auto i=0u;i<schema.attribute_num();++i)
    if(v2_fields.find(schema.attribute_name(i)) != v2_fields.end())
    {
      m_GT_schema_version = GT_SCHEMA_V2;
      break;
    }
  //Last queried attribute, keep incrementing
  GT_IDX_COUNT = 0u;
  GT_END_IDX = GT_IDX_COUNT++;
  GT_REF_IDX = GT_IDX_COUNT++;
  GT_ALT_IDX = GT_IDX_COUNT++;
  GT_PL_IDX = GT_IDX_COUNT++;
  if(m_GT_schema_version >= GT_SCHEMA_V2)
  {
    GT_AF_IDX = GT_IDX_COUNT++;
    GT_AN_IDX = GT_IDX_COUNT++;
    GT_AC_IDX = GT_IDX_COUNT++;
  }
  GT_NULL_IDX = GT_IDX_COUNT++;
  GT_OFFSETS_IDX = GT_IDX_COUNT++;
  GT_COORDINATES_IDX = GT_IDX_COUNT++;
  //Last attribute, keep incrementing
  GVCF_COORDINATES_IDX = 0u;
  GVCF_END_IDX = GVCF_COORDINATES_IDX++;
  GVCF_REF_IDX = GVCF_COORDINATES_IDX++;
  GVCF_ALT_IDX = GVCF_COORDINATES_IDX++;
  GVCF_QUAL_IDX = GVCF_COORDINATES_IDX++;
  GVCF_FILTER_IDX = GVCF_COORDINATES_IDX++;
  GVCF_BASEQRANKSUM_IDX = GVCF_COORDINATES_IDX++;
  GVCF_CLIPPINGRANKSUM_IDX = GVCF_COORDINATES_IDX++;
  GVCF_MQRANKSUM_IDX = GVCF_COORDINATES_IDX++;
  GVCF_READPOSRANKSUM_IDX = GVCF_COORDINATES_IDX++;
  GVCF_DP_IDX = GVCF_COORDINATES_IDX++;
  GVCF_MQ_IDX = GVCF_COORDINATES_IDX++;
  GVCF_MQ0_IDX = GVCF_COORDINATES_IDX++;
  GVCF_DP_FMT_IDX = GVCF_COORDINATES_IDX++;
  GVCF_MIN_DP_IDX = GVCF_COORDINATES_IDX++;
  GVCF_GQ_IDX = GVCF_COORDINATES_IDX++;
  GVCF_SB_1_IDX = GVCF_COORDINATES_IDX++;
  GVCF_SB_2_IDX = GVCF_COORDINATES_IDX++;
  GVCF_SB_3_IDX = GVCF_COORDINATES_IDX++;
  GVCF_SB_4_IDX = GVCF_COORDINATES_IDX++;
  GVCF_AD_IDX = GVCF_COORDINATES_IDX++;
  GVCF_PL_IDX = GVCF_COORDINATES_IDX++;
  if(m_GT_schema_version >= GT_SCHEMA_V2)
  {
    GVCF_AF_IDX = GVCF_COORDINATES_IDX++;
    GVCF_AN_IDX = GVCF_COORDINATES_IDX++;
    GVCF_AC_IDX = GVCF_COORDINATES_IDX++;
  }
  GVCF_NULL_IDX = GVCF_COORDINATES_IDX++;
  GVCF_OFFSETS_IDX = GVCF_COORDINATES_IDX++;
  //Change bitidx in NULL field for PL
  if(m_GT_schema_version >= GT_SCHEMA_V2)
    m_PL_NULL_bitidx = ((GVCF_NULL_IDX - 1) - GVCF_PL_IDX);
}

void VariantQueryProcessor::handle_gvcf_ranges(VariantIntervalPQ& end_pq, std::vector<PQStruct>& PQ_end_vec, GTColumn* gt_column,
      std::unordered_map<uint64_t, GTTileIteratorsTracker>& tile_idx_2_iters, std::ostream& output_stream,
    int64_t current_start_position, int64_t next_start_position, bool is_last_call)
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
        gt_fill_row<StorageManager::const_iterator>(gt_column, i, curr_struct.m_array_column, curr_struct.m_cell_pos, 
	    &((*find_iter).second.m_iter_vector[0]), &num_deref_tile_iters, false);
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

void VariantQueryProcessor::iterate_over_all_cells(const StorageManager::ArrayDescriptor* ad)
{
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
  
  unsigned num_queried_attributes = GT_COORDINATES_IDX + 1;
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  std::unordered_map<uint64_t, GTTileIteratorsTracker> tile_idx_2_iters;
  // END
  tile_its[GT_END_IDX] = get_storage_manager().begin(ad, GVCF_END_IDX);
  // REF
  tile_its[GT_REF_IDX] = get_storage_manager().begin(ad, GVCF_REF_IDX);
  // ALT
  tile_its[GT_ALT_IDX] = get_storage_manager().begin(ad, GVCF_ALT_IDX);
  // PL
  tile_its[GT_PL_IDX] = get_storage_manager().begin(ad, GVCF_PL_IDX);
  // NULL
  tile_its[GT_NULL_IDX] = get_storage_manager().begin(ad, GVCF_NULL_IDX);
  // OFFSETS
  tile_its[GT_OFFSETS_IDX] = get_storage_manager().begin(ad, GVCF_OFFSETS_IDX);
  // coordinates
  tile_its[GT_COORDINATES_IDX] = get_storage_manager().begin(ad, GVCF_COORDINATES_IDX);
  
  StorageManager::const_iterator tile_it_end = get_storage_manager().end(ad, GVCF_COORDINATES_IDX);
  GTColumn* gt_column = new GTColumn(0, total_num_samples);
  uint64_t num_deref_tile_iters = 0; 
  for(;tile_its[num_queried_attributes-1] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))
  {
    Tile::const_iterator cell_it = (*tile_its[num_queried_attributes-1]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[num_queried_attributes-1]).end();
    std::vector<int64_t> next_coord = *cell_it;
    gt_column->reset();
    gt_column->col_ = next_coord[1];
    gt_fill_row<StorageManager::const_iterator>(gt_column, next_coord[0], next_coord[1], cell_it.pos(), tile_its,
	&num_deref_tile_iters, false);
  }
}

void VariantQueryProcessor::scan_and_operate(const StorageManager::ArrayDescriptor* ad, std::ostream& output_stream)
{
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t total_num_samples = dim_domains[0].second - dim_domains[0].first + 1;
  
  unsigned num_queried_attributes = GT_COORDINATES_IDX + 1;
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  std::unordered_map<uint64_t, GTTileIteratorsTracker> tile_idx_2_iters;
  // END
  tile_its[GT_END_IDX] = get_storage_manager().begin(ad, GVCF_END_IDX);
  // REF
  tile_its[GT_REF_IDX] = get_storage_manager().begin(ad, GVCF_REF_IDX);
  // ALT
  tile_its[GT_ALT_IDX] = get_storage_manager().begin(ad, GVCF_ALT_IDX);
  // PL
  tile_its[GT_PL_IDX] = get_storage_manager().begin(ad, GVCF_PL_IDX);
  // NULL
  tile_its[GT_NULL_IDX] = get_storage_manager().begin(ad, GVCF_NULL_IDX);
  // OFFSETS
  tile_its[GT_OFFSETS_IDX] = get_storage_manager().begin(ad, GVCF_OFFSETS_IDX);
  // coordinates
  tile_its[GT_COORDINATES_IDX] = get_storage_manager().begin(ad, GVCF_COORDINATES_IDX);
  
  StorageManager::const_iterator tile_it_end = get_storage_manager().end(ad, GVCF_COORDINATES_IDX);
  
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
      curr_struct.m_end_point = static_cast<const AttributeTile<int64_t>& >(*tile_its[GT_END_IDX]).cell(cell_it.pos());
      assert(curr_struct.m_end_point >= current_start_position);
      //Store tile idx
      curr_struct.m_tile_idx = tile_idx;
      ++(current_iterator_tracker.m_reference_counter);
      //Store position of cell wrt tile
      curr_struct.m_cell_pos = cell_it.pos();
      assert(cell_it.pos() < (*tile_its[GT_COORDINATES_IDX]).cell_num());
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

GTColumn* VariantQueryProcessor::gt_get_column(
    const StorageManager::ArrayDescriptor* ad, uint64_t col, GTProfileStats* stats) const {
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
  while(tile_its[gt_attribute_num] != tile_it_end && filled_rows < row_num) {
    // Initialize cell iterators for the coordinates
    cell_it = (*tile_its[gt_attribute_num]).rbegin();
    cell_it_end = (*tile_its[gt_attribute_num]).rend();
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
        gt_fill_row<StorageManager::const_reverse_iterator>(gt_column, next_coord[0], next_coord[1], cell_it.pos(), tile_its, &num_deref_tile_iters,
                true);
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
    const ITER* tile_its, uint64_t* num_deref_tile_iters, bool load_v2_attributes) const {
  // First check if the row is NULL
  int64_t END_v = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[GT_END_IDX]).cell(pos);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  if(END_v < gt_column->col_) {
    gt_column->REF_[row] = "$";
    return;
  }

  // Retrieve the offsets
  const AttributeTile<int64_t>& OFFSETS_tile = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[GT_OFFSETS_IDX]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  int64_t REF_offset = OFFSETS_tile.cell(pos*5);
  int64_t ALT_offset = OFFSETS_tile.cell(pos*5+1);
  int64_t PL_offset = OFFSETS_tile.cell(pos*5+4);

  // Retrieve the NULL bitmap
  const AttributeTile<int>& NULL_tile = 
      static_cast<const AttributeTile<int>& >(*tile_its[GT_NULL_IDX]);
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
    const AttributeTile<char>& REF_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[GT_REF_IDX]);
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
  const AttributeTile<char>& ALT_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[GT_ALT_IDX]);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  i = 0;
  std::string ALT_s = "";
  while(1) {    //exit if empty ALT_s string or '&' seen
    c = ALT_tile.cell(ALT_offset+i);
    if(c == '\0') {
      if(ALT_s.length() == 0u) //end of ALT list, last string is empty
        break;
      gt_column->ALT_[row].push_back(ALT_s);
      ALT_s = "";
    } else {
      if(c == '&')
      {
        gt_column->ALT_[row].push_back("&");
        break;
      }
      else
        ALT_s.push_back(c);
    }
    i++;
  }
  assert(ALT_s == "");

  // Fill the PL values
  // NULL IDX is the last IDX after all the attributes so shifting from that offset
  if(((NULL_bitmap >> m_PL_NULL_bitidx)  & 1) == 0) { // If the PL values are
    const AttributeTile<int>& PL_tile = 
        static_cast<const AttributeTile<int>& >(*tile_its[GT_PL_IDX]);
#ifdef DO_PROFILING
    ++(*num_deref_tile_iters);
#endif
    int ALT_num = gt_column->ALT_[row].size(); 
    int PL_num = (ALT_num+1)*(ALT_num+2)/2;
    for(int i=0; i<PL_num; i++) 
      gt_column->PL_[row].push_back(PL_tile.cell(PL_offset+i));
  }

  if(m_GT_schema_version >= GT_SCHEMA_V2 && load_v2_attributes)
  {
    // Fill AF, AN, and AC
    fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, GT_AF_IDX, &(gt_column->AF_[row]));
    fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, GT_AN_IDX, &(gt_column->AN_[row]));
    fill_cell_attribute<ITER>(pos, tile_its, num_deref_tile_iters, GT_AC_IDX, &(gt_column->AC_[row]));
  }
}

inline
unsigned int VariantQueryProcessor::gt_initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    StorageManager::const_reverse_iterator*& tile_its, 
    StorageManager::const_reverse_iterator& tile_it_end,
    uint64_t col) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Create reverse iterators
  tile_its = new StorageManager::const_reverse_iterator[GT_IDX_COUNT];
  // Find the rank of the tile the left sweep starts from.
  auto start_rank = get_storage_manager().get_left_sweep_start_rank(ad, col);

  // END
  tile_its[GT_END_IDX] = get_storage_manager().rbegin(ad, GVCF_END_IDX, start_rank);
  // REF
  tile_its[GT_REF_IDX] = get_storage_manager().rbegin(ad, GVCF_REF_IDX, start_rank);
  // ALT
  tile_its[GT_ALT_IDX] = get_storage_manager().rbegin(ad, GVCF_ALT_IDX, start_rank);
  // PL
  tile_its[GT_PL_IDX] = get_storage_manager().rbegin(ad, GVCF_PL_IDX, start_rank);
  if(m_GT_schema_version >= GT_SCHEMA_V2)
  {
    // AF
    tile_its[GT_AF_IDX] = get_storage_manager().rbegin(ad, GVCF_AF_IDX, start_rank);
    // AN
    tile_its[GT_AN_IDX] = get_storage_manager().rbegin(ad, GVCF_AN_IDX, start_rank);
    // AC
    tile_its[GT_AC_IDX] = get_storage_manager().rbegin(ad, GVCF_AC_IDX, start_rank);
  }
  // NULL
  tile_its[GT_NULL_IDX] = get_storage_manager().rbegin(ad, GVCF_NULL_IDX, start_rank);
  // OFFSETS
  tile_its[GT_OFFSETS_IDX] = get_storage_manager().rbegin(ad, GVCF_OFFSETS_IDX, start_rank);
  // coordinates
  tile_its[GT_COORDINATES_IDX] = get_storage_manager().rbegin(ad, GVCF_COORDINATES_IDX, start_rank);
  tile_it_end = get_storage_manager().rend(ad, attribute_num);

  // The number of attributes is 6, and the coordinates is the extra one
  return GT_COORDINATES_IDX;
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

