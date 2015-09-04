#include "gt_common.h"
#include "query_variants.h"

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

unsigned KnownFieldInfo::get_num_elements_for_known_field_enum(unsigned num_ALT_alleles, unsigned ploidy) const
{
  unsigned length = 0u;
  unsigned num_alleles = num_ALT_alleles + 1u;
  switch(m_length_descriptor)
  {
    case BCF_VL_FIXED:
      length = m_num_elements;
      break;
    case BCF_VL_VAR:
      length = 1u;      //function that reads from tile will know what to do
      break;
    case BCF_VL_A:
      length = num_ALT_alleles;
      break;
    case BCF_VL_R:
      length = num_alleles;
      break;
    case BCF_VL_G:
      length = (num_alleles*(num_alleles+1))/2;
      break;
    case BCF_VL_P:
      length = ploidy;
      break;
    default:
      cerr << "Unknown length descriptor "<<m_length_descriptor<<" - ignoring\n";
      break;
  }
  return length;
}
//Static members
bool VariantQueryProcessor::m_are_static_members_initialized = false;
vector<string> VariantQueryProcessor::m_known_variant_field_names = vector<string>{
    "END",
    "REF",
    "ALT",
    "QUAL",
    "FILTER_ID",
    "BaseQRankSum",
    "ClippingRankSum",
    "MQRankSum",
    "ReadPosRankSum",
    "DP",
    "MQ",
    "MQ0",
    "DP_FMT",
    "MIN_DP",
    "GQ",
    "SB",
    "AD",
    "PL",
    "AF",
    "AN",
    "AC",
    "GT",
    "PS"
};
unordered_map<string, unsigned> VariantQueryProcessor::m_known_variant_field_name_to_enum;
unordered_map<type_index, shared_ptr<VariantFieldCreatorBase>> VariantQueryProcessor::m_type_index_to_creator;

//Initialize static members function
void VariantQueryProcessor::initialize_static_members()
{
  for(auto i=0u;i<VariantQueryProcessor::m_known_variant_field_names.size();++i)
    VariantQueryProcessor::m_known_variant_field_name_to_enum[VariantQueryProcessor::m_known_variant_field_names[i]] = i;
  VariantQueryProcessor::m_type_index_to_creator.clear();
  //Map type_index to creator functions
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

void VariantQueryProcessor::initialize_known(const ArraySchema& schema)
{
  //Ploidy requirements
  m_known_field_enum_to_info[GVCF_GT_IDX].m_ploidy_required = true;
  //Initialize schema idx <--> known_field_enum mapping
  m_schema_idx_to_known_variant_field_enum_LUT.resize_luts_if_needed(schema.attribute_num(), GVCF_NUM_KNOWN_FIELDS);
  for(auto i=0u;i<schema.attribute_num();++i)
  {
    auto iter = VariantQueryProcessor::m_known_variant_field_name_to_enum.find(schema.attribute_name(i));
    if(iter != VariantQueryProcessor::m_known_variant_field_name_to_enum.end())
      m_schema_idx_to_known_variant_field_enum_LUT.add_schema_idx_known_field_mapping(i, (*iter).second);
  }
}

void VariantQueryProcessor::initialize_v0(const ArraySchema& schema)
{
  m_GT_schema_version = GT_SCHEMA_V0;
}

//Added AF,AN and AC fields
void VariantQueryProcessor::initialize_v1(const ArraySchema& schema)
{
  //Check if any attributes in V2 schema
  const auto v1_fields = std::unordered_set<std::string>{ "AF", "AN", "AC" };
  for(auto i=0u;i<schema.attribute_num();++i)
    if(v1_fields.find(schema.attribute_name(i)) != v1_fields.end())
    {
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V1;
      break;
    }
}

//Added GT and PS fields
void VariantQueryProcessor::initialize_v2(const ArraySchema& schema)
{
  //Check if any attributes in V2 schema
  const auto v2_fields = std::unordered_set<std::string>{ "GT", "PS" };
  for(auto i=0u;i<schema.attribute_num();++i)
    if(v2_fields.find(schema.attribute_name(i)) != v2_fields.end())
    {
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V2;
      break;
    }
}

void VariantQueryProcessor::initialize_length_descriptor(unsigned idx)
{
  switch(idx)
  {
    case GVCF_REF_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      break;
    case GVCF_ALT_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      m_known_field_enum_to_info[idx].m_field_creator = std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldALTData>());
      break;
    case GVCF_FILTER_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      //TODO: How to get this?
      break;
    case GVCF_AD_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_R;
      break;
    case GVCF_PL_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_G;
      break;
    case GVCF_GT_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_P;
      break;
    case GVCF_SB_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      m_known_field_enum_to_info[idx].m_num_elements = 4u;
      break;
    default:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      m_known_field_enum_to_info[idx].m_num_elements = 1u;
      break;
  }
}

void VariantQueryProcessor::initialize_version(const ArraySchema& schema)
{
  initialize_known(schema);
  //Initialize to v0 schema by default
  initialize_v0(schema);
  initialize_v1(schema);
  initialize_v2(schema);
}

void VariantQueryProcessor::register_field_creators(const ArraySchema& schema)
{
  m_field_factory.resize(schema.attribute_num());
  for(auto i=0u;i<schema.attribute_num();++i)
  {
    const type_info* curr_type_info = (schema.type(i));
    type_index t = type_index(*curr_type_info);
    auto iter = VariantQueryProcessor::m_type_index_to_creator.find(t);
    if(iter == VariantQueryProcessor::m_type_index_to_creator.end())
      throw UnknownAttributeTypeException("Unknown type of schema attribute "+std::string(curr_type_info->name()));
    //For known fields, check for special creators
    unsigned enumIdx = m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(i);
    if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(enumIdx) && requires_special_creator(enumIdx))
      m_field_factory.Register(i, m_known_field_enum_to_info[enumIdx].m_field_creator);
    else
      m_field_factory.Register(i, (*iter).second);
  }
}

VariantQueryProcessor::VariantQueryProcessor(StorageManager* storage_manager, const std::string& array_name)
: QueryProcessor(storage_manager)
{
  //initialize static members
  if(!VariantQueryProcessor::m_are_static_members_initialized)
    VariantQueryProcessor::initialize_static_members();
  clear();
  m_ad = storage_manager->open_array(array_name, "r");
  auto status = storage_manager->get_array_schema(m_ad, m_array_schema);
  assert(status == TILEDB_OK);
  //Mapping from known_field enum
  m_known_field_enum_to_info.resize(GVCF_NUM_KNOWN_FIELDS);
  //Initialize versioning information
  initialize_version(*m_array_schema);
  //set length descriptors and creator objects for special attributes
  for(auto i=0u;i<m_known_field_enum_to_info.size();++i)
    initialize_length_descriptor(i);
  //Register creators in factory
  register_field_creators(*m_array_schema);
}

void VariantQueryProcessor::handle_gvcf_ranges(VariantCallEndPQ& end_pq,
    const VariantQueryConfig& query_config, Variant& variant,
    SingleVariantOperatorBase& variant_operator,
    int64_t current_start_position, int64_t next_start_position, bool is_last_call) const
{
  while(!end_pq.empty() && (current_start_position < next_start_position || is_last_call))
  {
    int64_t top_end_pq = end_pq.top()->get_column_end();
    int64_t min_end_point = (is_last_call || (top_end_pq < (next_start_position - 1))) ? top_end_pq : (next_start_position-1);
    //Prepare variant for aligned column interval
    variant.set_column_interval(current_start_position, min_end_point);
    //Set REF to N if the call interval is split in the middle
    if(query_config.is_defined_query_idx_for_known_field_enum(GVCF_REF_IDX))
      for(Variant::valid_calls_iterator iter=variant.begin();iter!=variant.end();++iter)
      {
        auto& curr_call = *iter;
        assert(curr_call.get_column_begin() <= current_start_position);
        modify_reference_if_in_middle(curr_call, query_config, current_start_position);
      }
    variant_operator.operate(variant, query_config);
    //The following intervals have been completely processed
    while(!end_pq.empty() && end_pq.top()->get_column_end() == min_end_point)
    {
      auto top_element = end_pq.top();
      top_element->mark_valid(false);
      end_pq.pop();
    }
    current_start_position = min_end_point + 1;   //next start position, after the end
  }
}

#if 0
void VariantQueryProcessor::iterate_over_all_tiles(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config) const
{
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
  //Variant object 
  Variant variant(&query_config);
  variant.resize_based_on_query();
  uint64_t num_deref_tile_iters = 0; 
  for(;tile_its[COORDS_query_idx] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))
  {
    Tile::const_iterator cell_it = (*tile_its[COORDS_query_idx]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[COORDS_query_idx]).end();
    std::vector<int64_t> next_coord = *cell_it;
    variant.set_column_interval(next_coord[1], next_coord[1]);
    gt_fill_row<StorageManager::const_iterator>(variant, next_coord[0], next_coord[1], cell_it.pos(), query_config, tile_its,
	&num_deref_tile_iters);
  }
  delete[] tile_its;
}

void VariantQueryProcessor::scan_and_operate(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config,
    SingleVariantOperatorBase& variant_operator, unsigned column_interval_idx) const
{
  assert(query_config.is_bookkeeping_done());
  //Priority queue of VariantCalls ordered by END positions
  VariantCallEndPQ end_pq;
  //Rank of tile from which scan should start
  int64_t start_rank = 0;
  //Current gVCF position being operated on
  int64_t current_start_position = -1ll;
  //Variant object
  Variant variant(&query_config);
  variant.resize_based_on_query();
  //Scan only queried interval, not whole array
  if(query_config.get_num_column_intervals() > 0)
  {
    //If the queried interval is [100:200], then the first part of this function gets a Variant object
    //containing Calls that intersect with 100. Some of them could start before 100 and extend beyond. This
    //information is recorded in the Call
    //This part of the code accumulates such Calls, sets the current_start_position to query column interval begin
    //and lets the code in the for loop nest (forward scan) handle calling handle_gvcf_ranges()
    //Row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
    //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
    gt_get_column(ad, query_config, column_interval_idx, variant, nullptr);
    //Insert valid calls produced by gt_get_column into the priority queue
    for(Variant::valid_calls_iterator iter=variant.begin();iter != variant.end();++iter)
    {
      auto& curr_call = *iter;
      end_pq.push(&curr_call);
      assert(end_pq.size() <= query_config.get_num_rows_to_query());
    }
    //Valid calls were found, start position == query colum interval begin
    if(end_pq.size() > 0)
      current_start_position = query_config.get_column_begin(column_interval_idx);
    //For forward scan, get the tile with lowest rank that contains a column > query_column. All cells with column == query_column
    //(or intersecting with query_column) will be handled by gt_get_column(). Hence, must start from next column
    start_rank = get_storage_manager().get_right_sweep_start_rank(ad, query_config.get_column_begin(column_interval_idx)+1);
  }
  //Forward scan iterators
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    tile_its[i] = get_storage_manager().begin(ad, schema_idx, start_rank);
  }
  //Get tile end iter for COORDS
  StorageManager::const_iterator tile_it_end = get_storage_manager().end(ad,
      m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  //Get query idx for COORDS 
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  //Iterator that points to the first cell to be accessed in the forward scan
  Tile::const_iterator first_cell_it;
  //Get first position in the array from which forward scan should start 
  if(tile_its[COORDS_query_idx] != tile_it_end)
  {
    //If a specific interval is being queried, then the first cell that should be accessed is the one
    //whose column is > queried column interval begin, since gt_get_column() would take care of all cells 
    //with column == queried column interval begin (or intersecting with)
    first_cell_it = (query_config.get_num_column_intervals() > 0) ? 
      get_first_cell_after((*tile_its[COORDS_query_idx]), query_config.get_column_begin(column_interval_idx))
      : (*(tile_its[COORDS_query_idx])).begin();
    Tile::const_iterator cell_it_end = (*(tile_its[COORDS_query_idx])).end();
    //Initialize current_start_position only if gt_get_column() could not find any valid Calls intersecting
    //with queried column interval begin
    if(first_cell_it != cell_it_end && current_start_position < 0)
    {
      std::vector<int64_t> current_coord = *first_cell_it;
      current_start_position = current_coord[1];
    }
  }
  //Set current column for variant (end is un-important as Calls are used to track end of intervals)
  variant.set_column_interval(current_start_position, current_start_position);
  //Next co-ordinate to consider
  int64_t next_start_position = -1ll;
  uint64_t tile_idx = 0ull;
  uint64_t num_deref_tile_iters = 0ull; //profiling
  bool first_tile = true;
  bool break_out = false;
  for(;tile_its[COORDS_query_idx] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))  //why -1, advance uses i<=N for loop
  {
    // Initialize cell iterators for the coordinates
    //For the first tile, start at cell after the interval begin position, since gt_get_column would have
    //handled everything <= begin
    Tile::const_iterator cell_it = first_tile ? first_cell_it : (*tile_its[COORDS_query_idx]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[COORDS_query_idx]).end();
    for(;cell_it != cell_it_end;++cell_it) {
      std::vector<int64_t> next_coord = *cell_it;
      //If only interval requested and end of interval crossed, exit
      if(query_config.get_num_column_intervals() > 0 && next_coord[1] > query_config.get_column_end(column_interval_idx))
      {
        break_out = true;
        break;
      }
      if(next_coord[1] != current_start_position) //have found cell with next gVCF position, handle accumulated values
      {
        next_start_position = next_coord[1];
        assert(next_coord[1] > current_start_position);
        handle_gvcf_ranges(end_pq, query_config, variant, variant_operator, current_start_position,
            next_start_position, false);
        assert(end_pq.empty() || end_pq.top()->get_column_end() >= next_start_position);  //invariant
        //Set new start for next interval
        current_start_position = next_start_position;
        variant.set_column_interval(current_start_position, current_start_position);
        //Do not reset variant as some of the Calls that are long intervals might still be valid 
      }
      //Accumulate cells with position == current_start_position
      //Include only if row is part of query
      if(query_config.is_queried_array_row_idx(next_coord[0]))
      {
        gt_fill_row<StorageManager::const_iterator>(variant, next_coord[0], next_coord[1], cell_it.pos(), query_config,
            tile_its, &num_deref_tile_iters);
        auto& curr_call = variant.get_call(query_config.get_query_row_idx_for_array_row_idx(next_coord[0]));
        end_pq.push(&curr_call);
        assert(end_pq.size() <= query_config.get_num_rows_to_query());
      }
    }
    ++tile_idx;
    first_tile = false;
    if(break_out)
      break;
  }
  //handle last interval
  handle_gvcf_ranges(end_pq, query_config, variant, variant_operator, current_start_position, 0, true);
  delete[] tile_its;
}
#endif

void VariantQueryProcessor::do_query_bookkeeping(const ArraySchema& array_schema,
    VariantQueryConfig& query_config) const
{
  QueryProcessor::obtain_TileDB_attribute_idxs(array_schema, query_config);
  //Add END as a query attribute by default
  unsigned END_schema_idx = 
          m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_END_IDX);
  assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(END_schema_idx));
  query_config.add_attribute_to_query("END", END_schema_idx);
  //Check if ALT needs to be added as part of queried attributes
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  unsigned ALT_schema_idx =
    m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_ALT_IDX);
  for(auto i=0u;i<num_queried_attributes;++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    unsigned schema_idx =  query_config.get_schema_idx_for_query_idx(i);
    unsigned known_variant_field_enum = 
      m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(schema_idx);
    //known field
    if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(known_variant_field_enum))
    {
      //Does the length of the field depend on the number of alleles
      if(is_length_allele_dependent(known_variant_field_enum))
      {
        assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(ALT_schema_idx));
        query_config.add_attribute_to_query("ALT", ALT_schema_idx);
      }
    }
  }
  //Re-order query fields so that special fields are first
  query_config.reorder_query_fields();
  //Resize info vector
  query_config.resize_known_field_info_vector();
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
      query_config.set_info_for_query_idx(i, &(m_known_field_enum_to_info[known_variant_field_enum]));
    }
  }
  //Set number of rows in the array
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t row_num = dim_domains[0].second - dim_domains[0].first;
  query_config.set_num_rows_in_array(row_num);
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
    vector<Variant>& variants, GTProfileStats* stats) const {
  uint64_t num_deref_tile_iters = 0ull;
  uint64_t start_variant_idx = variants.size();
  //Will be used later in the function to produce Variants with one CallSet
  VariantQueryConfig subset_query_config(query_config);
  vector<int64_t> subset_rows = vector<int64_t>(1u, 0);
  subset_query_config.update_rows_to_query(subset_rows);  //only 1 row, row 0
  //If the queried interval is [100:200], then the first part of the function gets a Variant object
  //containing Calls that intersect with 100. Some of them could start before 100 and extend beyond. This
  //information is recorded in the Call
  Variant interval_begin_variant(&query_config);
  interval_begin_variant.resize_based_on_query();
  //Row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
  //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
  vector<uint64_t> query_row_idx_in_order = vector<uint64_t>(query_config.get_num_rows_to_query(), UNDEFINED_NUM_ROWS_VALUE);
  #if VERBOSE>0
    std::cout << "[query_variants:gt_get_column_interval] Getting " << query_config.get_num_rows_to_query() << " rows" << std::endl;
  #endif
  gt_get_column(ad, query_config, column_interval_idx, interval_begin_variant, stats, &query_row_idx_in_order);
  //This interval contains many Calls, likely un-aligned (no common start/end). Split this variant 
  //into multiple Variants, each containing a single call
  GA4GHCallInfoToVariantIdx call_info_2_variant;
  interval_begin_variant.move_calls_to_separate_variants(query_config, variants, query_row_idx_in_order,
      call_info_2_variant);
  //If this is not a single position query, need to fetch more cells
  if(query_config.get_column_end(column_interval_idx) >
      query_config.get_column_begin(column_interval_idx))
  {
    //Get the iterator to the first cell that has column > query_column_start. All cells with column == query_column_start
    //(or intersecting) would have been handled by gt_get_column(). Hence, must start from next column
    ArrayConstCellIterator<int64_t>* forward_iter = 0;
    gt_initialize_forward_iter(ad, query_config, query_config.get_column_interval(column_interval_idx).first+1, forward_iter);
    uint64_t num_deref_tile_iters = 0; 
    //Used to store single call variants  - one variant per cell
    //Multiple variants could be merged later on
    Variant tmp_variant(&subset_query_config);
    tmp_variant.resize_based_on_query();
    #if VERBOSE>0
      std::cout << "[query_variants:gt_get_column_interval] Fetching columns from " << query_config.get_column_begin(column_interval_idx) + 1;
      std::cout << " to " << query_config.get_column_end(column_interval_idx) << std::endl;
    #endif
    //Cell object that will be used for iterating over attributes
    Cell cell(m_array_schema, query_config.get_query_attributes_schema_idxs(), 0, true);
    for(;!(forward_iter->end());++(*forward_iter))
    {
      auto* cell_ptr = **forward_iter;
      //Coordinates are at the start of the cell
      auto next_coord = reinterpret_cast<const int64_t*>(cell_ptr);
      if(next_coord[1] > query_config.get_column_end(column_interval_idx))    //Genomic interval begins after end of query interval
        break;
      if(query_config.is_queried_array_row_idx(next_coord[0]))       //If row is part of query, process cell
      {
        //Create Variant with single Call (subset_query_config contains one row)
        subset_rows[0] = next_coord[0];
        subset_query_config.update_rows_to_query(subset_rows);
        tmp_variant.resize_based_on_query();
        assert(tmp_variant.get_num_calls() == 1u);      //exactly 1 call
        tmp_variant.get_call(0u).set_row_idx(next_coord[0]); //set row idx
        tmp_variant.set_column_interval(next_coord[1], next_coord[1]);
        //Set contents of cell 
        cell.set_cell(cell_ptr);
        gt_fill_row(tmp_variant, next_coord[0], next_coord[1], subset_query_config, cell, &num_deref_tile_iters);
        //Set correct end for the variant 
        assert(tmp_variant.get_num_calls() == 1u);      //exactly 1 call
        //Move call to variants vector, creating new Variant if necessary
        move_call_to_variant_vector(subset_query_config, tmp_variant.get_call(0), variants, call_info_2_variant);
      }
    }
    #if VERBOSE>0
      std::cout << "[query_variants:gt_get_column_interval] Fetching columns complete " << std::endl;
    #endif
    delete forward_iter;
  }
  #if VERBOSE>0
    std::cout << "[query_variants:gt_get_column_interval] re-arrangement of variants " << std::endl;
  #endif
  GA4GHOperator variant_operator;
  for(auto i=start_variant_idx;i<variants.size();++i)
    if(variants[i].get_num_calls() > 1u) //possible re-arrangement of PL/AD/GT fields needed
    {
      variant_operator.clear();
      variant_operator.operate(variants[i], query_config);
      assert(variant_operator.get_variants().size() == 1u);     //exactly one variant
      variants[i] = std::move(variant_operator.get_variants()[0]);
    }
  #if VERBOSE>0
    std::cout << "[query_variants:gt_get_column_interval] query complete " << std::endl;
  #endif
}

void VariantQueryProcessor::gt_get_column(
    const int ad,
    const VariantQueryConfig& query_config, unsigned column_interval_idx,
    Variant& variant, GTProfileStats* stats, std::vector<uint64_t>* query_row_idx_in_order) const {

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
    std::cout << "[query_variants:gt_get_column] Fetching column : " << col << std::endl;
  #endif
  // Initialize reverse tile iterators
  // The reverse iterator will start with the cells
  // of the various attributes that have the largest
  // id that either intersect with col, or precede col.
  ArrayConstReverseCellIterator<int64_t>* reverse_iter = 0;
  gt_initialize_reverse_iter(ad, query_config, query_config.get_column_interval(column_interval_idx).first, reverse_iter);
  uint64_t num_deref_tile_iters = 0;
#ifdef DO_PROFILING
  uint64_t num_cells_touched = 0;
  uint64_t num_tiles_touched = 0;
  bool first_sample = true;
#endif
  // Indicates how many rows have been filled.
  uint64_t filled_rows = 0;
  uint64_t num_valid_rows = 0;
  //Cell object that will be re-used
  Cell cell(m_array_schema, query_config.get_query_attributes_schema_idxs(), 0, true);
  // Fill the genotyping column
  while(!(reverse_iter->end()) && filled_rows < query_config.get_num_rows_to_query()) {
#ifdef DO_PROFILING
    //FIXME: incorrect #tiles
    num_deref_tile_iters += 2;	//why 2, cell_it and cell_it_end
    ++num_tiles_touched;
    ++num_cells_touched;
#endif
    auto* cell_ptr = **reverse_iter;
    //Coordinates are at the start of the cell
    auto next_coord = reinterpret_cast<const int64_t*>(cell_ptr);
    // If next cell is not on the right of col, and 
    // The rowIdx is being queried and
    // The row/call is uninitialized (uninvestigated) in the Variant
    if(next_coord[1] <= col && query_config.is_queried_array_row_idx(next_coord[0]))
    {
      auto curr_query_row_idx = query_config.get_query_row_idx_for_array_row_idx(next_coord[0]);
      auto& curr_call = variant.get_call(curr_query_row_idx);
      if(!(curr_call.is_initialized()))
      {
        //Set contents of cell 
        cell.set_cell(cell_ptr);
        gt_fill_row(variant, next_coord[0], next_coord[1], query_config, cell, &num_deref_tile_iters);
        ++filled_rows;
        if(curr_call.is_valid() && query_row_idx_in_order)
        {
          assert(num_valid_rows < query_row_idx_in_order->size());
          (*query_row_idx_in_order)[num_valid_rows++] = curr_query_row_idx;
        }
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
    }
    ++(*reverse_iter);
  }
  //Free memory 
  //if(cell.cell())
  //free(const_cast<void*>(cell.cell()));
  delete reverse_iter;

  if(query_row_idx_in_order)
    query_row_idx_in_order->resize(num_valid_rows);

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
}

unsigned VariantQueryProcessor::get_num_elements_for_known_field_enum(unsigned known_field_enum,
    unsigned num_ALT_alleles, unsigned ploidy) const
{
  assert(known_field_enum < m_known_field_enum_to_info.size());
  return m_known_field_enum_to_info[known_field_enum].get_num_elements_for_known_field_enum(num_ALT_alleles, ploidy);
}

unsigned VariantQueryProcessor::get_length_descriptor_for_known_field_enum(unsigned known_field_enum) const
{
  assert(known_field_enum < m_known_field_enum_to_info.size());
  return m_known_field_enum_to_info[known_field_enum].get_length_descriptor();
}

void VariantQueryProcessor::fill_field(std::unique_ptr<VariantFieldBase>& field_ptr,
    const CellConstAttrIterator& attr_iter,
    const unsigned num_ALT_alleles, const unsigned ploidy,
    unsigned schema_idx, uint64_t* num_deref_tile_iters
    ) const
{
  if(field_ptr.get() == nullptr)       //Allocate only if null
    field_ptr = std::move(m_field_factory.Create(schema_idx)); 
  unsigned known_field_enum = m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(schema_idx);
  //For known fields, check length descriptors - default FIXED
  unsigned length_descriptor = BCF_VL_FIXED;
  unsigned num_elements = 1u;
  if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(known_field_enum))
  {
    length_descriptor = get_length_descriptor_for_known_field_enum(known_field_enum);
    num_elements = get_num_elements_for_known_field_enum(known_field_enum, num_ALT_alleles, ploidy);
  }
  field_ptr->set_valid(true);  //mark as valid, since tile is actually accessed
  //This function might mark the field as invalid - some fields are  determined to be invalid only
  //after accessing the data and comparing to NULL_* values
  field_ptr->copy_data_from_tile(attr_iter, length_descriptor, num_elements);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
}

void VariantQueryProcessor::gt_fill_row(
    Variant& variant, int64_t row, int64_t column,
    const VariantQueryConfig& query_config,
    const Cell& cell, uint64_t* num_deref_tile_iters) const {
  #if VERBOSE>1
    std::cout << "[query_variants:gt_fill_row] Fill Row " << row << " column " << column << std::endl;
  #endif
  //Current row should be part of query
  assert(query_config.is_queried_array_row_idx(row));
  VariantCall& curr_call = variant.get_call(query_config.get_query_row_idx_for_array_row_idx(row));
  //Curr call will be initialized, one way or the other
  curr_call.mark_initialized(true);
  // First check if the row contains valid data, i.e., check whether the interval intersects with the current queried interval
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  auto* END_ptr = cell[query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX)].operator const int64_t*();
  auto END_v = *END_ptr; 
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  if(END_v < variant.get_column_begin()) {
    curr_call.mark_valid(false);  //The interval in this cell stops before the current variant's start column
    return;                       //Implies, this row has no valid data for this query range. Mark Call as invalid
  }
  curr_call.mark_valid(true);   //contains valid data for this query
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
        num_ALT_alleles, ploidy,
        query_config.get_schema_idx_for_query_idx(i), num_deref_tile_iters
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
        num_ALT_alleles, ploidy,
        query_config.get_schema_idx_for_query_idx(i), num_deref_tile_iters
        );     
  }
}

inline
unsigned int VariantQueryProcessor::gt_initialize_reverse_iter(
    const int ad,
    const VariantQueryConfig& query_config, const int64_t column,
    ArrayConstReverseCellIterator<int64_t>*& reverse_iter) const {
  assert(query_config.is_bookkeeping_done());
  //Num attributes in query
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  //Assign reverse iterator
  //FIXME: No longer binary search?
  vector<int64_t> query_range = { 0ll, static_cast<int64_t>(query_config.get_num_rows_in_array()),
    0ll, column };
  reverse_iter = get_storage_manager()->rbegin<int64_t>(ad, &(query_range[0]), query_config.get_query_attributes_schema_idxs());
  return num_queried_attributes - 1;
}


inline
unsigned int VariantQueryProcessor::gt_initialize_forward_iter(
    const int ad,
    const VariantQueryConfig& query_config, const int64_t column,
    ArrayConstCellIterator<int64_t>*& forward_iter) const {
  assert(query_config.is_bookkeeping_done());
  //Num attributes in query
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  //Assign forward iterator
  //FIXME: No longer binary search?
  vector<int64_t> query_range = { 0ll, static_cast<int64_t>(query_config.get_num_rows_in_array()),
    column, 10000000000ull };
  forward_iter = get_storage_manager()->begin<int64_t>(ad, &(query_range[0]), query_config.get_query_attributes_schema_idxs());
  return num_queried_attributes - 1;
}

void VariantQueryProcessor::clear()
{
  m_known_field_enum_to_info.clear();
  m_schema_idx_to_known_variant_field_enum_LUT.reset_luts();
  m_field_factory.clear();
}

bool VariantQueryProcessor::get_known_field_enum_for_name(const std::string& field_name, unsigned& known_field_enum)
{
  assert(VariantQueryProcessor::m_are_static_members_initialized);
  auto iter = VariantQueryProcessor::m_known_variant_field_name_to_enum.find(field_name);
  if(iter == VariantQueryProcessor::m_known_variant_field_name_to_enum.end())
    return false;
  known_field_enum = (*iter).second;
  return true;
}

string VariantQueryProcessor::get_known_field_name_for_enum(unsigned known_field_enum)
{
  assert(VariantQueryProcessor::m_are_static_members_initialized);
  assert(known_field_enum < GVCF_NUM_KNOWN_FIELDS);
  return VariantQueryProcessor::m_known_variant_field_names[known_field_enum];
}
