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
    std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldData<std::string, const AttributeTile<char>>>());
  //Set initialized flag
  VariantQueryProcessor::m_are_static_members_initialized = true;
}

void VariantQueryProcessor::initialize_known(const StorageManager::ArrayDescriptor* ad)
{
  //Initialize NULL bit idx for all known variant fields
  m_known_field_enum_to_info[GVCF_AC_IDX].m_NULL_bitidx =                        GVCF_AC_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_AN_IDX].m_NULL_bitidx =                        GVCF_AN_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_AF_IDX].m_NULL_bitidx =                        GVCF_AF_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_PL_IDX].m_NULL_bitidx =                        GVCF_PL_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_AD_IDX].m_NULL_bitidx =                        GVCF_AD_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_SB_4_IDX].m_NULL_bitidx =                      GVCF_SB_4_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_SB_3_IDX].m_NULL_bitidx =                      GVCF_SB_3_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_SB_2_IDX].m_NULL_bitidx =                      GVCF_SB_2_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_SB_1_IDX].m_NULL_bitidx =                      GVCF_SB_1_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_GQ_IDX].m_NULL_bitidx =                        GVCF_GQ_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_MIN_DP_IDX].m_NULL_bitidx =                    GVCF_MIN_DP_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_DP_FMT_IDX].m_NULL_bitidx =                    GVCF_DP_FMT_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_MQ0_IDX].m_NULL_bitidx =                       GVCF_MQ0_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_MQ_IDX].m_NULL_bitidx =                        GVCF_MQ_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_DP_IDX].m_NULL_bitidx =                        GVCF_DP_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_READPOSRANKSUM_IDX].m_NULL_bitidx =            GVCF_READPOSRANKSUM_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_MQRANKSUM_IDX].m_NULL_bitidx =                 GVCF_MQRANKSUM_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_CLIPPINGRANKSUM_IDX].m_NULL_bitidx =           GVCF_CLIPPINGRANKSUM_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_BASEQRANKSUM_IDX].m_NULL_bitidx =              GVCF_BASEQRANKSUM_NULL_BITIDX;
  m_known_field_enum_to_info[GVCF_QUAL_IDX].m_NULL_bitidx =                      GVCF_QUAL_NULL_BITIDX;
  //Initialize OFFSETS idx for known fields
  m_known_field_enum_to_info[GVCF_REF_IDX].m_OFFSETS_idx =  GVCF_REF_OFFSET_IDX;
  m_known_field_enum_to_info[GVCF_ALT_IDX].m_OFFSETS_idx =  GVCF_ALT_OFFSET_IDX;
  m_known_field_enum_to_info[GVCF_FILTER_IDX].m_OFFSETS_idx =  GVCF_FILTER_OFFSET_IDX;
  m_known_field_enum_to_info[GVCF_AD_IDX].m_OFFSETS_idx =  GVCF_AD_OFFSET_IDX;
  m_known_field_enum_to_info[GVCF_PL_IDX].m_OFFSETS_idx =  GVCF_PL_OFFSET_IDX;
  //Could change based on version as well
  m_num_elements_per_offset_cell = GVCF_NUM_KNOWN_OFFSET_ELEMENTS_PER_CELL;
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
      assert(m_known_field_enum_to_info[i].m_NULL_bitidx >= GVCF_PL_NULL_BITIDX);
      m_known_field_enum_to_info[i].m_NULL_bitidx -= GVCF_PL_NULL_BITIDX;
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
          m_known_field_enum_to_info[i].m_NULL_bitidx +=  diff;
      }
      diff = (GVCF_AF_NULL_BITIDX - GVCF_AC_NULL_BITIDX);
      //validate null bit idx for AF, AC, AN
      for(unsigned i=GVCF_AF_IDX;i<=GVCF_AC_IDX;++i)
        m_known_field_enum_to_info[i].m_NULL_bitidx =  diff - (i - GVCF_AF_IDX);
      //set schema version
      m_GT_schema_version = GT_SCHEMA_V1;
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
    case GVCF_OFFSETS_IDX:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      m_known_field_enum_to_info[idx].m_num_elements = m_num_elements_per_offset_cell;
      break;
    default:
      m_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      m_known_field_enum_to_info[idx].m_num_elements = 1u;
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

void VariantQueryProcessor::register_field_creators(const StorageManager::ArrayDescriptor* ad)
{
  const auto& schema = ad->array_schema();
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

VariantQueryProcessor::VariantQueryProcessor(const std::string& workspace, StorageManager& storage_manager,
    const StorageManager::ArrayDescriptor* ad)
: QueryProcessor(workspace, storage_manager)
{
  if(!VariantQueryProcessor::m_are_static_members_initialized)
    VariantQueryProcessor::initialize_static_members();
  clear();
  //Mapping from known_field enum
  m_known_field_enum_to_info.resize(GVCF_NUM_KNOWN_FIELDS);
  //Invalidate everything
  for(auto i=0u;i<m_known_field_enum_to_info.size();++i)
  {
    //Invalidate m_known_field_enum_to_info 
    invalidate_NULL_bitidx(i);
    //set all fields NOT to use offset
    m_known_field_enum_to_info[i].m_OFFSETS_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
  }
  //Initialize versioning information
  initialize_version(ad);
  //set length descriptors and creator objects for special attributes
  for(auto i=0u;i<m_known_field_enum_to_info.size();++i)
    initialize_length_descriptor(i);
  //Register creators in factory
  register_field_creators(ad);
}

void VariantQueryProcessor::handle_gvcf_ranges(VariantCallEndPQ& end_pq,
    const VariantQueryConfig& query_config, Variant& variant,
    std::ostream& output_stream,
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
      for(VariantCall& curr_call : variant.get_calls())
      {
        if(curr_call.is_valid())
        {
          assert(curr_call.get_column_begin() <= current_start_position);
          if(curr_call.get_column_begin() < current_start_position) 
          {
            auto* REF_ptr = get_known_field<VariantFieldString,true>
              (curr_call, query_config, GVCF_REF_IDX);
            REF_ptr->get() = "N";
          }
        }
      }
    VariantOperations::do_dummy_genotyping(variant, output_stream);
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
    std::ostream& output_stream) const
{
  assert(query_config.is_bookkeeping_done());
  
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  StorageManager::const_iterator* tile_its = new StorageManager::const_iterator[num_queried_attributes];
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
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  unsigned END_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX);
  //Priority queue of VariantCalls ordered by END positions
  VariantCallEndPQ end_pq;
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
  //Variant object
  Variant variant(&query_config);
  variant.resize_based_on_query();
  variant.set_column_interval(current_start_position, current_start_position);
  //Next co-ordinate to consider
  int64_t next_start_position = -1ll;
  uint64_t tile_idx = 0ull;
  uint64_t num_deref_tile_iters = 0ull; //profiling
  for(;tile_its[COORDS_query_idx] != tile_it_end;advance_tile_its(num_queried_attributes-1, tile_its))  //why -1, advance uses i<=N for loop
  {
    // Initialize cell iterators for the coordinates
    Tile::const_iterator cell_it = (*tile_its[COORDS_query_idx]).begin();
    Tile::const_iterator cell_it_end = (*tile_its[COORDS_query_idx]).end();
    for(;cell_it != cell_it_end;++cell_it) {
      std::vector<int64_t> next_coord = *cell_it;
      if(next_coord[1] != current_start_position) //have found cell with next gVCF position, handle accumulated values
      {
        next_start_position = next_coord[1];
        assert(next_coord[1] > current_start_position);
        handle_gvcf_ranges(end_pq, query_config, variant, output_stream, current_start_position,
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
  }
  //handle last interval
  handle_gvcf_ranges(end_pq, query_config, variant, output_stream, current_start_position, 0, true);
  delete[] tile_its;
}

void VariantQueryProcessor::do_query_bookkeeping(const StorageManager::ArrayDescriptor* array_descriptor,
    VariantQueryConfig& query_config) const
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
      //Does the field require OFFSETS?
      if(uses_OFFSETS_field_for_known_field_enum(known_variant_field_enum))
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
      //Does the length of the field depend on the number of alleles
      if(is_length_allele_dependent(known_variant_field_enum))
      {
        assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(ALT_schema_idx));
        query_config.add_attribute_to_query("ALT", ALT_schema_idx);
        //ALT uses OFFSETS
        assert(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(OFFSETS_schema_idx));
        query_config.add_attribute_to_query("OFFSETS", OFFSETS_schema_idx);
      }
    }
  }
  //Re-order query fields so that special fields are first
  query_config.reorder_query_fields();
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
  //Set number of rows in the array
  const ArraySchema& array_schema = array_descriptor->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t row_num = dim_domains[0].second - dim_domains[0].first + 1;
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

void VariantQueryProcessor::gt_get_column(
    const StorageManager::ArrayDescriptor* ad,
    const VariantQueryConfig& query_config, unsigned column_interval_idx,
    Variant& variant, GTProfileStats* stats) const {

  //New interval starts
  variant.reset_for_new_interval();

  assert(query_config.get_num_column_intervals() > 0 && column_interval_idx < query_config.get_num_column_intervals());
  assert(query_config.is_bookkeeping_done() && "Before calling gt_get_column(), do_query_bookkeeping() function must be called");
  assert(query_config.get_first_normal_field_query_idx() >= 2u); //COORDINATES, END are required by default
  
  variant.set_column_interval(query_config.get_column_interval(column_interval_idx).first,
          query_config.get_column_interval(column_interval_idx).second);
  //TODO: Still single position query
  uint64_t col = query_config.get_column_interval(column_interval_idx).first;
  // Initialize reverse tile iterators for 
  // queried fields 
  // The reverse tile iterator will start with the tiles
  // of the various attributes that have the largest
  // id that either intersect with col, or precede col.
  StorageManager::const_reverse_iterator* tile_its;
  StorageManager::const_reverse_iterator tile_it_end;  
  unsigned int gt_attribute_num = 
      gt_initialize_tile_its(ad, query_config, column_interval_idx, tile_its, tile_it_end);
  //Get query idx for COORDS
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  unsigned COORDS_query_idx = query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX);
  // Create cell iterators
  Tile::const_reverse_iterator cell_it, cell_it_end;
  uint64_t num_deref_tile_iters = 0;
#ifdef DO_PROFILING
  uint64_t num_cells_touched = 0;
  uint64_t num_tiles_touched = 0;
  bool first_sample = true;
#endif
  // Indicates how many rows have been filled.
  uint64_t filled_rows = 0;
  // Fill the genotyping column
  while(tile_its[COORDS_query_idx] != tile_it_end && filled_rows < query_config.get_num_rows_to_query()) {
    // Initialize cell iterators for the coordinates
    cell_it = (*tile_its[COORDS_query_idx]).rbegin();
    cell_it_end = (*tile_its[COORDS_query_idx]).rend();
#ifdef DO_PROFILING
    num_deref_tile_iters += 2;	//why 2, cell_it and cell_it_end
    ++num_tiles_touched;
#endif
    while(cell_it != cell_it_end && filled_rows < query_config.get_num_rows_to_query()) {
      std::vector<int64_t> next_coord = *cell_it;
#ifdef DO_PROFILING
      ++num_cells_touched;
#endif
      // If next cell is not on the right of col, and 
      // The rowIdx is being queried and
      // The row/call is uninitialized (uninvestigated) in the Variant
      if(next_coord[1] <= col && query_config.is_queried_array_row_idx(next_coord[0]) &&
          !(variant.get_call(query_config.get_query_row_idx_for_array_row_idx(next_coord[0])).is_initialized())
        ) {
        gt_fill_row<StorageManager::const_reverse_iterator>(variant, next_coord[0], next_coord[1], cell_it.pos(), query_config,
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
}

unsigned VariantQueryProcessor::get_num_elements_for_known_field_enum(unsigned known_field_enum,
    unsigned num_ALT_alleles) const
{
  assert(known_field_enum < m_known_field_enum_to_info.size());
  unsigned length = 0u;
  unsigned num_alleles = num_ALT_alleles + 1;
  switch(m_known_field_enum_to_info[known_field_enum].m_length_descriptor)
  {
    case BCF_VL_FIXED:
      length = m_known_field_enum_to_info[known_field_enum].m_num_elements;
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
    default:
      cerr << "Unknown length descriptor "<<m_known_field_enum_to_info[known_field_enum].m_length_descriptor<<" - ignoring\n";
      break;
  }
  return length;
}

void VariantQueryProcessor::fill_field(std::unique_ptr<VariantFieldBase>& field_ptr,
    const Tile& tile, int64_t pos,
    const vector<int64_t>* OFFSETS_values, const int NULL_bitmap, const unsigned num_ALT_alleles,
    unsigned schema_idx, uint64_t* num_deref_tile_iters
    ) const
{
  if(field_ptr.get() == nullptr)       //Allocate only if null
    field_ptr = std::move(m_field_factory.Create(schema_idx)); 
  unsigned known_field_enum = m_schema_idx_to_known_variant_field_enum_LUT.get_known_field_enum_for_schema_idx(schema_idx);
  uint64_t num_elements = 0ull;
  //Default value of offset == offset of cell in co-ordinates tile
  //Except for OFFSETS field
  uint64_t field_offset = (known_field_enum == GVCF_OFFSETS_IDX) ? m_num_elements_per_offset_cell*pos : pos;
  //For known fields, check NULL, OFFSETS, length descriptors
  if(m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value(known_field_enum))
  {
    //NULL bit
    if(is_NULL_bitidx_defined_for_known_field_enum(known_field_enum))
    {
      auto NULL_bitidx = get_NULL_bitidx_for_known_field_enum(known_field_enum);
      if((NULL_bitmap >> NULL_bitidx) & 1) //is null, nothing to do
        return;
    }
    num_elements = get_num_elements_for_known_field_enum(known_field_enum, num_ALT_alleles);
    //Uses OFFSETS
    if(uses_OFFSETS_field_for_known_field_enum(known_field_enum))
    {
      assert(OFFSETS_values);
      auto OFFSETS_idx = get_OFFSETS_idx_for_known_field_enum(known_field_enum);
      assert(OFFSETS_idx < OFFSETS_values->size());
      field_offset = OFFSETS_values->operator[](OFFSETS_idx);
    }
  }
  field_ptr->copy_data_from_tile(tile, field_offset, num_elements);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
}

template<class ITER>
void VariantQueryProcessor::gt_fill_row(
    Variant& variant, int64_t row, int64_t column, int64_t pos,
    const VariantQueryConfig& query_config,
    const ITER* tile_its, uint64_t* num_deref_tile_iters) const {
  //Current row should be part of query
  assert(query_config.is_queried_array_row_idx(row));
  VariantCall& curr_call = variant.get_call(query_config.get_query_row_idx_for_array_row_idx(row));
  //Curr call will be initialized, one way or the other
  curr_call.mark_initialized(true);
  // First check if the row contains valid data, i.e., check whether the interval intersects with the current queried interval
  assert(query_config.is_defined_query_idx_for_known_field_enum(GVCF_END_IDX));
  int64_t END_v = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX)]).cell(pos);
#ifdef DO_PROFILING
  ++(*num_deref_tile_iters);
#endif
  if(END_v < variant.get_column_begin()) {
    curr_call.mark_valid(false);  //The interval in this cell stops before the current variant's start column
    return;                     //Implies, this row has no valid data for this query range. Mark Call as invalid
  }
  curr_call.mark_valid(true);   //contains valid data for this query
  //Set begin,end of the Call - NOTE: need not be same as Variant's begin,end
  curr_call.set_column_interval(column, END_v);
  //Coordinates and END should be the first 2 queried attributes - see reorder_query_fields()
  assert(query_config.get_query_idx_for_known_field_enum(GVCF_COORDINATES_IDX) < 2u);
  assert(query_config.get_query_idx_for_known_field_enum(GVCF_END_IDX) < 2u);
  assert(query_config.get_first_normal_field_query_idx() >= 2u);
  assert(query_config.get_first_normal_field_query_idx() != UNDEFINED_ATTRIBUTE_IDX_VALUE);
  assert(query_config.get_after_OFFSETS_field_query_idx() >= 2u);
  assert(query_config.get_after_OFFSETS_field_query_idx() != UNDEFINED_ATTRIBUTE_IDX_VALUE);
  //Variables to store special fields
  //OFFSETS values - initialized later
  const vector<int64_t>* OFFSETS_values = 0;
  // NULL bitmap, if needed
  int NULL_bitmap = 0;
  //Num alternate alleles
  unsigned num_ALT_alleles = 0u;
  //Load special fields up to and including OFFSETS
  //NOTE, none of the special fields should use NULL bitmap or OFFSETS field
  for(auto i=2u;i<query_config.get_after_OFFSETS_field_query_idx();++i)
  {
    //Read from Tile
    fill_field(curr_call.get_field(i),(*tile_its[i]), pos,
        OFFSETS_values, NULL_bitmap, num_ALT_alleles,
        query_config.get_schema_idx_for_query_idx(i), num_deref_tile_iters
        );
  }
  //Initialize NULL_bitmap, if needed
  const auto* NULL_field_ptr =
    get_known_field_if_queried<VariantFieldPrimitiveVectorData<int>, true>(curr_call, query_config, GVCF_NULL_IDX); 
  if(NULL_field_ptr)
    NULL_bitmap = NULL_field_ptr->get()[0];        //NULL_field data is vector<int>
  const auto* OFFSETS_field_ptr =
    get_known_field_if_queried<VariantFieldPrimitiveVectorData<int64_t>, true>(curr_call, query_config, GVCF_OFFSETS_IDX); 
  if(OFFSETS_field_ptr)
    OFFSETS_values = &(OFFSETS_field_ptr->get());
  for(auto i=query_config.get_after_OFFSETS_field_query_idx();
      i<query_config.get_first_normal_field_query_idx();++i)
  {
    //Read from Tile
    fill_field(curr_call.get_field(i),(*tile_its[i]), pos,
        OFFSETS_values, NULL_bitmap, num_ALT_alleles,
        query_config.get_schema_idx_for_query_idx(i), num_deref_tile_iters
        );
  }
  //Initialize ALT field, if needed
  const auto* ALT_field_ptr = get_known_field_if_queried<VariantFieldALTData, true>(curr_call, query_config, GVCF_ALT_IDX); 
  num_ALT_alleles = ALT_field_ptr->get().size();   //ALT field data is vector<string>

  //Go over all normal query fields and fetch data
  for(auto i=query_config.get_first_normal_field_query_idx();i<query_config.get_num_queried_attributes();++i)
  {
    //Read from Tile
    fill_field(curr_call.get_field(i),(*tile_its[i]), pos,
        OFFSETS_values, NULL_bitmap, num_ALT_alleles,
        query_config.get_schema_idx_for_query_idx(i), num_deref_tile_iters
        );     
  }
  VariantFieldString* REF_field_ptr =
    get_known_field_if_queried<VariantFieldString, true>(curr_call, query_config, GVCF_REF_IDX);
  if(REF_field_ptr)
  {
    if(column < variant.get_column_begin())
      REF_field_ptr->get() = "N";
  }
}

inline
unsigned int VariantQueryProcessor::gt_initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    const VariantQueryConfig& query_config, const unsigned column_interval_idx,
    StorageManager::const_reverse_iterator*& tile_its, 
    StorageManager::const_reverse_iterator& tile_it_end) const {
  //Num attributes in query
  unsigned num_queried_attributes = query_config.get_num_queried_attributes();
  // Create reverse iterators
  tile_its = new StorageManager::const_reverse_iterator[num_queried_attributes];
  //TODO: still assumes single position query
  // Find the rank of the tile the left sweep starts from.
  auto start_rank = get_storage_manager().get_left_sweep_start_rank(ad, query_config.get_column_interval(column_interval_idx).first);
  for(auto i=0u;i<query_config.get_num_queried_attributes();++i)
  {
    assert(query_config.is_schema_idx_defined_for_query_idx(i));
    auto schema_idx = query_config.get_schema_idx_for_query_idx(i);
    tile_its[i] = get_storage_manager().rbegin(ad, schema_idx, start_rank);
  }
  tile_it_end = get_storage_manager().rend(ad, m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(GVCF_COORDINATES_IDX));
  return num_queried_attributes - 1;
}

void VariantQueryProcessor::clear()
{
  m_known_field_enum_to_info.clear();
  m_schema_idx_to_known_variant_field_enum_LUT.reset_luts();
  m_field_factory.clear();
}
