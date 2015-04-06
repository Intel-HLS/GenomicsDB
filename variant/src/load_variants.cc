#include "load_variants.h"

VariantLoader::VariantLoader(const std::string& workspace, StorageManager& storage_manager)
: Loader(workspace, storage_manager)
{ 
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
  GVCF_AF_IDX = GVCF_COORDINATES_IDX++;
  GVCF_AN_IDX = GVCF_COORDINATES_IDX++;
  GVCF_AC_IDX = GVCF_COORDINATES_IDX++;
  GVCF_NULL_IDX = GVCF_COORDINATES_IDX++;
  GVCF_OFFSETS_IDX = GVCF_COORDINATES_IDX++;
}

void VariantLoader::load_CSV_gVCF(const std::string& filename, const char* array_name, const uint64_t max_sample_idx, const bool is_input_sorted, const std::string temp_space) const {
  // Create gVCF array schema
  ArraySchema* array_schema = create_gVCF_array_schema(array_name, max_sample_idx);

  // Open array in CREATE mode
  StorageManager::ArrayDescriptor* ad = 
      get_storage_manager().open_array(*array_schema);

  // Prepare filenames
  std::string to_be_sorted_filename = filename;
  if(to_be_sorted_filename[0] == '~') 
    to_be_sorted_filename = std::string(getenv("HOME")) +
        to_be_sorted_filename.substr(1, get_workspace().size()-1);
  assert(check_on_load(to_be_sorted_filename));
  std::string sorted_filename = get_workspace() + "/sorted_" +
                                array_schema->array_name() + ".csv";

  // Sort CSV
  if(is_input_sorted)
    sorted_filename = to_be_sorted_filename;
  else
    sort_csv_file(to_be_sorted_filename, sorted_filename, *array_schema, temp_space);

  // Make tiles
  try {
    make_tiles_irregular_CSV_gVCF(sorted_filename, ad, *array_schema);
  } catch(LoaderException& le) {
    if(!is_input_sorted)
      remove(sorted_filename.c_str());
    get_storage_manager().delete_array(array_schema->array_name());
    throw LoaderException("[Loader] Cannot load CSV file '" + filename + 
                          "'.\n " + le.what());
  } 

  // Clean up and close array 
  if(!is_input_sorted)
    remove(sorted_filename.c_str());
  get_storage_manager().close_array(ad);
  delete array_schema;
}

void VariantLoader::append_cell_gVCF(const ArraySchema& array_schema, 
                              CSVLine* csv_line, Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = array_schema.attribute_num();

  // Bitmap that indicates the NULL attribute values
  int NULL_bitmap = 0;  
 
  // Variables for retrieving the CSV values
  int int_v;
  int64_t int64_t_v;
  float float_v;
  std::string string_v;

  // Offsets
  int64_t REF_offset = tiles[GVCF_REF_IDX]->cell_num(); 
  int64_t ALT_offset = tiles[GVCF_ALT_IDX]->cell_num(); 
  int64_t FILTER_ID_offset = tiles[GVCF_FILTER_IDX]->cell_num(); 
  int64_t AD_offset = tiles[GVCF_AD_IDX]->cell_num(); 
  int64_t PL_offset = tiles[GVCF_PL_IDX]->cell_num(); 
  
  // ALT_num and FILTER_num variables
  int ALT_num;
  int FILTER_num;

  // Append coordinates first
  if(!(*csv_line >> *tiles[attribute_num]))
    throw LoaderException("Cannot read coordinates from CSV file.");

  // END -- Attribute 0
  if(!(*csv_line >> *tiles[GVCF_END_IDX]))
    throw LoaderException("Cannot read END value from CSV file.");

  // REF -- Attribute GVCF_REF_IDX
  if(!(*csv_line >> string_v))
    throw LoaderException("Cannot read REF value from CSV file.");
  for(int i=0; i<string_v.size(); i++)
    *tiles[GVCF_REF_IDX] << string_v[i];
  *tiles[GVCF_REF_IDX] << '\0';

  // ALT -- Attribute GVCF_ALT_IDX
  if(!(*csv_line >> ALT_num))
    throw LoaderException("Cannot read ALT_num value from CSV file.");
  for(int j=0; j<ALT_num; j++) {
    if(!(*csv_line >> string_v))
      throw LoaderException("Cannot read ALT value from CSV file.");
    for(int i=0; i<string_v.size(); i++)
      *tiles[GVCF_ALT_IDX] << string_v[i];
    if(j < ALT_num-1) // Do not store '\0' after last allele '&' for <NON_REF>
      *tiles[GVCF_ALT_IDX] << '\0';
    else
      assert(string_v == "&");
  }
  
  // QUAL -- Attribute GVCF_QUAL_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read QUAL value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_QUAL_IDX] << 0;
  } else {  
    *tiles[GVCF_QUAL_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // FILTER_ID -- Attribute GVCF_FILTER_IDX
  if(!(*csv_line >> FILTER_num))
    throw LoaderException("Cannot read FILTER_num value from CSV file.");
  *tiles[GVCF_FILTER_IDX] << FILTER_num;  
  for(int j=0; j<FILTER_num; j++) {
    if(!(*csv_line >> *tiles[GVCF_FILTER_IDX]))
      throw LoaderException("Cannot read FILTER_ID value from CSV file.");
  }

  // BaseQRankSum -- Attribute GVCF_BASEQRANKSUM_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read BaseQRankSum value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_BASEQRANKSUM_IDX] << 0;
  } else {
    *tiles[GVCF_BASEQRANKSUM_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // ClippingRankSum -- Attribute GVCF_CLIPPINGRANKSUM_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read ClippingSum value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_CLIPPINGRANKSUM_IDX] << 0;
  } else {
    *tiles[GVCF_CLIPPINGRANKSUM_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // MQRankSum -- Attribute GVCF_MQRANKSUM_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read MQRankSum value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_MQRANKSUM_IDX] << 0;
  } else {
    *tiles[GVCF_MQRANKSUM_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // ReadPosRankSum -- Attribute GVCF_READPOSRANKSUM_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read ReadPosRankSum value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_READPOSRANKSUM_IDX] << 0;
  } else {
    *tiles[GVCF_READPOSRANKSUM_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // DP -- Attribute GVCF_DP_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read DP value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_DP_IDX] << 0;
  } else {
    *tiles[GVCF_DP_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // MQ -- Attribute GVCF_MQ_IDX
  if(!(*csv_line >> float_v))
    throw LoaderException("Cannot read MQ value from CSV file.");
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_MQ_IDX] << 0;
  } else {
    *tiles[GVCF_MQ_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // MQ0 -- Attribute GVCF_MQ0_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read MQ0 value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_MQ0_IDX] << 0;
  } else {
    *tiles[GVCF_MQ0_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // DP_FMT -- Attribute GVCF_DP_FMT_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read DP_FMT value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_DP_FMT_IDX] << 0;
  } else {
    *tiles[GVCF_DP_FMT_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // MIN_DP -- Attribute GVCF_MIN_DP_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read MIN_DP value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_MIN_DP_IDX] << 0;
  } else {
    *tiles[GVCF_MIN_DP_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // GQ -- Attribute GVCF_GQ_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read GQ value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_GQ_IDX] << 0;
  } else {
    *tiles[GVCF_GQ_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // SB_1 -- Attribute GVCF_SB_1_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read SB_1 value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_SB_1_IDX] << 0;
  } else {
    *tiles[GVCF_SB_1_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // SB_2 -- Attribute GVCF_SB_2_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read SB_2 value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_SB_2_IDX] << 0;
  } else {
    *tiles[GVCF_SB_2_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // SB_3 -- Attribute GVCF_SB_3_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read SB_3 value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_SB_3_IDX] << 0;
  } else {
    *tiles[GVCF_SB_3_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // SB_4 -- Attribute GVCF_SB_4_IDX
  if(!(*csv_line >> int_v))
    throw LoaderException("Cannot read SB_4 value from CSV file.");
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_SB_4_IDX] << 0;
  } else {
    *tiles[GVCF_SB_4_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1;

  // AD -- Attribute GVCF_AD_IDX
  for(int i=0; i<ALT_num+1; i++) {
    if(!(*csv_line >> int_v))
      throw LoaderException("Cannot read AD value from CSV file.");
    if(int_v == CSV_NULL_INT)
      *tiles[GVCF_AD_IDX] << 0;
    else 
      *tiles[GVCF_AD_IDX] << int_v;
  }
  if(int_v == CSV_NULL_INT) // We assume that if one AD value is NULL, all are
    NULL_bitmap += 1;
  NULL_bitmap = NULL_bitmap << 1;
  
  // PL -- Attribute GVCF_PL_IDX
  for(int i=0; i<(ALT_num+1)*(ALT_num+2)/2; i++) {
    if(!(*csv_line >> int_v))
      throw LoaderException("Cannot read PL value from CSV file.");
    if(int_v == CSV_NULL_INT)
      *tiles[GVCF_PL_IDX] << 0;
    else 
      *tiles[GVCF_PL_IDX] << int_v;
  }
  if(int_v == CSV_NULL_INT) // We assume that if one PL value is NULL, all are
    NULL_bitmap += 1;
  NULL_bitmap = NULL_bitmap << 1;

  // AF -- Attribute GVCF_AF_IDX
  if(!(*csv_line >> float_v))
    // throw LoaderException("Cannot read AF value from CSV file.");
    // We are not throwing an exception here because we want to be backward 
    // compatible with VCF Files
    // Assign float_v to CSV_NULL_FLOAT and fall through rest of the logic
    float_v = CSV_NULL_FLOAT;
  if(float_v == CSV_NULL_FLOAT) {
    NULL_bitmap += 1;
    *tiles[GVCF_AF_IDX] << 0;
  } else {  
    *tiles[GVCF_AF_IDX] << float_v;
  }
  NULL_bitmap = NULL_bitmap << 1; 

  // AN -- Attribute GVCF_AN_IDX
  if(!(*csv_line >> int_v))
    // throw LoaderException("Cannot read AN value from CSV file.");
    // We are not throwing an exception here because we want to be backward 
    // compatible with VCF Files
    // Assign int_v to CSV_NULL_INT and fall through rest of the logic
    int_v = CSV_NULL_INT;
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_AN_IDX] << 0;
  } else {
    *tiles[GVCF_AN_IDX] << int_v;
  }
  NULL_bitmap = NULL_bitmap << 1;

  // AC -- Attribute GVCF_AC_IDX
  if(!(*csv_line >> int_v))
    // throw LoaderException("Cannot read AC value from CSV file.");
    // We are not throwing an exception here because we want to be backward 
    // compatible with VCF Files
    // Assign int_v to CSV_NULL_INT and fall through rest of the logic
    int_v = CSV_NULL_INT;
  if(int_v == CSV_NULL_INT) {
    NULL_bitmap += 1;
    *tiles[GVCF_AC_IDX] << 0;
  } else {
    *tiles[GVCF_AC_IDX] << int_v;
  }
  
  // NULL -- Attribute GVCF_NULL_IDX
  *tiles[GVCF_NULL_IDX] << NULL_bitmap;

  // OFFSETS -- Attribute GVCF_OFFSETS_IDX
  *tiles[GVCF_OFFSETS_IDX] << REF_offset;
  *tiles[GVCF_OFFSETS_IDX] << ALT_offset;
  *tiles[GVCF_OFFSETS_IDX] << FILTER_ID_offset;
  *tiles[GVCF_OFFSETS_IDX] << AD_offset;
  *tiles[GVCF_OFFSETS_IDX] << PL_offset;
}

ArraySchema* VariantLoader::create_gVCF_array_schema(const char* array_name, const uint64_t max_sample_idx) const {
  // Set dimension names
  std::vector<std::string> dim_names;
  dim_names.push_back("SampleID"); 
  dim_names.push_back("POS"); 
  
  // Set dimension domains
  std::vector<std::pair<double,double> > dim_domains;
  dim_domains.push_back(std::pair<double,double>(0, max_sample_idx));
  dim_domains.push_back(std::pair<double,double>(0, 6000000000));

  // Set attribute names
  std::vector<std::string> attribute_names;
  attribute_names.push_back("END"); 
  attribute_names.push_back("REF"); 
  attribute_names.push_back("ALT"); 
  attribute_names.push_back("QUAL"); 
  attribute_names.push_back("FILTER_ID"); 
  attribute_names.push_back("BaseQRankSum"); 
  attribute_names.push_back("ClippingRankSum"); 
  attribute_names.push_back("MQRankSum"); 
  attribute_names.push_back("ReadPosRankSum"); 
  attribute_names.push_back("DP"); 
  attribute_names.push_back("MQ"); 
  attribute_names.push_back("MQ0"); 
  attribute_names.push_back("DP_FMT"); 
  attribute_names.push_back("MIN_DP"); 
  attribute_names.push_back("GQ"); 
  attribute_names.push_back("SB_1"); 
  attribute_names.push_back("SB_2"); 
  attribute_names.push_back("SB_3"); 
  attribute_names.push_back("SB_4"); 
  attribute_names.push_back("AD"); 
  attribute_names.push_back("PL"); 
  attribute_names.push_back("AF"); 
  attribute_names.push_back("AN"); 
  attribute_names.push_back("AC"); 
  attribute_names.push_back("NULL"); 
  attribute_names.push_back("OFFSETS"); 
 
  // Set attribute types. 
  // The first types are for the attributes, whereas the very last type is 
  // for all the dimensions collectively. 
  std::vector<const std::type_info*> types;
  types.push_back(&typeid(int64_t));  // END
  types.push_back(&typeid(char));     // REF 
  types.push_back(&typeid(char));     // ALT
  types.push_back(&typeid(float));    // QUAL
  types.push_back(&typeid(int));      // FILTER_ID
  types.push_back(&typeid(float));    // BaseQRankSum
  types.push_back(&typeid(float));    // ClippingRankSum
  types.push_back(&typeid(float));    // MQRankSum
  types.push_back(&typeid(float));    // ReadPosRankSum
  types.push_back(&typeid(int));      // DP
  types.push_back(&typeid(float));    // MQ
  types.push_back(&typeid(int));      // MQ0
  types.push_back(&typeid(int));      // DP_FMT
  types.push_back(&typeid(int));      // MIN_DP
  types.push_back(&typeid(int));      // GQ
  types.push_back(&typeid(int));      // SB_1
  types.push_back(&typeid(int));      // SB_2
  types.push_back(&typeid(int));      // SB_3
  types.push_back(&typeid(int));      // SB_4
  types.push_back(&typeid(int));      // AD
  types.push_back(&typeid(int));      // PL
  types.push_back(&typeid(float));    // AF
  types.push_back(&typeid(int));      // AN
  types.push_back(&typeid(int));      // AC
  types.push_back(&typeid(int));      // NULL
  types.push_back(&typeid(int64_t));  // OFFSETS 
  types.push_back(&typeid(int64_t));  // Coordinates (SampleID, POS)

  // Set order and capacity
  ArraySchema::Order order = ArraySchema::COLUMN_MAJOR;
  uint64_t capacity = 1000;

  // Return schema with regular tiles
  return new ArraySchema(array_name, attribute_names, dim_names, dim_domains, 
                         types, order, capacity);
}

void VariantLoader::make_tiles_irregular_CSV_gVCF(
    const std::string& filename,
    const StorageManager::ArrayDescriptor* ad, 
    const ArraySchema& array_schema) const {
  // For easy reference
  ArraySchema::Order order = array_schema.order();
  uint64_t capacity = array_schema.capacity();
 
  // Initialization 
  CSVFile csv_file(filename, CSVFile::READ);
  CSVLine csv_line;
  Tile** tiles = new Tile*[array_schema.attribute_num() + 1]; 
  uint64_t tile_id = 0;
  uint64_t cell_num = 0;

  new_tiles(array_schema, tile_id, tiles);

  while(csv_file >> csv_line) {
    if(cell_num == capacity) {
      store_tiles(ad, tiles);
      new_tiles(array_schema, ++tile_id, tiles);
      cell_num = 0;
    }
   
    try {
      // Every CSV line is a logical cell
      append_cell_gVCF(array_schema, &csv_line, tiles);     
    } catch(LoaderException& le) {
      delete [] tiles;
      throw LoaderException("[Make tiles] " + le.what()); 
    }
    cell_num++;
  }

  // Store the lastly created tiles
  store_tiles(ad, tiles);

  delete [] tiles; 
}


