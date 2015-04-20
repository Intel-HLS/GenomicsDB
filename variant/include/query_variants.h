#ifndef QUERY_VARIANTS_H
#define QUERY_VARIANTS_H

#include "query_processor.h"
#include "gt_common.h"
#include "variant.h"
#include "variant_query_config.h"
#include "variant_query_field_data.h"

enum GTSchemaVersionEnum
{
    GT_SCHEMA_V0=0,
    GT_SCHEMA_V1
};

//Bit positions in the NULL bitmap of all known field enums
//Note: this MAY not be the real positions in a given array schema, the real positions
//are initialized in the query_variants.cc::initialize_version() functions. This enum 
//simply lists the bitidx in the latest version of the Loader (or the latest version of
//the variant schema)
//Should be in the reverse order of fields as used in load_variants.cc
enum KnownFieldsNULLBitidxEnum
{
  GVCF_AC_NULL_BITIDX=0,
  GVCF_AN_NULL_BITIDX,
  GVCF_AF_NULL_BITIDX,
  GVCF_PL_NULL_BITIDX,
  GVCF_AD_NULL_BITIDX,
  GVCF_SB_4_NULL_BITIDX,
  GVCF_SB_3_NULL_BITIDX,
  GVCF_SB_2_NULL_BITIDX,
  GVCF_SB_1_NULL_BITIDX,
  GVCF_GQ_NULL_BITIDX,
  GVCF_MIN_DP_NULL_BITIDX,
  GVCF_DP_FMT_NULL_BITIDX,
  GVCF_MQ0_NULL_BITIDX,
  GVCF_MQ_NULL_BITIDX,
  GVCF_DP_NULL_BITIDX,
  GVCF_READPOSRANKSUM_NULL_BITIDX,
  GVCF_MQRANKSUM_NULL_BITIDX,
  GVCF_CLIPPINGRANKSUM_NULL_BITIDX,
  GVCF_BASEQRANKSUM_NULL_BITIDX,
  GVCF_QUAL_NULL_BITIDX,
  GVCF_NUM_NULL_BITS_USED
};
/* Structure to store profiling information */
class GTProfileStats {
  public:
    GTProfileStats()
    {
      m_sum_num_cells_touched = 0;
      m_sum_num_deref_tile_iters = 0;
      m_sum_num_tiles_touched = 0;
      m_sum_num_cells_last_iter = 0;
      m_sum_num_cells_first_sample = 0;
      m_sum_sq_num_cells_touched = 0;
      m_sum_sq_num_deref_tile_iters = 0;      
      m_sum_sq_num_tiles_touched = 0;
      m_sum_sq_num_cells_last_iter = 0;
      m_sum_sq_num_cells_first_sample = 0;
      m_num_samples = 0;
    }
    uint64_t m_sum_num_cells_touched;
    uint64_t m_sum_num_deref_tile_iters;
    uint64_t m_sum_num_tiles_touched;
    uint64_t m_sum_num_cells_last_iter;
    uint64_t m_sum_num_cells_first_sample;
    uint64_t m_sum_sq_num_cells_touched;
    uint64_t m_sum_sq_num_deref_tile_iters;
    uint64_t m_sum_sq_num_tiles_touched;
    uint64_t m_sum_sq_num_cells_last_iter;
    uint64_t m_sum_sq_num_cells_first_sample;
    uint64_t m_num_samples;
};

/*
 * Stores the tile iterator value for a given tile.
 * Used by scan_and_operate to track multiple tile iterators when variant intervals
 * span large ranges
 */
class GTTileIteratorsTracker
{
  public:
    GTTileIteratorsTracker(unsigned num_attributes)
    {
      m_iter_vector.resize(num_attributes);
      m_reference_counter = 0ull;
    }
    uint64_t m_reference_counter;
    std::vector<StorageManager::const_iterator> m_iter_vector;
};

/*
 * Child class of QueryProcessor customized to handle variants
 */
class VariantQueryProcessor : public QueryProcessor {

  public:
    /** 
     * Simple constructor. The workspace is where the query processor will create 
     * its data. The storage manager is the module the query processor interefaces 
     * with.
     */
    VariantQueryProcessor(const std::string& workspace, StorageManager& storage_manager,
        const StorageManager::ArrayDescriptor* ad);
    void clear();
    /**
     * When querying, setup bookkeeping structures first 
     */
    void do_query_bookkeeping(const StorageManager::ArrayDescriptor* array_descriptor,
        VariantQueryConfig& query_config);
    /** Returns the genotyping info for column col from the input array. */
    GTColumn* gt_get_column(
        const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config, GTProfileStats* stats=0) const;
    void scan_and_operate(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config, std::ostream& output_stream) const;
    void iterate_over_all_tiles(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config) const;
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum
     * returns true if the field requires NULL bitidx field to be used 
     */
    inline bool is_NULL_bitidx_defined_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size());
      return (m_known_field_enum_to_info[enumIdx].m_NULL_bitidx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
    }
    inline unsigned get_NULL_bitidx_for_known_field_enum(unsigned enumIdx) const
    {
        assert(enumIdx < m_known_field_enum_to_info.size());
        return (m_known_field_enum_to_info[enumIdx].m_NULL_bitidx);
    }
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum
     * returns true if the field requires the OFFSETS field 
     */
    inline bool uses_OFFSETS_field(unsigned enumIdx)
    {
      assert(enumIdx < m_known_field_enum_to_info.size()); 
      return (m_known_field_enum_to_info[enumIdx].m_OFFSETS_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
    }
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum
     * returns the schema idx for the given array 
     */
    inline unsigned get_schema_idx_for_known_field_enum(unsigned enumIdx)
    {
      assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
      return m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(enumIdx);
    }
  private:
    /*initialize all known info about variants*/
    void initialize_known(const StorageManager::ArrayDescriptor* ad);
    /*Initialize schema version v1 info*/
    void initialize_v0(const StorageManager::ArrayDescriptor* ad);
    /*Check and initialize schema version v2 info*/
    void initialize_v1(const StorageManager::ArrayDescriptor* ad);
    /*Initialize versioning information based on schema*/
    void initialize_version(const StorageManager::ArrayDescriptor* ad);
    /**
     * Function invalidating enum to NULL bitidx
     */
    void invalidate_NULL_bitidx(unsigned enumIdx)
    {
      assert(enumIdx < m_known_field_enum_to_info.size());
      m_known_field_enum_to_info[enumIdx].m_NULL_bitidx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    /**
     * Initialized field length info
     */
    void initialize_length_descriptor_and_fetch(unsigned idx);
    /** Called by scan_and_operate to handle all ranges for given set of cells */
    void handle_gvcf_ranges(VariantIntervalPQ& end_pq, std::vector<PQStruct>& PQ_end_vec,
            const VariantQueryConfig& queryConfig, GTColumn* gt_column,
        std::unordered_map<uint64_t, GTTileIteratorsTracker>& tile_idx_2_iters, std::ostream& output_stream,
        int64_t current_start_position, int64_t next_start_position, bool is_last_call) const;
    /** Fills a row of the input genotyping column with the proper info. */
    template<class ITER>
    void gt_fill_row(
        GTColumn* gt_column, int64_t row, int64_t column, int64_t pos, const VariantQueryConfig& query_config,
        const ITER* tile_its, uint64_t* num_deref_tile_iters) const;
    /** 
     * Initializes tile iterators for joint genotyping for column col. 
     * Returns the number of attributes used in joint genotyping.
     */
    unsigned int gt_initialize_tile_its(
        const StorageManager::ArrayDescriptor* ad,
        const VariantQueryConfig& query_config, const unsigned col_range_idx,
        StorageManager::const_reverse_iterator*& tile_its,
        StorageManager::const_reverse_iterator& tile_it_end ) const;
    /**
     * Helper function to fill the attribute given the tile pointer,
     * position, and index
     */
    template<class ITER>
    void fill_cell_attribute(const int64_t& pos, const ITER* tile_its, 
            uint64_t* num_deref_tile_iters,
            const unsigned IDX, int *p_int_v) const;
    /**
     * Override of the function above for float type
     */
    template<class ITER>
    void fill_cell_attribute(int64_t pos, const ITER* tile_its, 
            uint64_t* num_deref_tile_iters,
            unsigned IDX, float *p_float_v) const;
    /**
     * Variables to store versioning information about array schema
     */
    unsigned m_GT_schema_version;
    /**
     * Map the known field enum to cell attribute idx for the given schema
     */
    SchemaIdxToKnownVariantFieldsEnumLUT m_schema_idx_to_known_variant_field_enum_LUT;
    /**
     * Vector that stores information about the known fields - NULL bitidx, OFFSETS bitidx, length etc
     */
    std::vector<VariantQueryFieldInfo> m_known_field_enum_to_info;
    /*
     * Static members that track information known about variant data
     */
    //All known field names specific to variant data
    static std::vector<std::string> m_known_variant_field_names;
    //Mapping from field name to enum idx
    static std::unordered_map<std::string, unsigned> m_known_variant_field_name_to_enum;
    //Flag to check whether static members are initialized
    static bool m_are_static_members_initialized;
    //Function that initializes static members
    static void initialize_static_members();
};


#endif
