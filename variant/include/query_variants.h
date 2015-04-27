#ifndef QUERY_VARIANTS_H
#define QUERY_VARIANTS_H

#include "query_processor.h"
#include "gt_common.h"
#include "variant_query_config.h"
#include "variant.h"

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
//Stores index in OFFSETS field for each field
enum KnownVariantFieldOffsetsEnum
{
  GVCF_REF_OFFSET_IDX=0,
  GVCF_ALT_OFFSET_IDX,
  GVCF_FILTER_OFFSET_IDX,
  GVCF_AD_OFFSET_IDX,
  GVCF_PL_OFFSET_IDX,
  GVCF_NUM_KNOWN_OFFSET_ELEMENTS_PER_CELL
};
/*
 * Class that stores info about the some of the known fields
 */
class KnownFieldInfo
{
  public:
    KnownFieldInfo()
    {
      m_NULL_bitidx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_OFFSETS_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_length_descriptor = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_num_elements = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    unsigned m_NULL_bitidx;
    unsigned m_OFFSETS_idx;
    unsigned m_length_descriptor;
    unsigned m_num_elements;
    std::shared_ptr<VariantFieldCreatorBase> m_field_creator;
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
        VariantQueryConfig& query_config) const;
    /*
     * Equivalent of gt_get_column, but for interval
     */
    void gt_get_column_interval(
        const StorageManager::ArrayDescriptor* ad,
        const VariantQueryConfig& query_config, unsigned column_interval_idx,
        std::vector<Variant>& variants, GTProfileStats* stats) const;
    /** Fills genotyping info for column col from the input array. */
    //Row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
    //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
    void gt_get_column(
        const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config, unsigned column_interval_idx,
        Variant& variant, GTProfileStats* stats=0, std::vector<uint64_t>* query_row_idx_in_order=0) const;
    void scan_and_operate(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config, std::ostream& output_stream) const;
    void iterate_over_all_tiles(const StorageManager::ArrayDescriptor* ad, const VariantQueryConfig& query_config) const;
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum,
     * checks if the field requires NULL bitidx field to be used 
     */
    inline bool is_NULL_bitidx_defined_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size());
      return (m_known_field_enum_to_info[enumIdx].m_NULL_bitidx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
    }
    /*
     * Client code MUST check is_NULL_bitidx_defined_for_known_field_enum() before using the
     * value returned by this function
     */
    inline unsigned get_NULL_bitidx_for_known_field_enum(unsigned enumIdx) const
    {
        assert(enumIdx < m_known_field_enum_to_info.size());
        return (m_known_field_enum_to_info[enumIdx].m_NULL_bitidx);
    }
    /*
     * Functions that determine number of elements for known fields
     */
    inline bool is_length_allele_dependent(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size());
      unsigned length_descriptor = m_known_field_enum_to_info[enumIdx].m_length_descriptor;
      return (length_descriptor == BCF_VL_A || length_descriptor == BCF_VL_R || length_descriptor == BCF_VL_G);
    }
    unsigned get_num_elements_for_known_field_enum(unsigned enumIdx, unsigned num_ALT_alleles) const;
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum
     * returns true if the field requires the OFFSETS field 
     */
    inline bool uses_OFFSETS_field_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size()); 
      return (m_known_field_enum_to_info[enumIdx].m_OFFSETS_idx != UNDEFINED_ATTRIBUTE_IDX_VALUE);
    }
    /*
     * Client code MUST check uses_OFFSETS_field_for_known_field_enum() before using the
     * value returned by this function
     */
    inline unsigned get_OFFSETS_idx_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size()); 
      return (m_known_field_enum_to_info[enumIdx].m_OFFSETS_idx);
    }
    /*
     * Function that, given an enum value from KnownVariantFieldsEnum
     * returns the schema idx for the given array 
     * NOTE: returned value may be invalid, client code MUST check validity using
     * m_schema_idx_to_known_variant_field_enum_LUT.is_defined_value()
     */
    inline unsigned get_schema_idx_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
      return m_schema_idx_to_known_variant_field_enum_LUT.get_schema_idx_for_known_field_enum(enumIdx);
    }
    /*
     * Check whether the known field requires a special creator
     */
    inline bool requires_special_creator(unsigned enumIdx) const
    {
      assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
      return (m_known_field_enum_to_info[enumIdx].m_field_creator.get() != 0);
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
    /*Register field creator pointers with the factory object*/
    void register_field_creators(const StorageManager::ArrayDescriptor* ad);
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
    void initialize_length_descriptor(unsigned idx);
    /** Called by scan_and_operate to handle all ranges for given set of cells */
    void handle_gvcf_ranges(VariantCallEndPQ& end_pq, 
        const VariantQueryConfig& queryConfig, Variant& variant,
        std::ostream& output_stream,
        int64_t current_start_position, int64_t next_start_position, bool is_last_call) const;
    /** Fills a row of the input genotyping column with the proper info. */
    template<class ITER>
    void gt_fill_row(
        Variant& variant, int64_t row, int64_t column, int64_t pos, const VariantQueryConfig& query_config,
        const ITER* tile_its, uint64_t* num_deref_tile_iters) const;
    /** 
     * Initializes tile iterators for joint genotyping for column col. 
     * Returns the number of attributes used in joint genotyping.
     */
    unsigned int gt_initialize_tile_its(
        const StorageManager::ArrayDescriptor* ad,
        const VariantQueryConfig& query_config, const unsigned column_interval_idx,
        StorageManager::const_reverse_iterator*& tile_its,
        StorageManager::const_reverse_iterator& tile_it_end ) const;
    /*
     * Fill data from tile for attribute query_idx into curr_call
     * @param curr_call  VariantCall object in which data will be stored
     * @param tile AttributeTile object from which data needs to be copied
     * @param pos Cell position in the co-ordinates tile
     * The following 3 parameters MUST be valid if the attribute being accessed needs them
     * Else, it's ok to leave them NULL, 0 etc
     * @param OFFSETS_values- can be nullptr, if the current field does not use OFFSETS
     * @param NULL_bitmap - can be 0, if the current field does not use NULL field
     * @param num_ALT_alleles - number of ALT alleles:can be 0, if the current field does not use this info
     * @schema_idx The idx of the attribute in the schema of the current array 
     * @num_deref_tile_iters - useful only if compiled in PROFILING mode
     */
    void fill_field(std::unique_ptr<VariantFieldBase>& field_ptr,
        const Tile& tile, int64_t pos,
        const std::vector<int64_t>* OFFSETS_values, const int NULL_bitmap, const unsigned num_ALT_alleles,
        const unsigned schema_idx, uint64_t* num_deref_tile_iters=0
        ) const;

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
    std::vector<KnownFieldInfo> m_known_field_enum_to_info;
    /**
     * Number of elements per cell in the OFFSETS field - could change with different schema versions
     */
    unsigned m_num_elements_per_offset_cell;
    /**
     * Factory object that creates variant fields as and when needed
     */
    VariantFieldFactory m_field_factory;
    /*
     * Static members that track information known about variant data
     */
    //All known field names specific to variant data
    static std::vector<std::string> m_known_variant_field_names;
    //Mapping from field name to enum idx
    static std::unordered_map<std::string, unsigned> m_known_variant_field_name_to_enum;
    //Mapping from std::type_index to VariantFieldCreator pointers, used when schema loaded to set creators for each attribute
    static std::unordered_map<std::type_index, std::shared_ptr<VariantFieldCreatorBase>> m_type_index_to_creator;
    //Flag to check whether static members are initialized
    static bool m_are_static_members_initialized;
    //Function that initializes static members
    static void initialize_static_members();
};


#endif
