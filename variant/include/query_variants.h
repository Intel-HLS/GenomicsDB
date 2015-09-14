#ifndef QUERY_VARIANTS_H
#define QUERY_VARIANTS_H

#include "query_processor.h"
#include "gt_common.h"
#include "variant_query_config.h"
#include "variant.h"
#include "variant_operations.h"

enum GTSchemaVersionEnum
{
    GT_SCHEMA_V0=0,
    GT_SCHEMA_V1,
    GT_SCHEMA_V2
};

/*
 * Class that stores info about the some of the known fields
 */
class KnownFieldInfo
{
  public:
    KnownFieldInfo()
    {
      m_ploidy_required = false;
      m_length_descriptor = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_num_elements = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    bool m_ploidy_required;
    unsigned m_length_descriptor;
    unsigned m_num_elements;
    std::shared_ptr<VariantFieldCreatorBase> m_field_creator;
  public:
    inline bool is_length_allele_dependent() const
    {
      unsigned length_descriptor = m_length_descriptor;
      return (length_descriptor == BCF_VL_A || length_descriptor == BCF_VL_R || length_descriptor == BCF_VL_G);
    }
    inline unsigned get_length_descriptor() const { return m_length_descriptor; }
    inline bool is_length_genotype_dependent() const { return m_length_descriptor == BCF_VL_G; }
    inline bool is_length_only_ALT_alleles_dependent() const { return m_length_descriptor == BCF_VL_A; }
    unsigned get_num_elements_for_known_field_enum(unsigned num_ALT_alleles, unsigned ploidy) const;
    inline bool ploidy_required_for_known_field_enum() const { return m_ploidy_required; }
};

/* Structure to store profiling information */
class GTProfileStats {
  public:
    enum GTStatIdx
    {
      GT_NUM_CELLS=0,   //total #cells traversed in the query, coord cells are accessed for every time
      GT_NUM_CELLS_IN_LEFT_SWEEP,//total #cells traversed in the query left sweep, coord cells are accessed for every time
      GT_NUM_VALID_CELLS_IN_QUERY,//#valid cells actually returned in query 
      GT_NUM_ATTR_CELLS_ACCESSED,//#attribute cells accessed in the query
      GT_NUM_STATS
    };
    GTProfileStats();
    inline void update_stat(unsigned stat_idx, uint64_t value)
    {
      assert(stat_idx < m_stats_sum_vector.size());
      m_stats_sum_vector[stat_idx] += value;
      double v = value;
      m_stats_sum_sq_vector[stat_idx] += (v*v);
    }
    inline void update_stat_from_tmp_counter(unsigned stat_idx)
    {
      assert(stat_idx < m_stats_tmp_count_vector.size());
      update_stat(stat_idx, m_stats_tmp_count_vector[stat_idx]);
      m_stats_tmp_count_vector[stat_idx] = 0;
    }
    inline void update_all_from_tmp_counters()
    {
      for(auto stat_idx=0u;stat_idx<GT_NUM_STATS;++stat_idx)
      {
        assert(stat_idx < m_stats_tmp_count_vector.size());
        update_stat(stat_idx, m_stats_tmp_count_vector[stat_idx]);
        m_stats_tmp_count_vector[stat_idx] = 0;
      }
    }
    inline void increment_tmp_counter(unsigned stat_idx, uint64_t value=1u)
    {
      assert(stat_idx < m_stats_tmp_count_vector.size());
      m_stats_tmp_count_vector[stat_idx] += value;
    }
    inline void reset_tmp_counter(unsigned stat_idx)
    {
      assert(stat_idx < m_stats_tmp_count_vector.size());
      m_stats_tmp_count_vector[stat_idx] = 0;
    }
    void print_stats(std::ostream& fptr=std::cout) const;
    void increment_num_queries() { ++m_num_queries; }
  private:
    std::vector<uint64_t> m_stats_sum_vector;
    std::vector<uint64_t> m_stats_tmp_count_vector;
    std::vector<double> m_stats_sum_sq_vector;
    std::vector<std::string> m_stats_name_vector;
    uint64_t m_num_queries;
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
    VariantQueryProcessor(StorageManager* storage_manager, const std::string& array_name);
    void clear();
    /**
     * When querying, setup bookkeeping structures first 
     */
    void do_query_bookkeeping(const ArraySchema& array_schema,
        VariantQueryConfig& query_config) const;
    /*
     * Equivalent of gt_get_column, but for interval
     */
    void gt_get_column_interval(
        const int ad,
        const VariantQueryConfig& query_config, unsigned column_interval_idx,
        std::vector<Variant>& variants, GTProfileStats* stats=0) const;
#if 0
    void scan_and_operate(const ArraySchema& array_schema, const VariantQueryConfig& query_config,
        SingleVariantOperatorBase& variant_operator,
        unsigned column_interval_idx=0u) const;
    void iterate_over_all_tiles(const ArraySchema& array_schema, const VariantQueryConfig& query_config) const;
#endif
    /** Fills genotyping info for column col from the input array. */
    //Row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
    //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
    void gt_get_column(
        const int ad, const VariantQueryConfig& query_config, unsigned column_interval_idx,
        Variant& variant, GTProfileStats* stats=0, std::vector<uint64_t>* query_row_idx_in_order=0) const;
    /*
     * Functions that determine number of elements for known fields
     */
    inline bool is_length_allele_dependent(unsigned enumIdx) const
    {
      assert(enumIdx < m_known_field_enum_to_info.size());
      return m_known_field_enum_to_info[enumIdx].is_length_allele_dependent();
    }
    unsigned get_num_elements_for_known_field_enum(unsigned enumIdx, unsigned num_ALT_alleles, unsigned ploidy) const;
    unsigned get_length_descriptor_for_known_field_enum(unsigned known_field_enum) const;
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
    /*
     * Check whether the known field requires ploidy - e.g. GT, GQ etc
     */
    inline bool ploidy_required_for_known_field_enum(unsigned enumIdx) const
    {
      assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
      return m_known_field_enum_to_info[enumIdx].ploidy_required_for_known_field_enum();
    }
    int get_array_descriptor() const { return m_ad; }
    const ArraySchema& get_array_schema() const { return *m_array_schema; }
    /*
     * Given a field name, checks for m_known_variant_field_name_to_enum to see if valid entry exists.
     * If yes, fills known_field_enum and returns true
     * Else returns false. Leaves known_field_enum untouched.
     */
    static bool get_known_field_enum_for_name(const std::string& field_name, unsigned& known_field_enum);
    /*
     * Get name for known field enum
     */
    static std::string get_known_field_name_for_enum(unsigned known_field_enum);
  private:
    /*initialize all known info about variants*/
    void initialize_known(const ArraySchema& array_schema);
    /*Initialize schema version v0 info*/
    void initialize_v0(const ArraySchema& array_schema);
    /*Check and initialize schema version v1 info*/
    void initialize_v1(const ArraySchema& array_schema);
    /*Check and initialize schema version v2 info*/
    void initialize_v2(const ArraySchema& array_schema);
    /*Initialize versioning information based on schema*/
    void initialize_version(const ArraySchema& array_schema);
    /*Register field creator pointers with the factory object*/
    void register_field_creators(const ArraySchema& array_schema);
    /**
     * Initialized field length info
     */
    void initialize_length_descriptor(unsigned idx);
    /** Called by scan_and_operate to handle all ranges for given set of cells */
    void handle_gvcf_ranges(VariantCallEndPQ& end_pq, 
        const VariantQueryConfig& queryConfig, Variant& variant,
        SingleVariantOperatorBase& variant_operator,
        int64_t current_start_position, int64_t next_start_position, bool is_last_call) const;
    /** Fills a row of the input genotyping column with the proper info. */
    void gt_fill_row(
        Variant& variant, int64_t row, int64_t column, const VariantQueryConfig& query_config,
        const Cell& cell, GTProfileStats* stats) const;
    /** 
     * Initializes reverse iterators for joint genotyping for column col. 
     * Returns the number of attributes used in joint genotyping.
     */
    unsigned int gt_initialize_reverse_iter(
        const int ad,
        const VariantQueryConfig& query_config, const int64_t column,
        ArrayConstReverseCellIterator<int64_t>*& reverse_iter) const;
    /** 
     * Initializes forward iterators for joint genotyping for column col. 
     * Returns the number of attributes used in joint genotyping.
     */
    unsigned int gt_initialize_forward_iter(
        const int ad,
        const VariantQueryConfig& query_config, const int64_t column,
        ArrayConstCellIterator<int64_t>*& forward_iter) const;
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
     */
    void fill_field(std::unique_ptr<VariantFieldBase>& field_ptr, const CellConstAttrIterator& attr_iter,
        const unsigned num_ALT_alleles, const unsigned ploidy,
        const unsigned schema_idx
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
     * Factory object that creates variant fields as and when needed
     */
    VariantFieldFactory m_field_factory;
    /*
     * Array descriptor and schema
     */
    int m_ad;
    const ArraySchema* m_array_schema;
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
