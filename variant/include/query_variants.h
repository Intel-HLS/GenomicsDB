#ifndef QUERY_VARIANTS_H
#define QUERY_VARIANTS_H

#include "gt_common.h"
#include "variant_storage_manager.h"
#include "variant_array_schema.h"
#include "variant_query_config.h"
#include "variant.h"
#include "variant_operations.h"
#include "variant_cell.h"

enum GTSchemaVersionEnum
{
    GT_SCHEMA_V0=0,
    GT_SCHEMA_V1,
    GT_SCHEMA_V2
};

//Exceptions thrown 
class InconsistentQueryOptionsException : public std::exception {
  public:
    InconsistentQueryOptionsException(const std::string m="") : msg_("InconsistentQueryOptionsException : "+m) { ; }
    ~InconsistentQueryOptionsException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
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
class VariantQueryProcessor {

  public:
    /** 
     * Simple constructor. The workspace is where the query processor will create 
     * its data. The storage manager is the module the query processor interefaces 
     * with.
     */
    VariantQueryProcessor(VariantStorageManager* storage_manager, const std::string& array_name);
    /*
     * This constructor is useful when there is no VariantStorageManager. This happens when loading and querying are 
     * combined and the data does not exist on disk as TileDB arrays
     */
    VariantQueryProcessor(const VariantArraySchema& array_schema);
    ~VariantQueryProcessor()
    {
      if(m_array_schema)
        delete m_array_schema;
      m_array_schema = 0;
    }
    void clear();
    /**
     * When querying, setup bookkeeping structures first 
     */
    void do_query_bookkeeping(const VariantArraySchema& array_schema,
        VariantQueryConfig& query_config) const;
    /*
     * Equivalent of gt_get_column, but for interval
     */
    void gt_get_column_interval(
        const int ad,
        const VariantQueryConfig& query_config, unsigned column_interval_idx,
        std::vector<Variant>& variants, GA4GHPagingInfo* paging_info=0, GTProfileStats* stats=0) const;
    /*
     * Scans column interval, aligns intervals and runs operate
     */
    void scan_and_operate(const int ad, const VariantQueryConfig& query_config,
        SingleVariantOperatorBase& variant_operator,
        unsigned column_interval_idx=0u, bool handle_spanning_deletions=false) const;
    /*
     * Deal with next cell in forward iteration in a scan
     * */
    bool scan_handle_cell(const VariantQueryConfig& query_config, unsigned column_interval_idx,
        Variant& variant, SingleVariantOperatorBase& variant_operator,
        BufferVariantCell& cell, const void* cell_ptr,
        VariantCallEndPQ& end_pq, std::vector<VariantCall*>& tmp_pq_buffer,
        int64_t& current_start_position, int64_t& next_start_position,
        uint64_t& num_calls_with_deletions, bool handle_spanning_deletions,
        GTProfileStats* stats_ptr) const;
    /** Called by scan_and_operate to handle all ranges for given set of cells */
    void handle_gvcf_ranges(VariantCallEndPQ& end_pq, 
        const VariantQueryConfig& queryConfig, Variant& variant,
        SingleVariantOperatorBase& variant_operator,
        int64_t current_start_position, int64_t next_start_position, bool is_last_call, uint64_t& num_calls_with_deletions) const;
    //while scan breaks up the intervals, iterate does not
    void iterate_over_cells(
        const int ad,
        const VariantQueryConfig& query_config, 
        SingleCellOperatorBase& variant_operator, unsigned column_interval_idx) const;
    /** Fills genotyping info for column col from the input array. */
    //Row ordering vector stores the query row idx in the order in which rows were filled by gt_get_column function
    //This is the reverse of the cell position order (as reverse iterators are used in gt_get_column)
    void gt_get_column(
        const int ad, const VariantQueryConfig& query_config, unsigned column_interval_idx,
        Variant& variant, GTProfileStats* stats=0, std::vector<uint64_t>* query_row_idx_in_order=0) const;
    /*
     * Create Variant from a buffer produced by the binary_serialize() functions
     * This is a member of VariantQueryProcessor because the Factory methods are already setup for 
     * constructing various types of fields.
     */
    void binary_deserialize(Variant& variant, const VariantQueryConfig& query_config,
        const std::vector<uint8_t>& buffer, uint64_t& offset) const;
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
    int get_array_descriptor() const { return m_ad; }
    const VariantArraySchema& get_array_schema() const { return *m_array_schema; }
    /**
     * A function that obtains cell attribute idxs for queried attribute names in the queryConfig object
     */
    void obtain_TileDB_attribute_idxs(const VariantArraySchema& array_schema, VariantQueryConfig& queryConfig) const;
    /**
     * Return VariantStorageManager object
     */
    const VariantStorageManager* get_storage_manager() const { return m_storage_manager; }
  private:
    /*initialize all known info about variants*/
    void initialize_known(const VariantArraySchema& array_schema);
    /*Initialize schema version v0 info*/
    void initialize_v0(const VariantArraySchema& array_schema);
    /*Check and initialize schema version v1 info*/
    void initialize_v1(const VariantArraySchema& array_schema);
    /*Check and initialize schema version v2 info*/
    void initialize_v2(const VariantArraySchema& array_schema);
    /*Initialize versioning information based on schema*/
    void initialize_version(const VariantArraySchema& array_schema);
    /*Register field creator pointers with the factory object*/
    void register_field_creators(const VariantArraySchema& array_schema);
    /*Wrapper for initialize functions - assumes m_array_schema is initialized correctly*/
    void initialize();
    /**
     * Initialized field length info
     */
    static void initialize_length_descriptor(unsigned idx);
    /*
     * Fills a row of the input genotyping column with the proper info.
     * When cells are duplicated at the END, the flag traverse_end_copies controls whether END copies should be
     * traversed or begin copies
     */
    void gt_fill_row(
        Variant& variant, int64_t row, int64_t column, const VariantQueryConfig& query_config,
        const BufferVariantCell& cell, GTProfileStats* stats
#ifdef DUPLICATE_CELL_AT_END
        , bool traverse_end_copies=false
#endif
        ) const;
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
    void fill_field(std::unique_ptr<VariantFieldBase>& field_ptr, const BufferVariantCell::FieldsIter& attr_iter,
        const unsigned num_ALT_alleles, const unsigned ploidy,
        const unsigned schema_idx
        ) const;
    /*
     * Prep work for filling out a field_ptr
     * Allocates field_ptr using the Factory
     * Marks field_ptr as valid
     * Gets length_descriptor, num_elements for known fields
     */
    void fill_field_prep(std::unique_ptr<VariantFieldBase>& field_ptr, unsigned schema_idx,
        unsigned& length_descriptor, unsigned& num_elements) const;
    /*
     * VariantStorage manager
     */
    const VariantStorageManager* m_storage_manager;
    /**
     * Variables to store versioning information about array schema
     */
    unsigned m_GT_schema_version;
    /**
     * Map the known field enum to cell attribute idx for the given schema
     */
    SchemaIdxToKnownVariantFieldsEnumLUT m_schema_idx_to_known_variant_field_enum_LUT;
    /**
     * Factory object that creates variant fields as and when needed
     */
    VariantFieldFactory m_field_factory;
    /*
     * Array descriptor and schema
     */
    int m_ad;
    VariantArraySchema* m_array_schema;
    //Mapping from std::type_index to VariantFieldCreator pointers, used when schema loaded to set creators for each attribute
    static std::unordered_map<std::type_index, std::shared_ptr<VariantFieldCreatorBase>> m_type_index_to_creator;
    //Flag to check whether static members are initialized
    static bool m_are_static_members_initialized; 
    //Function that initializes static members
    static void initialize_static_members();
};


#endif
