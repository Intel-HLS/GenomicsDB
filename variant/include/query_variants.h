#ifndef QUERY_VARIANTS_H
#define QUERY_VARIANTS_H

#include "gt_common.h"
#include "variant.h"
#include "query_processor.h"

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
        const StorageManager::ArrayDescriptor* ad)
      : QueryProcessor(workspace, storage_manager)
    {
      initialize_version(ad);
    }
    void initialize_version(const StorageManager::ArrayDescriptor* ad);
    /** Returns the genotyping info for column col from the input array. */
    GTColumn* gt_get_column(
        const StorageManager::ArrayDescriptor* ad, uint64_t col, GTProfileStats* stats=0) const;
    void scan_and_operate(const StorageManager::ArrayDescriptor* ad, std::ostream& output_stream);
    void iterate_over_all_cells(const StorageManager::ArrayDescriptor* ad);
    /**
     * Field attribute idx in schema
     */
    unsigned GVCF_END_IDX;
    unsigned GVCF_REF_IDX;
    unsigned GVCF_ALT_IDX;
    unsigned GVCF_QUAL_IDX;
    unsigned GVCF_FILTER_IDX;
    unsigned GVCF_BASEQRANKSUM_IDX;
    unsigned GVCF_CLIPPINGRANKSUM_IDX;
    unsigned GVCF_MQRANKSUM_IDX;
    unsigned GVCF_READPOSRANKSUM_IDX;
    unsigned GVCF_DP_IDX;
    unsigned GVCF_MQ_IDX;
    unsigned GVCF_MQ0_IDX;
    unsigned GVCF_DP_FMT_IDX;
    unsigned GVCF_MIN_DP_IDX;
    unsigned GVCF_GQ_IDX;
    unsigned GVCF_SB_1_IDX;
    unsigned GVCF_SB_2_IDX;
    unsigned GVCF_SB_3_IDX;
    unsigned GVCF_SB_4_IDX;
    unsigned GVCF_AD_IDX;
    unsigned GVCF_PL_IDX;
    unsigned GVCF_AF_IDX;
    unsigned GVCF_AN_IDX;
    unsigned GVCF_AC_IDX;
    unsigned GVCF_NULL_IDX;
    unsigned GVCF_OFFSETS_IDX;
    unsigned GVCF_COORDINATES_IDX;

  private:
    /** Called by scan_and_operate to handle all ranges for given set of cells */
    void handle_gvcf_ranges(VariantIntervalPQ& end_pq, std::vector<PQStruct>& PQ_end_vec, GTColumn* gt_column,
        std::unordered_map<uint64_t, GTTileIteratorsTracker>& tile_idx_2_iters, std::ostream& output_stream,
        int64_t current_start_position, int64_t next_start_position, bool is_last_call);
    /** Fills a row of the input genotyping column with the proper info. */
    template<class ITER>
    void gt_fill_row(
        GTColumn* gt_column, int64_t row, int64_t column, int64_t pos,
        const ITER* tile_its, uint64_t* num_deref_tile_iters) const;
    /** 
     * Initializes tile iterators for joint genotyping for column col. 
     * Returns the number of attributes used in joint genotyping.
     */
    unsigned int gt_initialize_tile_its(
        const StorageManager::ArrayDescriptor* ad,
        StorageManager::const_reverse_iterator*& tile_its,
        StorageManager::const_reverse_iterator& tile_it_end,
        uint64_t col) const;
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
    unsigned m_PL_NULL_bitidx;
    /**
     * Queried attributes order
     */
    unsigned  GT_END_IDX;
    unsigned  GT_REF_IDX;
    unsigned  GT_ALT_IDX;
    unsigned  GT_PL_IDX;
    unsigned  GT_AF_IDX;
    unsigned  GT_AN_IDX;
    unsigned  GT_AC_IDX;
    unsigned  GT_NULL_IDX;
    unsigned  GT_OFFSETS_IDX;
    unsigned  GT_COORDINATES_IDX;
    unsigned  GT_IDX_COUNT;
};


#endif
