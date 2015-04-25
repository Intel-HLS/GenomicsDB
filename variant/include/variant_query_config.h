#ifndef VARIANT_QUERY_CONFIG_H
#define VARIANT_QUERY_CONFIG_H

#include "query_config.h"
#include "lut.h"

#define UNDEFINED_NUM_ROWS_VALUE 0xFFFFFFFFFFFFFFFFull

typedef std::pair<int64_t, int64_t> ColumnRange;
bool ColumnRangeCompare(const ColumnRange& x, const ColumnRange& y);

//Out of bounds query exception
class OutOfBoundsQueryException {
  public:
    OutOfBoundsQueryException(const std::string m="Out of bounds query exception") : msg_(m) { ; }
    ~OutOfBoundsQueryException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const std::string& what() const { return msg_; }
  private:
    std::string msg_;
};

class VariantQueryConfig : public QueryConfig
{
  public:
    VariantQueryConfig()
      : QueryConfig()
    {
      clear();
      m_query_idx_known_variant_field_enum_LUT.reset_luts();
      m_done_bookkeeping = false;
      m_query_all_rows = true;
      m_num_rows_in_array = UNDEFINED_NUM_ROWS_VALUE;
      m_first_normal_field_query_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    void clear()
    {
      QueryConfig::clear();
      m_query_rows.clear();
      m_query_column_intervals.clear();
      m_array_row_idx_to_query_row_idx.clear();
    }
    /*
     * Re-order query fields so that special fields like COORDS,END,NULL,OFFSET,ALT are first
     */
    void reorder_query_fields();
    unsigned get_first_normal_field_query_idx() const { return m_first_normal_field_query_idx; }
    //Query idx <--> known fields mapping
    void resize_LUT(unsigned num_known_fields)
    { m_query_idx_known_variant_field_enum_LUT.resize_luts_if_needed(get_num_queried_attributes(), num_known_fields); }
    //Map queryIdx <--> knownEnumIdx
    inline void add_query_idx_known_field_enum_mapping(unsigned queryIdx, unsigned knownEnumIdx)
    {
      assert(queryIdx < get_num_queried_attributes());
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      m_query_idx_known_variant_field_enum_LUT.add_query_idx_known_field_enum_mapping(queryIdx, knownEnumIdx);
    }
    //Get query idx for given knownEnumIdx
    inline unsigned get_query_idx_for_known_field_enum(unsigned knownEnumIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      return m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx);
    }
    //Check whether query contains given knownEnumIdx
    inline bool is_defined_query_idx_for_known_field_enum(unsigned knownEnumIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
        m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx));
    }
    inline bool is_bookkeeping_done() const { return m_done_bookkeeping; }
    inline void set_done_bookkeeping(bool value) { m_done_bookkeeping = value; }
    /*
     * Function that specifies which rows to query
     */
    void set_rows_to_query(const std::vector<int64_t>& rowIdx);
    /*
     * Used by QueryProcessor to set number of rows if all rows need to be queried.
     */
    void set_num_rows_in_array(uint64_t num_rows) { m_num_rows_in_array = num_rows; }
    /*
     * Initialize map between array row to query row
     */
    void setup_array_row_idx_to_query_row_idx_map();
    /**
     * Rows to query
     */
    inline bool query_all_rows() const { return m_query_all_rows; }
    inline const std::vector<int64_t>& get_rows_to_query() const { return m_query_rows; }
    /**
     * If all rows are queried, return m_num_rows_in_array (set by QueryProcessor)
     * Else return size of m_query_rows vector
     */
    inline uint64_t get_num_rows_to_query() const
    {
      /*Either query subset of rows (in m_query_rows) or set the variable m_num_rows_in_array correctly*/
      assert(!m_query_all_rows || (m_num_rows_in_array != UNDEFINED_NUM_ROWS_VALUE));
      return m_query_all_rows ? m_num_rows_in_array : m_query_rows.size();
    }
    inline uint64_t get_num_rows_in_array() const
    {
      assert(m_num_rows_in_array != UNDEFINED_NUM_ROWS_VALUE);
      return m_num_rows_in_array;
    }
    /**
     * If all rows are queried, return idx
     * Else return m_query_rows[idx]
     */
    inline int64_t get_array_row_idx_for_query_row_idx(uint64_t idx) const
    {
      assert(idx < get_num_rows_to_query());
      return m_query_all_rows ? idx : m_query_rows[idx];
    }
    /*
     * Index in m_query_rows for given array row idx
     */
    inline uint64_t get_query_row_idx_for_array_row_idx(int64_t row_idx) const
    {
      assert(row_idx >= 0 && row_idx < get_num_rows_in_array());
      if(m_query_all_rows)
        return row_idx;
      assert(row_idx < m_array_row_idx_to_query_row_idx.size());
      return m_array_row_idx_to_query_row_idx[row_idx];
    }
    /*
     * Check if this row is being queried or no
     */
    inline bool is_queried_array_row_idx(int64_t row_idx) const
    {
      return m_query_all_rows ? true : (get_query_row_idx_for_array_row_idx(row_idx) != UNDEFINED_NUM_ROWS_VALUE);
    }
    /*
     * Function that specifies which column ranges to query
     */
    void add_column_interval_to_query(const int64_t colBegin, const int64_t colEnd);
    /**
     * Returns number of ranges queried
     */
    inline unsigned get_num_column_intervals() const
    {
      return  m_query_column_intervals.size();
    }
    inline ColumnRange get_column_interval(unsigned idx) const
    {
      assert(idx < m_query_column_intervals.size());
      return m_query_column_intervals[idx];
    }
  private:
    //Flag that tracks whether book-keeping is done
    bool m_done_bookkeeping;
    //Idx in m_query_attributes_names with the first common attribute - see reorder_query_fields();
    unsigned m_first_normal_field_query_idx;
    //Mapping between queried idx and known fields enum
    QueryIdxToKnownVariantFieldsEnumLUT m_query_idx_known_variant_field_enum_LUT;
    /*Rows to query*/
    bool m_query_all_rows;
    /*vector of specific row idxs to query*/
    std::vector<int64_t> m_query_rows;
    /*vector mapping array row_idx to query row idx*/
    std::vector<uint64_t> m_array_row_idx_to_query_row_idx;
    /*Set by query processor*/
    uint64_t m_num_rows_in_array;
    /*Column ranges to query*/
    std::vector<ColumnRange> m_query_column_intervals;
};

#endif
