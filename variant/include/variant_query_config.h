#ifndef VARIANT_QUERY_CONFIG_H
#define VARIANT_QUERY_CONFIG_H

#include "headers.h"
#include "lut.h"
#include "known_field_info.h"

//Out of bounds query exception
class OutOfBoundsQueryException : public std::exception {
  public:
    OutOfBoundsQueryException(const std::string m="Out of bounds query exception") : msg_(m) { ; }
    ~OutOfBoundsQueryException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class UnknownQueryAttributeException {
  public:
    UnknownQueryAttributeException(const std::string m="Invalid queried attribute") : msg_(m) { ; }
    ~UnknownQueryAttributeException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const std::string& what() const { return msg_; }
  private:
    std::string msg_;
};

class VariantQueryConfig
{
  public:
    VariantQueryConfig()
    {
      clear();
      m_query_idx_known_variant_field_enum_LUT.reset_luts();
      m_done_bookkeeping = false;
      m_query_all_rows = true;
      m_num_rows_in_array = UNDEFINED_NUM_ROWS_VALUE;
      m_smallest_row_idx = 0;
      m_first_normal_field_query_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    void clear()
    {
      m_query_attributes_names.clear();
      m_query_attributes_schema_idxs.clear();
      m_query_attribute_name_to_query_idx.clear();
      m_query_rows.clear();
      m_query_column_intervals.clear();
      m_array_row_idx_to_query_row_idx.clear();
    }
    /**
     * Function that specifies which attributes to query from each cell
     */
    void set_attributes_to_query(const std::vector<std::string>& attributeNames);
    /**
     * Function used by query processor to add extra attributes to query
     */
    void add_attribute_to_query(const std::string& name, unsigned schema_idx);
    /**
     * Check whether attribute is defined
     */
    inline bool is_schema_idx_defined_for_query_idx(unsigned idx) const
    {
      assert(idx < m_query_attributes_schema_idxs.size());
      return (m_query_attributes_schema_idxs[idx] != static_cast<int>(UNDEFINED_ATTRIBUTE_IDX_VALUE));
    }
    /**
     * Set TileDB array schema attribute idx (schemaIdx) for the queried attribute idx
     */
    void set_schema_idx_for_query_idx(unsigned idx, unsigned schemaIdx)
    {
      assert(idx < m_query_attributes_schema_idxs.size());
      m_query_attributes_schema_idxs[idx] = schemaIdx;
    }
    /**
     * Get TileDB array schema attribute idx for the queried attribute idx
     */
    inline unsigned get_schema_idx_for_query_idx(unsigned idx) const
    {
      assert(idx < m_query_attributes_schema_idxs.size());
      return m_query_attributes_schema_idxs[idx];
    }
    /**
     * Get idx in the query for given attribute name
     */
    inline bool get_query_idx_for_name(const std::string& name, unsigned& idx) const
    {
      auto iter = m_query_attribute_name_to_query_idx.find(name);
      if(iter != m_query_attribute_name_to_query_idx.end())
      {
        idx = (*iter).second;
        return true;
      }
      else
        return false;
    }
    /**
     * Get name for idx
     */
    inline std::string get_query_attribute_name(unsigned idx) const
    {
      assert(idx < m_query_attributes_schema_idxs.size());
      return m_query_attributes_names[idx];
    }
    /**
     * Get number of attributes in query
     */
    inline unsigned get_num_queried_attributes() const { return m_query_attributes_names.size(); }
    inline const std::vector<int>& get_query_attributes_schema_idxs() const { return m_query_attributes_schema_idxs; }
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
    //Get known field enum for given queryIdx
    inline unsigned get_known_field_enum_for_query_idx(unsigned queryIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(queryIdx));
      return m_query_idx_known_variant_field_enum_LUT.get_known_field_enum_for_query_idx(queryIdx);
    }
    //Check whether query contains given knownEnumIdx
    inline bool is_defined_query_idx_for_known_field_enum(unsigned knownEnumIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
        m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx));
    }
    //Check whether query idx is known field
    inline bool is_defined_known_field_enum_for_query_idx(unsigned queryIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(queryIdx));
      return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
        m_query_idx_known_variant_field_enum_LUT.get_known_field_enum_for_query_idx(queryIdx));
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
    void set_num_rows_in_array(uint64_t num_rows, const uint64_t smallest_row_idx)
    {
      m_num_rows_in_array = num_rows;
      m_smallest_row_idx = smallest_row_idx;
    }
    /*
     * Initialize map between array row to query row
     */
    void setup_array_row_idx_to_query_row_idx_map();
    /*
     * Function that specifies which rows to query and also updates 
     * m_array_row_idx_to_query_row_idx
     * Pre-requisite: query bookkeeping should be done before calling this function
     */
    void update_rows_to_query(const std::vector<int64_t>& rowIdx);
    /*
     * Function that specifies all rows should be queried. 
     * m_array_row_idx_to_query_row_idx
     * Pre-requisite: query bookkeeping should be done before calling this function
     */
    void update_rows_to_query_to_all_rows();
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
    inline int64_t get_smallest_row_idx_in_array() const { return m_smallest_row_idx; }
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
      return m_query_all_rows ? (idx+m_smallest_row_idx) : m_query_rows[idx];
    }
    /*
     * Index in m_query_rows for given array row idx
     */
    inline uint64_t get_query_row_idx_for_array_row_idx(int64_t row_idx) const
    {
      assert(row_idx >= m_smallest_row_idx && (row_idx-m_smallest_row_idx) < static_cast<int64_t>(get_num_rows_in_array()));
      if(m_query_all_rows)
        return row_idx - m_smallest_row_idx;
      assert((row_idx-m_smallest_row_idx) < static_cast<int64_t>(m_array_row_idx_to_query_row_idx.size()));
      return m_array_row_idx_to_query_row_idx[row_idx-m_smallest_row_idx];
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
     * Note that the order in which these intervals are queried is NOT the same order in which
     * the intervals are added. Book-keeping sorts the intervals in ascending order so that
     * intervals close by have a chance of re-using cached tiles in memory 
     */
    void add_column_interval_to_query(const int64_t colBegin, const int64_t colEnd);
    /*
     * Function that sets the interval as the only interval to be queried
     */
    void set_column_interval_to_query(const int64_t colBegin, const int64_t colEnd);
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
    inline uint64_t get_column_begin(unsigned idx) const { return get_column_interval(idx).first; }
    inline uint64_t get_column_end(unsigned idx) const { return get_column_interval(idx).second; }
    /*
     * Functions for dealing with known field info
     * Returns pointer in vector if known field, else null
     */
    const KnownFieldInfo* get_info_for_query_idx(unsigned idx) const
    {
      assert(idx < get_num_queried_attributes());
      unsigned known_field_idx = m_query_idx_known_variant_field_enum_LUT.get_known_field_enum_for_query_idx(idx);
      if(m_query_idx_known_variant_field_enum_LUT.is_defined_value(known_field_idx))
        return &(g_known_field_enum_to_info[known_field_idx]);
      else
        return 0;
    }
  private:
    /*
     * Function to invalid TileDB array row idx -> query row idx mapping
     * @param all_rows if true, invalidates all mappings, else invalidates mapping for rows in m_query_rows only
     */
    void invalidate_array_row_idx_to_query_row_idx_map(bool all_rows);
    //Query attribute names
    std::vector<std::string> m_query_attributes_names;
    //Query attribute schema idxs
    std::vector<int> m_query_attributes_schema_idxs;
    //Map from query name to index in m_query_attributes_names/m_query_attributes_schema_idxs
    std::unordered_map<std::string, unsigned> m_query_attribute_name_to_query_idx;
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
    int64_t m_smallest_row_idx;
    /*Column ranges to query*/
    std::vector<ColumnRange> m_query_column_intervals;
};

#endif
