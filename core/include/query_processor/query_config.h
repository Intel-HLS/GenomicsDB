#ifndef QUERY_CONFIG_H
#define QUERY_CONFIG_H

#include <vector>
#include <string>
#include <unordered_map>
#include <assert.h>

#define UNDEFINED_ATTRIBUTE_IDX_VALUE 0xFFFFFFFFu

/**
 * Class to store the query request information
 */
class QueryConfig
{
  public:
    /**
     * Simple constructor
     */
    QueryConfig();
    void clear();
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
      return (m_query_attributes_schema_idxs[idx] != UNDEFINED_ATTRIBUTE_IDX_VALUE);
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
  protected:
    std::vector<std::string> m_query_attributes_names;
    std::vector<int> m_query_attributes_schema_idxs;
    std::unordered_map<std::string, unsigned> m_query_attribute_name_to_query_idx;
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

#endif
