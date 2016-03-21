/**
 * Copyright (c) 2016 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef VARIANT_H
#define VARIANT_H

#include <assert.h>
#include "gt_common.h"
#include "variant_field_data.h"
#include "variant_query_config.h"
#include "known_field_info.h"
#include "vid_mapper.h"

/*
 * Class that stores info that helps determine whether 2 VariantCalls should be
 * merged into a single Variant object.
 * 2 VariantCalls x,y must satisfy the following conditions to be collapsed:
 * x.m_col_begin == y.m_col_begin && x.m_col_end == y.m_col.end && x.REF == y.REF
 * && x.ALT == y.ALT
 * This condition is checked using a hierarchical structure:
 * unordered_map : m_col_begin ->  
 * unordered_map : m_col_end ->
 * unordered_map : REF(string)->
 * map : set<ALT(string)> -> variant_idx
 * The ALTs are stored as set<string> since the ordering of ALT in two Calls may be different,
 * yet they may be identical
 */
class VariantCall;
class GA4GHCallInfoToVariantIdx 
{
  typedef std::map<std::set<std::string>, uint64_t> ALTSetToVariantIdxTy;
  typedef std::unordered_map<std::string, ALTSetToVariantIdxTy> REFToVariantIdxTy;
  typedef std::unordered_map<uint64_t, REFToVariantIdxTy> EndToVariantIdxTy;
  public:
    GA4GHCallInfoToVariantIdx() { clear(); }
    /*
     * This function finds an index corresponding to the given info if exists, else creates an entry 
     * in the hierarchical map with the new info and assigns the value of variant_index passed
     * Hence, variant_index must be set to a new index value when the function is called.
     * If the function determines the existence of a record with the same info, then variant_index is modified to
     * point to the existing record.
     * @return True if new record created, false otherwise
     */
    bool find_or_insert(uint64_t begin, uint64_t end, const std::string& REF, 
        const std::vector<std::string>& ALT_vec, uint64_t& variant_idx);
    bool find_or_insert(const VariantQueryConfig& query_config, VariantCall& to_move_call,
      uint64_t& variant_idx);
    void clear();
  private:
    std::unordered_map<uint64_t, EndToVariantIdxTy> m_begin_to_variant;
};

/**
 * Class equivalent to GACall in GA4GH API. Stores info about 1 CallSet/row for a given position
 */
class VariantCall
{
  public:
    /**
     * simple constructor - makes invalid VariantCall object
     */
    VariantCall()
    {
      m_is_valid = false;
      m_is_initialized = false;
      m_contains_deletion = false;
      m_row_idx = UNDEFINED_NUM_ROWS_VALUE;
      clear();
    }
    /**
     * Row in the TileDB array, NOT the idx of the VariantCall object in Variant::m_calls vector
     */
    VariantCall(uint64_t rowIdx)
    {
      m_is_valid = false;
      m_is_initialized = false;
      m_contains_deletion = false;
      m_row_idx = rowIdx;
      clear();
    } 
    /**
     * Default move constructor is good enough 
     */
    VariantCall(VariantCall&& other) = default;
    /**
     * Member element contains unique_ptr, cannot copy from const VariantCall&
     */
    VariantCall(const VariantCall& other) = delete;
    VariantCall(VariantCall& other) = delete;
    /*
     * Assignment with move semantics is fine, copy assignment not allowed
     */
    VariantCall& operator=(VariantCall&& other)
    {
      move_in(other);
      return *this;
    }
    /*
    * Free allocated memory
    */ 
    void clear()  {  m_fields.clear(); }        //also frees memory associated with m_fields elements (unique_ptr)
    /*
     * Same query_config, but new interval is starting. Reset what needs to be reset
     */
    void reset_for_new_interval();
    /*
     * Set TileDB array row index
     */
    void set_row_idx(uint64_t rowIdx) { m_row_idx = rowIdx; }
    uint64_t get_row_idx() const { return m_row_idx; }
    /*
     * Set genomic interval associated with this VariantCall
     * Could be different from the begin, end of the Variant that this Call is part of
     * because this Call is an interval that overlaps with the queried Variant
     */
    void set_column_interval(uint64_t col_begin, uint64_t col_end)
    {
      m_col_begin = col_begin;
      m_col_end = col_end;
    }
    uint64_t get_column_begin() const { return m_col_begin; }
    uint64_t get_column_end() const { return m_col_end; }
    /** 
     * Sometimes VariantCall might be allocated, but not store any valid data. This might happen if Variant allocates N
     * VariantCall objects, where N == number of rows. However, for a given query, not all VariantCall objects might contain valid data.
     * This could happen, if for a given queried location, row X does not have any valid data.
     * The m_is_valid flag specifies whether a this object contains valid data or not
     * Variant might pre-allocate N objects to avoid frequent de-allocations and re-allocations at run-time
     */
    void mark_valid(bool val) { m_is_valid = val; }
    bool is_valid() const { return m_is_valid; }
    /*
     * This variant call is not yet initialized by the query (if false)
     */
    void mark_initialized(bool val) { m_is_initialized = val; }
    bool is_initialized() const { return m_is_initialized; }
    /**
     * Functions dealing with field vector - self explanatory
     */
    void resize(unsigned num_fields)
    { m_fields.resize(num_fields);  }
    /**
     * Set field does a move transfers ownership of field data to member unique ptr. 
     */
    inline void set_field(unsigned idx, std::unique_ptr<VariantFieldBase>& field)
    {
      assert(idx < m_fields.size());
      m_fields[idx] = std::move(field); //transfer ownership of pointer
    }
    /*
     * Set field through raw pointer
     */
    inline void set_field(unsigned idx, VariantFieldBase* field)
    {
      assert(idx < m_fields.size());
      m_fields[idx] = std::move(std::unique_ptr<VariantFieldBase>(field)); //transfer ownership of pointer
    }
    void add_field(std::unique_ptr<VariantFieldBase>& field)
    {
      m_fields.push_back(std::move(field));
    }
    /*
     * Get reference to vector of fields
     */
    inline std::vector<std::unique_ptr<VariantFieldBase>>& get_all_fields()  {  return m_fields; }
    inline const std::vector<std::unique_ptr<VariantFieldBase>>& get_all_fields()  const {  return m_fields; }
    inline unsigned get_num_fields() const { return m_fields.size(); }
    /*
     * Get field at idx
     */
    inline std::unique_ptr<VariantFieldBase>& get_field(unsigned idx)
    {
      assert(idx < m_fields.size());
      return m_fields[idx];
    }
    inline const std::unique_ptr<VariantFieldBase>& get_field(unsigned idx) const
    {
      assert(idx < m_fields.size());
      return m_fields[idx];
    }
    template<class VariantFieldTy>
    inline VariantFieldTy* get_field(unsigned idx)
    {
      std::unique_ptr<VariantFieldBase>& smart_ptr = get_field(idx);
      auto* raw_ptr = smart_ptr.get();
      //Either the pointer is NULL itself, else make sure the correct subclass is produced
      assert(raw_ptr == 0 || dynamic_cast<VariantFieldTy*>(raw_ptr));
      return static_cast<VariantFieldTy*>(raw_ptr);
    }
    template<class VariantFieldTy>
    inline const VariantFieldTy* get_field(unsigned idx) const
    {
      const std::unique_ptr<VariantFieldBase>& smart_ptr = get_field(idx);
      const auto* raw_ptr = smart_ptr.get();
      //Either the pointer is NULL itself, else make sure the correct subclass is produced
      assert(raw_ptr == 0 || dynamic_cast<const VariantFieldTy*>(raw_ptr));
      return static_cast<const VariantFieldTy*>(raw_ptr);
    }
    /** print **/
    void print(std::ostream& stream, const VariantQueryConfig* query_config=0, const std::string& indent_prefix="") const;
    /* 
     * Print JSON for Cotton
     */
    void print_Cotton_JSON(std::ostream& fptr, unsigned field_idx) const;
    /*
     * Binary serialize into buffer
     */
    void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const;
    /*
     * Deserialize header of VariantCall (column interval, num fields etc)
     */
    void binary_deserialize_header(const std::vector<uint8_t>& buffer, uint64_t& offset);
    /**
     * Deep copy VariantCall, avoid using as much as possible (performance)
     */
    void copy_from_call(const VariantCall& src);
    /*
     * Flag that determines whether the current call contains a deletion
     */
    void set_contains_deletion(bool val) { m_contains_deletion = val; }
    bool contains_deletion() const { return m_is_valid && m_contains_deletion; }
    /*
     * Only copy simple elements and resize m_fields vector
     * Do not copy the fields themselves
     */
    void deep_copy_simple_members(const VariantCall& src);
  private:
    /*
     * Performs move from other object
     */
    void move_in(VariantCall& other);
    /*
     * Copies the simple member elements
     */
    void copy_simple_members(const VariantCall& other);
    /*
     * Member data elements - check clear, copy, move_in,binary_serialize/deserialize functions while adding new members
     */
    //Could be initialized, but invalid (no data for this column interval)
    bool m_is_valid;
    //If false, not initialized (not yet considered in query)
    bool m_is_initialized;
    //whether the ALT contains a deletion
    bool m_contains_deletion;
    uint64_t m_row_idx;
    std::vector<std::unique_ptr<VariantFieldBase>> m_fields;
    /**
     * Begin,end of this VariantCall
     * Could be different from the begin, end of the Variant that this Call is part of
     * because this Call is an interval that overlaps with the queried Variant
     **/
    uint64_t m_col_begin;
    uint64_t m_col_end;
};

class GA4GHPagingInfo;
/*
 * Class equivalent to GAVariant in GA4GH API. Stores info about 1 position/interval
 */
class Variant
{
  public:
    /*
     * Iterator object that traverses over VariantCalls vector, skipping invalid/un-initialized Calls
     * Templated for const, non-const iterators
     */
    template<class VariantCallTy, class IteratorTy>
    class ValidVariantCallIter
    {
      public:
        ValidVariantCallIter(IteratorTy x, IteratorTy end, 
            uint64_t call_idx_in_variant) 
          : m_iter_position(x), m_end(end), m_call_idx_in_variant(call_idx_in_variant)
        { 
          //If iter points to invalid Call, move forward
          if(m_iter_position != m_end && !((*m_iter_position).is_valid()))
            operator++();
        }
        bool operator!=(const ValidVariantCallIter& other) const { return m_iter_position != other.m_iter_position; }
        VariantCallTy& operator*() const { return *m_iter_position; }
        const ValidVariantCallIter& operator++()
        {
          ++m_iter_position;
          ++m_call_idx_in_variant;
          //Increment till end or next valid record
          for(;m_iter_position != m_end && !((*m_iter_position).is_valid());++m_iter_position,++m_call_idx_in_variant);
          return *this;
        }
        uint64_t get_call_idx_in_variant() const { return m_call_idx_in_variant; } 
      private:
        IteratorTy m_iter_position;
        IteratorTy m_end;
        uint64_t m_call_idx_in_variant;
    };
    using const_valid_calls_iterator = ValidVariantCallIter<const VariantCall, std::vector<VariantCall>::const_iterator>;
    using valid_calls_iterator = ValidVariantCallIter<VariantCall, std::vector<VariantCall>::iterator>;
  public:
    /*
     * Simple constructor
     */
    Variant()
    {
      m_query_config = 0;
      m_col_begin = m_col_end = UNDEFINED_NUM_ROWS_VALUE;
      clear();
    }
    /*
     * Constructor that takes the query_config object in whose response
     * this Variant will be created
     */
    Variant(const VariantQueryConfig* query_config)
    {
      m_query_config = query_config;
      m_col_begin = m_col_end = UNDEFINED_NUM_ROWS_VALUE;
      clear();
    }
    /*
     * Set genomic interval associated with this variant
     */
    void set_column_interval(uint64_t col_begin, uint64_t col_end)
    {
      m_col_begin = col_begin;
      m_col_end = col_end;
    }
    uint64_t get_column_begin() const { return m_col_begin; }
    uint64_t get_column_end() const { return m_col_end; }
    /**
     * Members of m_calls contains unique_ptr, cannot copy from const VariantCall&
     */
    Variant(const Variant& other) = delete;
    Variant(Variant& other) = delete;
    /*
     * Default move constructor is fine
     */
    Variant(Variant&& other) = default;
    /*
     * Move assignment operator is fine, no copy assignment allowed
     */
    Variant& operator=(Variant&& other)
    {
      move_in(other);
      return *this;
    } 
    /*
     * Memory de-allocation
     */
    void clear();
    /*
     * Same query_config, but new interval is starting. Reset what needs to be reset
     */
    void reset_for_new_interval();
    /* Set query config ptr */
    void set_query_config(const VariantQueryConfig* query_config) { m_query_config = query_config; }
    /**
     * Allocate maximum possible #calls and #fields per call based on m_query_config
     */
    void resize_based_on_query();
    /*
     * Resize call vector
     */
    void resize(uint64_t num_calls, unsigned num_query_call_fields)
    {
      m_calls.resize(num_calls);
      for(uint64_t i=0ull;i<num_calls;++i)
        m_calls[i].resize(num_query_call_fields);
    }
    /**
     * Append call to m_calls
     * @param call Only rvalues allowed as move done 
     */
    void add_call(VariantCall&& call)
    {
      m_calls.emplace_back(std::move(call));
    }
    /**
     * Create call object with m_row_idx = rowIdx
     */
    void add_call(const uint64_t rowIdx)
    {
      m_calls.emplace_back(VariantCall(rowIdx));
    }
    inline uint64_t get_num_calls() const { return m_calls.size(); }
    /**
     * Return VariantCall at index call_idx
     */
    inline VariantCall& get_call(uint64_t call_idx)
    {
      assert(call_idx < m_calls.size());
      return m_calls[call_idx];
    }
    inline const VariantCall& get_call(uint64_t call_idx) const
    {
      assert(call_idx < m_calls.size());
      return m_calls[call_idx];
    }
    inline std::vector<VariantCall>& get_calls() { return m_calls; }
    inline const std::vector<VariantCall>& get_calls() const { return m_calls; }
    /*
     * If this Variant object has N valid VariantCall objects, then create
     * N variants each with a single valid VariantCall
     */
    void move_calls_to_separate_variants(const VariantQueryConfig& query_config,
        std::vector<Variant>& variants, std::vector<uint64_t>& query_row_idx_in_order,
        GA4GHCallInfoToVariantIdx& call_info_2_variant, GA4GHPagingInfo* paging_info=0);
    /*Non-const iterators for iterating over valid calls*/
    valid_calls_iterator begin() { return valid_calls_iterator(m_calls.begin(), m_calls.end(), 0ull); }
    valid_calls_iterator end() { return valid_calls_iterator(m_calls.end(), m_calls.end(), m_calls.size()); }
    /*const iterators for iterating over valid calls*/
    const_valid_calls_iterator begin() const { return const_valid_calls_iterator(m_calls.begin(), m_calls.end(), 0ull); }
    const_valid_calls_iterator end() const { return const_valid_calls_iterator(m_calls.end(), m_calls.end(), m_calls.size()); }
    /*
     * Set call field idx call_field_idx for call call_idx
     */
    void set_call_field(uint64_t call_idx, unsigned call_field_idx, std::unique_ptr<VariantFieldBase>& field)
    {
      assert(call_idx < m_calls.size());
      m_calls[call_idx].set_field(call_field_idx, field);
    }
    /*
     * Get call field idx call_field_idx for call call_idx
     */
    inline std::unique_ptr<VariantFieldBase>& get_call_field(uint64_t call_idx, unsigned call_field_idx)
    {
      assert(call_idx < m_calls.size());
      return m_calls[call_idx].get_field(call_field_idx);
    }
    /* Return query config */
    const VariantQueryConfig* get_query_config() const { return m_query_config; }
    /** print **/
    void print(std::ostream& stream=std::cout, const VariantQueryConfig* query_config=0, const std::string& indent_prefix="") const;
    /*
     * Binary serialize into buffer
     */
    void binary_serialize(std::vector<uint8_t>& buffer, uint64_t& offset) const;
    /*
     * Deserialize header of Variant (column interval, num calls etc)
     */
    void binary_deserialize_header(const std::vector<uint8_t>& buffer, uint64_t& offset, unsigned num_queried_attributes);
    /*
     * Deep copy variant from src to this. Avoid using as much as possible
     */
    void copy_from_variant(const Variant& src);
    /*
     * Functions dealing with fields in this object, common across all Calls
     */
    void resize_common_fields(size_t size)
    {
      m_fields.resize(size);
      m_common_fields_query_idxs.resize(size);
    }
    /*
     * Move unique ptr
     */
    void set_common_field(unsigned idx, unsigned query_idx, std::unique_ptr<VariantFieldBase>& field)
    {
      assert(idx < m_fields.size());
      m_common_fields_query_idxs[idx] = query_idx;
      m_fields[idx] = std::move(field); //transfer ownership of pointer
    }
    /*
     * Create unique_ptr around raw ptr
     */
    void set_common_field(unsigned idx, unsigned query_idx, VariantFieldBase* field)
    {
      assert(idx < m_fields.size());
      m_common_fields_query_idxs[idx] = query_idx;
      m_fields[idx] = std::move(std::unique_ptr<VariantFieldBase>(field)); //transfer ownership of pointer
    }
    /*
     * Return common field ptr
     */
    const std::unique_ptr<VariantFieldBase>& get_common_field(unsigned idx) const
    {
      assert(idx < m_fields.size());
      return m_fields[idx];
    }
    std::unique_ptr<VariantFieldBase>& get_common_field(unsigned idx)
    {
      assert(idx < m_fields.size());
      return m_fields[idx];
    }
    /*
     * Get query idx given idx in common field vector
     */
    unsigned get_query_idx_for_common_field(unsigned idx) const
    {
      assert(idx < m_common_fields_query_idxs.size());
      return m_common_fields_query_idxs[idx];
    }
    /*
     * Set query idx for given idx in common field vector
     */
    void set_query_idx_for_common_field(unsigned idx, unsigned query_idx)
    {
      assert(idx < m_common_fields_query_idxs.size());
      m_common_fields_query_idxs[idx] = query_idx;
    }
    std::vector<unsigned>& get_common_fields_query_idxs() { return m_common_fields_query_idxs; }
    const std::vector<unsigned>& get_common_fields_query_idxs() const { return m_common_fields_query_idxs; }
    unsigned get_num_common_fields() const { return m_fields.size(); }
    const std::vector<std::unique_ptr<VariantFieldBase>>& get_common_fields() const { return m_fields; }
    /*
     * Copies simple member elements from other as well as simple member
     * elements of other's VariantCall objects
     */
    void deep_copy_simple_members(const Variant& other);
    /*
     * Copy common fields and query idxs associated with those fields
     */
    void copy_common_fields(const Variant& other);
  private:
    //Function that moves from other to self
    void move_in(Variant& other); 
    /*
     * Copies simple member elements from other
     */
    void copy_simple_members(const Variant& other);
    /*
     * Member data elements - check clear, copy, move_in functions while adding new members
     */
    std::vector<VariantCall> m_calls;
    const VariantQueryConfig* m_query_config;
    uint64_t m_col_begin;
    uint64_t m_col_end;
    /*
     * Fields common across VariantCalls for this interval
     */
    std::vector<std::unique_ptr<VariantFieldBase>> m_fields;
    /*
     * Query idx of common fields
     */
    std::vector<unsigned> m_common_fields_query_idxs;
};

//Comparator for sorting VariantCall idxs in column major order
class VariantCallIdxColumnMajorLT
{
  public:
    VariantCallIdxColumnMajorLT(const Variant* variant) : m_variant(variant) { }
    bool operator() (const int64_t i, const int64_t j) const
    {
      auto& call_i = m_variant->get_call(i);
      auto& call_j = m_variant->get_call(j);
      return (call_i.get_column_begin() < call_j.get_column_begin() ||
          (call_i.get_column_begin() == call_j.get_column_begin() && call_i.get_row_idx() < call_j.get_row_idx()));
    }
  private:
    const Variant* m_variant;
};

//GA4GH page token exception
class InvalidGA4GHPageTokenException : public std::exception {
  public:
    InvalidGA4GHPageTokenException(const std::string m="Invalid GA4GH page token exception") : msg_(m) { ; }
    ~InvalidGA4GHPageTokenException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};


/*
 * Class to store information about paging in a GA4GH query
 */
class GA4GHPagingInfo
{
  public:
    GA4GHPagingInfo()
    {
      m_max_num_variants_per_page = 1000000000u;        //large number, 1B
      reset();
    }
    void reset()
    {
      m_is_query_completed = false;
      m_last_row_idx = 0;
      m_last_column_idx = 0;
      m_last_page_end_token = "";
      m_num_handled_variants_in_last_column = 0u;
      m_num_variants_to_shift_left = 0u;
      m_num_variants_in_curr_page = 0u;
    }
    //Call at the start of new page
    void init_page_query();
    //Page limit functions
    inline void set_page_size(unsigned page_size) { m_max_num_variants_per_page = page_size; }
    inline unsigned get_page_size() const { return m_max_num_variants_per_page; }
    /*
     * Page limit is hit if one of two conditions satisfied:
     * (a) #new variants inserted in this page >= page_size. This condition is important when the page limit
     * has not been hit previously and the variable m_num_variants_in_curr_page is un-initialized
     * (b) m_num_variants_in_curr_page is initialized and has hit the page limit. The variable is initialized
     * after set_last_cell_info() function is called. Hence, all subsequent calls to is_page_limit_hit()
     * will return true (until init_page_query() function is called, which resets the value of m_num_variants_in_curr_page)
     */
    inline bool is_page_limit_hit(unsigned variants_vec_size) const
    { return ( (variants_vec_size >= get_num_handled_variants_in_last_column() + get_page_size())
        || (m_num_variants_in_curr_page >= get_page_size()) ); }
    inline unsigned get_num_variants_in_curr_page() const { return m_num_variants_in_curr_page; }
    void set_last_cell_info(std::vector<Variant>& variants,
        const uint64_t row_idx, const uint64_t column_idx, const unsigned num_last_column_variants_handled_after_curr_page);
    void shift_left_variants(std::vector<Variant>& variants);
    inline bool handled_previously(const uint64_t row_idx, const uint64_t column_idx)
    {
      //Column major order storage
      //hence, if m_last_column_idx is less than column_idx, not handled
      if(column_idx > m_last_column_idx)
        return false;
      else
        if(column_idx < m_last_column_idx)      //if m_last_column_idx is greater than column_idx, handled
          return true;
        else    //m_last_column_idx == column_idx
        {
          return false; //need to check all cells with the same column_idx, because of how GA4GH variants are constructed
#if 0
          //sorted by rows
          if(row_idx > m_last_row_idx)
            return false;
          else
            return true;
#endif
        }
    }
    inline uint64_t get_last_column() const { return m_last_column_idx; }
    inline unsigned get_num_handled_variants_in_last_column()  const { return m_num_handled_variants_in_last_column; }
    /*
     * End of query, no more pages
     */
    inline bool is_query_completed() const { return m_is_query_completed; }
    inline void set_query_completed(std::vector<Variant>& variants)
    {
      m_is_query_completed = true;
      set_last_cell_info(variants, ULLONG_MAX, ULLONG_MAX, 0u);
    }
    //GA4GH string tokens
    void serialize_page_end(const std::string& array_name);
    void deserialize_page_end();
    const std::string& get_page_end_token() const { return m_last_page_end_token; }
    void set_page_end_token(const std::string& last_end_token) { m_last_page_end_token = last_end_token; }
  private:
    bool m_is_query_completed;
    unsigned m_max_num_variants_per_page;      //page size
    uint64_t m_last_row_idx;       //Set by query to mark end of curr page, i.e., next query will
    uint64_t m_last_column_idx;    //access cells with row,column > these values
    unsigned m_num_handled_variants_in_last_column;     //#variants starting at the last column that were returned in the last page
    unsigned m_num_variants_to_shift_left;      //helps determine how many elements need to be removed from the head of the vector before returning
    unsigned m_num_variants_in_curr_page;
    //GA4GH string tokens
    std::string m_last_page_end_token;
};

#define PAGE_END_CHECK_LOGIC \
if(paging_info && newly_inserted) \
{ \
  /*Moved to a different column, num_last_column_variants_handled_after_curr_page gets reset*/ \
  if(last_column_idx != curr_column_idx) \
    num_last_column_variants_handled_after_curr_page = 1u; \
  else \
    ++num_last_column_variants_handled_after_curr_page; \
  if(paging_info->is_page_limit_hit(variants.size())) /*Page limit hit*/ \
  { \
    /*First time page limit hit is detected*/ \
    if(paging_info->get_num_variants_in_curr_page() == 0u) \
    { \
      stop_inserting_new_variants = true; \
      paging_info->set_last_cell_info(variants, curr_row_idx, curr_column_idx, num_last_column_variants_handled_after_curr_page); \
    } \
    else /*page limit was hit in a previous iteration, have moved to a new column, break out of loop*/ \
         /*since no more Calls can be added to any Variant in the variants vector */ \
      if(paging_info->get_last_column() != curr_column_idx) \
        break; \
  } \
}

//Priority queue ordered by END position of intervals for VariantCall objects
//Ensures that interval with the smallest end is at the top of the PQ/min-heap
struct EndCmpVariantCallStruct
{
  bool operator()(const VariantCall* x, const VariantCall* y) { return x->get_column_end() > y->get_column_end(); }
};
typedef std::priority_queue<VariantCall*, std::vector<VariantCall*>, EndCmpVariantCallStruct> VariantCallEndPQ;

/*
 * Function that checks whether a ptr is NULL or not
 * The template parameter do_assert controls whether to do actual check or not
 * By default, do nothing
 */
template<bool do_assert>
inline void assert_not_null(const void* ptr)
{ }

/*Specialization, actually checks*/
/*Only in DEBUG compile mode*/
template<>
inline void assert_not_null<true>(const void* ptr)
{
  assert(ptr);
}

/*
 * Get VariantFieldTy* for given known field enum, if the query requested it
 */
template<class VariantFieldTy, bool do_assert>
VariantFieldTy* get_known_field_if_queried(VariantCall& curr_call, const VariantQueryConfig& query_config,
    unsigned known_field_enum)
{
  if(query_config.is_defined_query_idx_for_known_field_enum(known_field_enum))
  {
    auto* field_ptr = curr_call.get_field<VariantFieldTy>
      (query_config.get_query_idx_for_known_field_enum(known_field_enum));
    assert_not_null<do_assert>(static_cast<void*>(field_ptr));
    return field_ptr;
  }
  else
    return nullptr;
}
/*
 * Get VariantFieldTy* for given known field enum, if the query requested it (no checks)
 */
template<class VariantFieldTy, bool do_assert>
VariantFieldTy* get_known_field(VariantCall& curr_call, const VariantQueryConfig& query_config,
    unsigned known_field_enum)
{
  auto* field_ptr = curr_call.get_field<VariantFieldTy>
    (query_config.get_query_idx_for_known_field_enum(known_field_enum));
  assert_not_null<do_assert>(static_cast<void*>(field_ptr));
  return field_ptr;
}
/*
 * Get VariantFieldTy* for given known field enum, if the query requested it (no checks)
 */
template<class VariantFieldTy, bool do_assert>
const VariantFieldTy* get_known_field(const VariantCall& curr_call, const VariantQueryConfig& query_config,
    unsigned known_field_enum)
{
  auto* field_ptr = curr_call.get_field<VariantFieldTy>
    (query_config.get_query_idx_for_known_field_enum(known_field_enum));
  assert_not_null<do_assert>(static_cast<const void*>(field_ptr));
  return field_ptr;
}

/*
 * Move call to variants vector - create new Variant if necessary
 */
bool move_call_to_variant_vector(const VariantQueryConfig& query_config, VariantCall& to_move_call,
    std::vector<Variant>& variants, GA4GHCallInfoToVariantIdx& call_info_2_variant, bool stop_inserting_new_variants);

enum VariantOutputFormatEnum
{
  COTTON_JSON_OUTPUT_FORMAT_IDX=0,
  POSITIONS_JSON_OUTPUT_FORMAT_IDX,
  GA4GH_OUTPUT_FORMAT_IDX,
  DEFAULT_OUTPUT_FORMAT_IDX
};

enum VaritantPrintTypesEnum
{
    INDICES_IDX = 0,
    START_IDX,
    END_IDX,
    ATTRIBUTE_IDX
};

/*
 * JSON as required by John and Cotton
 */
void print_Cotton_JSON(std::ostream& fptr, const std::vector<Variant>& variants, const VariantQueryConfig& query_config, const VidMapper* id_mapper);
/*
 * Prints variants in requested format
 */
void print_variants(const std::vector<Variant>& variants,
                    const std::string& output_format,
                    const VariantQueryConfig& query_config,
                    std::ostream& fptr=std::cout,
                    const bool is_partitioned_by_column = true,
                    const VidMapper* id_mapper = 0,
                    const std::vector<uint64_t>& query_column_lengths = std::vector<uint64_t>(),
                    const std::vector<uint64_t>& num_column_intervals = std::vector<uint64_t>(),
                    const std::vector<uint64_t>& queried_column_positions = std::vector<uint64_t>(),
                    bool output_directly=false);

/*
 * Copies field from src to dst. Optimized to reduce #re-allocations
 * Handles the case where src and/or dst may be null
 */
void copy_field(std::unique_ptr<VariantFieldBase>& dst, const std::unique_ptr<VariantFieldBase>& src);
void copy_fields(std::vector<std::unique_ptr<VariantFieldBase>>& dst, const std::vector<std::unique_ptr<VariantFieldBase>>& src);
#endif
