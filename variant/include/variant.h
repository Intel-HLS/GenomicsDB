#ifndef VARIANT_H
#define VARIANT_H

#include <assert.h>
#include "gt_common.h"
#include "variant_field_data.h"
#include "variant_query_config.h"

/*a string to store <NON_REF> string (read-only)*/
extern std::string g_non_reference_allele;

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
     * Set field transfers ownership of field data to member unique ptr. The argument field is released from
     * managing more memory
     */
    inline void set_field(unsigned idx, std::unique_ptr<VariantFieldBase>& field)
    {
      assert(idx < m_fields.size());
      assert(m_fields[idx].get() == 0);        //should not be managing any memory
      VariantFieldBase* ptr = field.release();      //Release field from management
      m_fields[idx] = std::move(std::unique_ptr<VariantFieldBase>(ptr)); //transfer ownership of pointer
    }
    /*
     * Set field through raw pointer
     */
    inline void set_field(unsigned idx, VariantFieldBase* field)
    {
      assert(idx < m_fields.size());
      assert(m_fields[idx].get() == 0);        //should not be managing any memory
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
    void print(std::ostream& stream, const VariantQueryConfig* query_config=0) const;
    /**
     * Deep copy VariantCall, avoid using as much as possible (performance)
     */
    void copy_from_call(const VariantCall& src);
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
     * Member data elements - check clear, copy, move_in functions while adding new members
     */
    //Could be initialized, but invalid (no data for this column interval)
    bool m_is_valid;
    //If false, not initialized (not yet considered in query)
    bool m_is_initialized;
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
    void clear()
    {
      for(auto& call : m_calls)
        call.clear();
      m_calls.clear();
    }
    /*
     * Same query_config, but new interval is starting. Reset what needs to be reset
     */
    void reset_for_new_interval();
    /* Set query config ptr */
    void set_query_config(const VariantQueryConfig* query_config) { m_query_config = query_config; }
    /**
     * Allocate maximum possible #calls and #fields per call based on m_query_config
     */
    void resize_based_on_query()
    {
      assert(m_query_config);
      assert(m_query_config->is_bookkeeping_done());
      //Initialize VariantCall vector and pointer vector
      uint64_t num_rows = m_query_config->get_num_rows_to_query();
      resize(num_rows, m_query_config->get_num_queried_attributes());
      for(uint64_t i=0ull;i<num_rows;++i)
      {
        uint64_t row_idx = m_query_config->get_array_row_idx_for_query_row_idx(i);
        m_calls[i].set_row_idx(row_idx);
      }
    }
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
    void move_calls_to_separate_variants(std::vector<Variant>& variants, std::vector<uint64_t>& query_row_idx_in_order);
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
    void print(std::ostream& stream=std::cout, const VariantQueryConfig* query_config=0) const;
    /*
     * Deep copy variant from src to this. Avoid using as much as possible
     */
    void copy_from_variant(const Variant& src);
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
};

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
#endif
