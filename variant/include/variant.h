#ifndef VARIANT_H
#define VARIANT_H

#include <assert.h>
#include "gt_common.h"
#include "variant_field_data.h"
#include "variant_query_config.h"

/*a string to store <NON_REF> string (read-only)*/
extern std::string g_non_reference_allele;

/**
 * This class stores and processes information about a single column
 * of the gVCF array. A column corresponds to a unique position in
 * the genome. It holds the REF, ALT and PL values of every row (i.e.,
 * individual) in the array. If a REF value is equal to "$"
 * for a row, then it means that this row is NULL (i.e., it contains 
 * no useful info for joint genotyping). 
 */
class GTColumn {
  public:
    // CONSTRUCTORS AND DESTRUCTORS
    /** Simple constructor. */
    GTColumn(int64_t col, uint64_t row_num);
    ~GTColumn() {}

    // OPERATIONS
    // TODO: Implement function for deriving the final ALT and PL vailues
    // TODO: Implement the genotyping function

    /** Holds the (seqeuence of) ALT values for each row. */
    std::vector<std::vector<std::string> > ALT_;
    /** The id of the column. */
    int64_t col_;
    /** Holds the REF values for each row. */
    std::vector<std::string> REF_;
    /** Holds the (seqeuence of) PL values for each row. */
    std::vector<std::vector<int> > PL_;
    /** Holds the (seqeuence of) AF values for each row. */
    std::vector<float> AF_;
    /** Holds the (seqeuence of) AN values for each row. */
    std::vector<int> AN_;
    /** Holds the (seqeuence of) AC values for each row. */
    std::vector<int> AC_;
    void reset();
};	

//Structure stored as part of a priority queue (min-heap) to align genomic intervals
class PQStruct
{
  public:
    PQStruct()    { m_needs_to_be_processed = false; }
    bool m_needs_to_be_processed;
    int64_t m_end_point;
    int64_t m_sample_idx;
    int64_t m_array_column;
    uint64_t m_cell_pos;
    uint64_t m_tile_idx;
};

//Ensures that interval with smallest end is at the top of the PQ/min-heap
struct CmpPQStruct
{
  bool operator()(const PQStruct* x, const PQStruct* y) { return x->m_end_point > y->m_end_point; }
};

typedef std::priority_queue<PQStruct*, std::vector<PQStruct*>, CmpPQStruct> VariantIntervalPQ;

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
    void add_field(std::unique_ptr<VariantFieldBase>& field)
    {
      m_fields.push_back(std::move(field));
    }
    inline std::vector<std::unique_ptr<VariantFieldBase>>& get_all_fields()
    {
      return m_fields;
    }
    inline std::unique_ptr<VariantFieldBase>& get_field(unsigned idx)
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
    /** print **/
    void print(std::ostream& stream, const VariantQueryConfig* query_config=0) const;
  private:
    /*
     * Performs move from other object
     */
    void move_in(VariantCall& other)
    {
      m_is_valid = other.is_valid();
      m_is_initialized = other.is_initialized();
      m_row_idx = other.get_row_idx();
      clear();
      m_fields.resize(other.get_all_fields().size());
      unsigned idx = 0u;
      for(auto& other_field : other.get_all_fields())
      {
        set_field(idx, other_field);
        ++idx;
      }
    }
    //Could be initialized, but invalid (no data for this column interval)
    bool m_is_valid;
    //If false, not initialized (not yet considered in query)
    bool m_is_initialized;
    uint64_t m_row_idx;
    std::vector<std::unique_ptr<VariantFieldBase>> m_fields;
};

/*
 * Class equivalent to GAVariant in GA4GH API. Stores info about 1 position/interval
 */
class Variant
{
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
    /*
     * Resize call vector
     */
    void resize(uint64_t num_calls, unsigned num_query_call_fields)
    {
      m_calls.resize(num_calls);
      for(uint64_t i=0ull;i<num_calls;++i)
        m_calls[i].resize(num_query_call_fields);
    }
    inline uint64_t get_num_calls() { return m_calls.size(); }
    /**
     * Return VariantCall at index call_idx
     */
    inline VariantCall& get_call(uint64_t call_idx)
    {
      assert(call_idx < m_calls.size());
      return m_calls[call_idx];
    }
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
    /** print **/
    void print(std::ostream& stream) const;
  private:
    //Function that moves from other to self
    void move_in(Variant& other)
    {
      m_query_config = other.m_query_config;
      m_col_begin = other.m_col_begin;
      m_col_end = other.m_col_end;
      //De-allocates existing data
      clear();
      for(auto i=0ull;i<m_calls.size();++i)
        m_calls.emplace_back(std::move(other.get_call(i)));
    }
    std::vector<VariantCall> m_calls;
    const VariantQueryConfig* m_query_config;
    uint64_t m_col_begin;
    uint64_t m_col_end;
};

/*
 * Function that checks whether a ptr is NULL or not
 * The template parameter do_assert controls whether to do actual check or not
 * By default, do nothing
 */
template<bool do_assert>
inline void assert_not_null(void* ptr)
{ }

/*Specialization, actually checks*/
/*Only in DEBUG compile mode*/
template<>
inline void assert_not_null<true>(void* ptr)
{
  assert(ptr);
}

/*
 * Get unique_ptr<VariantFieldTy> for given known field enum
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


#endif
