/**
 * The MIT License (MIT)
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

#ifndef GENOMICSDB_COLUMNAR_FIELD_H
#define GENOMICSDB_COLUMNAR_FIELD_H

#include "headers.h"
#include "vcf.h"
#include "gt_common.h"

class GenomicsDBColumnarFieldException : public std::exception {
  public:
    GenomicsDBColumnarFieldException(const std::string m) : msg_(m) { ; }
    ~GenomicsDBColumnarFieldException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

//ceil(buffer_size/field_size)*field_size
#define GET_ALIGNED_BUFFER_SIZE(buffer_size, field_size) ((((buffer_size)+(field_size)-1u)/(field_size))*(field_size))

class GenomicsDBBuffer
{
  public:
    GenomicsDBBuffer(const size_t num_bytes)
    {
      m_is_in_live_list = false;
      m_buffer.resize(num_bytes);
      m_valid.resize(num_bytes);
      m_num_live_entries = 0ull;
      m_num_filled_entries = 0ull;
      m_num_unprocessed_entries = 0ull;
      m_next_buffer = 0;
      m_previous_buffer = 0;
    }
    //Delete both move and copy constructors
    GenomicsDBBuffer(const GenomicsDBBuffer& other) = delete;
    GenomicsDBBuffer(GenomicsDBBuffer&& other) = delete;
    //Functions
    /*
     * Set if this buffer object is in the live list
     */
    inline void set_is_in_live_list(const bool val) { m_is_in_live_list = val; }
    /*
     * Check if this buffer object is in the live list
     */
    inline bool is_in_live_list() const { return m_is_in_live_list; }
    //Live and filled entry control
    /*
     * Set number of live entries in this buffer
     */
    inline size_t get_num_live_entries() const { return m_num_live_entries; }
    /*
     * Get number of live entries in this buffer
     */
    inline void set_num_live_entries(const size_t n)
    {
      m_num_live_entries = n;
      m_num_filled_entries = n;
      m_num_unprocessed_entries = n;
    }
    /*
     * Decrement number of live entries in this buffer by value
     */
    inline void decrement_num_live_entries(const size_t value)
    {
      assert(m_num_live_entries >= value);
      m_num_live_entries -= value;
    }
    /*
     * Increment number of live entries in this buffer by value
     */
    inline void increment_num_live_entries(const size_t value=1u)
    {
      m_num_live_entries += value;
    }
    /*
     * Get number of filled entries in this buffer
     */
    inline size_t get_num_filled_entries() const
    {
      return m_num_filled_entries;
    }
    /*
     * Get number of unprocessed entries in this buffer
     */
    inline size_t get_num_unprocessed_entries() const
    {
      return m_num_unprocessed_entries;
    }
    /*
     * Decrement number of unprocessed entries in this buffer by value
     */
    inline void decrement_num_unprocessed_entries(const size_t value)
    {
      assert(m_num_unprocessed_entries >= value);
      m_num_unprocessed_entries -= value;
    }
    /*
     * Get raw byte buffer in read-write mode
     */
    inline std::vector<uint8_t>& get_buffer() { return m_buffer; }
    /*
     * Get raw byte buffer in read-only mode
     */
    inline const std::vector<uint8_t>& get_buffer() const { return m_buffer; }
    /*
     * Get pointer into raw byte buffer in read-write mode
     */
    inline uint8_t* get_buffer_pointer() { return &(m_buffer[0]); }
    /*
     * Get pointer into raw byte buffer in read-only mode
     */
    inline const uint8_t* get_buffer_pointer() const { return &(m_buffer[0]); }
    /*
     * Get raw byte buffer size in bytes
     */
    inline size_t get_buffer_size_in_bytes() const { return m_buffer.size(); }
    /*inline std::vector<size_t>& get_offsets() { return m_offsets; }*/
    /*
     * Get pointer into offsets vector in read-write mode
     */
    inline size_t* get_offsets_pointer() { return &(m_offsets[0]); }
    /*
     * Get number of bytes of offset elements that can be written when fetching from TileDB
     */
    inline size_t get_offsets_size_in_bytes() const { return ((m_offsets.size()-1u)*sizeof(size_t)); } //last offsets element holds the "size"
    /*
     * Get number of offset elements that can be written when fetching from TileDB
     */
    inline size_t get_offsets_length() const { return (m_offsets.size()-1u); } //last offsets element holds the "size"
    /*
     * Set offset in buffer of the idx-th element
     * Useful only for variable length fields and used by the iterators
     * This allows us to use m_offsets for computing the size of every cell data without
     * having to write "if(last_cell) else" expressions
     */
    inline void set_offset(const size_t idx, const size_t value)
    {
      assert(idx < m_offsets.size());
      m_offsets[idx] = value;
    }
    /*
     * Get offset of idx-th element
     */
    inline size_t get_offset(const size_t idx) const
    {
      assert(idx < m_offsets.size() && idx < m_num_filled_entries);
      return m_offsets[idx ];
    }
    /*
     * Get size in bytes of the idx-th element
     */
    inline size_t get_size_of_variable_length_field(const size_t idx) const
    {
      assert(idx+1u < m_offsets.size() && idx < m_num_filled_entries);
      return m_offsets[idx+1u] - m_offsets[idx];
    }
    /*
     * Return true if the idx-th cell in the buffer is valid (!= TILEDB_NULL)
     */
    inline bool is_valid(const size_t idx) const
    {
      assert(idx < m_valid.size() && idx < m_num_filled_entries);
      return m_valid[idx];
    }
    /*
     * Return bitmask vector specifying the validity of each cell of the buffer
     * in read-only mode
     */
    inline const std::vector<bool>& get_valid_vector() const { return m_valid; }
    /*
     * Return bitmask vector specifying the validity of each element of the buffer
     * in read-write mode
     */
    inline std::vector<bool>& get_valid_vector() { return m_valid; }
    //Buffer pointers in a doubly linked list
    /*
     * Set pointer to next buffer in the doubly linked list
     */
    inline void set_next_buffer(GenomicsDBBuffer* next) { m_next_buffer = next; }
    /*
     * Set pointer to previous buffer in the doubly linked list
     */
    inline void set_previous_buffer(GenomicsDBBuffer* previous) { m_previous_buffer = previous; }
    /*
     * Get pointer to next buffer in the doubly linked list in read-write mode
     */
    inline GenomicsDBBuffer* get_next_buffer() { return m_next_buffer; }
    /*
     * Get pointer to previous buffer in the doubly linked list in read-write mode
     */
    inline GenomicsDBBuffer* get_previous_buffer() { return m_previous_buffer; }
  protected:
    //bool used only in assertions - for debugging
    bool m_is_in_live_list;
    //Byte buffer for storing raw data from TileDB etc
    std::vector<uint8_t> m_buffer;
    //Bitmask specifying validity for each cell read from TileDB etc
    std::vector<bool> m_valid;
    //For variable length fields, stores the offsets
    //If the buffer has data from N-1 cells, this vector has N elements
    //The total size of the variable length field is always stored in the Nth element
    //This allows us to use m_offsets for computing the size of every cell data without
    //having to write "if(last_cell) else" expressions
    std::vector<size_t> m_offsets;
    //m_num_filled_entries - the #items filled in the buffer from TileDB
    //m_num_unprocessed_entries - a counter that's decremented by iterators
    //m_num_unprocessed_entries <= m_num_filled_entries
    size_t m_num_live_entries;
    size_t m_num_filled_entries;
    size_t m_num_unprocessed_entries;
    //Pointers in the linked list of buffers
    GenomicsDBBuffer* m_next_buffer;
    GenomicsDBBuffer* m_previous_buffer;
};

typedef GenomicsDBBuffer GenomicsDBFixedSizeFieldBuffer;

class GenomicsDBVariableSizeFieldBuffer: public GenomicsDBBuffer
{
  public:
    GenomicsDBVariableSizeFieldBuffer(const size_t num_bytes)
      : GenomicsDBBuffer(num_bytes)
    {
      m_offsets.resize(GET_ALIGNED_BUFFER_SIZE(num_bytes, sizeof(size_t))/sizeof(size_t), 0ull);
    }
};

//Printer
//Partial class specialization is allowed, but partial function template specialization isn't
template<typename T, bool print_as_list>
class GenomicsDBColumnarFieldPrintOperator
{
  public:
    static void print(std::ostream& fptr, const uint8_t* ptr, const size_t num_elements);
    static void print_csv(std::ostream& fptr, const uint8_t* ptr, const size_t num_elements,
        const bool is_variable_length_field, const bool is_valid);
};

/*
 * An STL RandomAccessIterator for iterating over a single coords buffer
 * Useful in STL algorithms like upper_bound etc
 */
class GenomicsDBCoordinatesIteratorInBuffer
: std::iterator<std::random_access_iterator_tag, const int64_t*, size_t, const int64_t*, const int64_t*>

{
  public:
    GenomicsDBCoordinatesIteratorInBuffer(const GenomicsDBBuffer* coords_buffer, const size_t cell_idx)
    {
      m_buffer = reinterpret_cast<const int64_t*>(coords_buffer->get_buffer_pointer());
      m_num_filled_entries = coords_buffer->get_num_filled_entries();
      m_cell_idx = cell_idx;
    }
    //default copy, move constructors and equality operators

    //pre-increment
    GenomicsDBCoordinatesIteratorInBuffer& operator++()
    {
      ++m_cell_idx;
      return *this;
    }
    //pre-decrement
    GenomicsDBCoordinatesIteratorInBuffer& operator--()
    {
      --m_cell_idx;
      return *this;
    }
    GenomicsDBCoordinatesIteratorInBuffer& operator+=(const size_t value)
    {
      m_cell_idx += value;
      return *this;
    }
    GenomicsDBCoordinatesIteratorInBuffer& operator-=(const size_t value)
    {
      m_cell_idx -= value;
      return *this;
    }
    const int64_t* operator*()
    {
      assert(m_cell_idx < m_num_filled_entries);
      return m_buffer + (m_cell_idx << 1); //each cell has two int64 coordinates
    }
  private:
    const int64_t* m_buffer;
    size_t m_num_filled_entries;
    size_t m_cell_idx;
};

/*
 * Class that stores data for a given field from multiple VariantCalls/TileDB cells together in a columnar fashion
 * Maintains a big buffer for data, valid bitvector and offsets array
 */
class GenomicsDBColumnarField
{
  public:
    GenomicsDBColumnarField(const std::type_index element_type, const int length_descriptor,
        const unsigned fixed_length_field_length, const size_t num_bytes);
    ~GenomicsDBColumnarField();
    //Delete copy constructor
    GenomicsDBColumnarField(const GenomicsDBColumnarField& other) = delete;
    //Define move constructor
    GenomicsDBColumnarField(GenomicsDBColumnarField&& other);
    void clear() { }
    //Static functions - function pointers hold addresses to one of these functions
    //Check validity
    template<typename T>
    static bool check_tiledb_valid_element(const uint8_t* ptr, const size_t num_elements)
    {
      auto data = reinterpret_cast<const T*>(ptr);
      for(auto i=0ull;i<num_elements;++i)
        if(!is_tiledb_missing_value<T>(data[i]))
          return true;
      return false;
    }
    //Buffer pointer handling
    GenomicsDBBuffer* get_or_allocate_free_buffer()
    {
      if(m_free_buffer_list_head_ptr == 0)
        add_new_buffer();
      return m_free_buffer_list_head_ptr;
    }
    GenomicsDBBuffer* get_live_buffer_list_head_ptr() { return m_live_buffer_list_head_ptr; }
    GenomicsDBBuffer* get_live_buffer_list_tail_ptr() { return m_live_buffer_list_tail_ptr; }
    const GenomicsDBBuffer* get_live_buffer_list_head_ptr() const { return m_live_buffer_list_head_ptr; }
    const GenomicsDBBuffer* get_live_buffer_list_tail_ptr() const { return m_live_buffer_list_tail_ptr; }
    void move_buffer_to_live_list(GenomicsDBBuffer* buffer);
    void move_buffer_to_free_list(GenomicsDBBuffer* buffer);
    void move_all_buffers_from_live_list_to_free_list();
    void set_valid_vector_in_live_buffer_list_tail_ptr();
    //Field lengths
    bool is_variable_length_field() const { return (m_length_descriptor != BCF_VL_FIXED); }
    inline size_t get_fixed_length_field_size_in_bytes() const { return m_fixed_length_field_size; }
    //Functions to update index in live list tail buffer
    inline void set_curr_index_in_live_list_tail(const size_t val) { m_curr_index_in_live_buffer_list_tail = val; }
    inline void advance_curr_index_in_live_list_tail(const size_t val=1u) { m_curr_index_in_live_buffer_list_tail += val; }
    inline size_t get_curr_index_in_live_list_tail() const { return m_curr_index_in_live_buffer_list_tail; }
    //Get pointer to data in GenomicsDBBuffer*
    inline const uint8_t* get_pointer_to_data_in_buffer_at_index(const GenomicsDBBuffer* buffer_ptr,
        const size_t index) const
    {
      assert(buffer_ptr);
      if(m_length_descriptor == BCF_VL_FIXED)
        return buffer_ptr->get_buffer_pointer() + (m_fixed_length_field_size*index);
      else
        return buffer_ptr->get_buffer_pointer() + buffer_ptr->get_offset(index);
    }
    inline uint8_t* get_pointer_to_data_in_buffer_at_index(GenomicsDBBuffer* buffer_ptr,
        const size_t index)
    {
      assert(buffer_ptr);
      if(m_length_descriptor == BCF_VL_FIXED)
        return buffer_ptr->get_buffer_pointer() + (m_fixed_length_field_size*index);
      else
        return buffer_ptr->get_buffer_pointer() + buffer_ptr->get_offset(index);
    }
    //Get num bytes for current idx in buffer
    inline size_t get_size_of_data_in_buffer_at_index(const GenomicsDBBuffer* buffer_ptr,
        const size_t index) const
    {
      assert(buffer_ptr);
      if(m_length_descriptor == BCF_VL_FIXED)
        return m_fixed_length_field_size;
      else
        return buffer_ptr->get_size_of_variable_length_field(index);
    }
    inline size_t get_length_of_data_in_buffer_at_index(const GenomicsDBBuffer* buffer_ptr,
        const size_t index) const
    {
      return get_size_of_data_in_buffer_at_index(buffer_ptr, index) >> m_log2_element_size;
    }
    inline bool is_valid_data_in_buffer_at_index(const GenomicsDBBuffer* buffer_ptr,
        const size_t index) const
    {
      return buffer_ptr->is_valid(index);
    }
    /*
     * Prints data in buffer as list (variable length non-string fields)
     * or element (single element or string fields) to fptr
     */
    void print_data_in_buffer_at_index(std::ostream& fptr,
        const GenomicsDBBuffer* buffer_ptr, const size_t index) const;
    /*
     * Print ALT alleles as a list of strings
     */
    void print_ALT_data_in_buffer_at_index(std::ostream& fptr,
        const GenomicsDBBuffer* buffer_ptr, const size_t index) const;
    /*
     * Print data in a CSV format that can be imported into GenomicsDB directly
     */
    void print_data_in_buffer_at_index_as_csv(std::ostream& fptr,
        const GenomicsDBBuffer* buffer_ptr, const size_t index) const;
    //Get list lengths
    size_t get_free_buffer_list_length() const { return m_free_buffer_list_length; }
    size_t get_live_buffer_list_length() const { return m_live_buffer_list_length; }
  private:
    void copy_simple_members(const GenomicsDBColumnarField& other);
    void assign_function_pointers();
    template<bool print_as_list>
    void assign_print_function_pointers(const int bcf_ht_type);
    void add_new_buffer()
    {
      auto buffer_ptr = create_new_buffer();
      if(m_free_buffer_list_head_ptr)
      {
        assert(m_free_buffer_list_head_ptr->get_previous_buffer() == 0);
        m_free_buffer_list_head_ptr->set_previous_buffer(buffer_ptr);
        buffer_ptr->set_next_buffer(m_free_buffer_list_head_ptr);
      }
      m_free_buffer_list_head_ptr = buffer_ptr;
      ++m_free_buffer_list_length;
    }
    GenomicsDBBuffer* create_new_buffer() const
    {
      return (m_length_descriptor == BCF_VL_FIXED) ? (new GenomicsDBFixedSizeFieldBuffer(m_buffer_size))
        : (new GenomicsDBVariableSizeFieldBuffer(m_buffer_size));
    }
  private:
    int m_length_descriptor;
    unsigned m_fixed_length_field_num_elements;
    unsigned m_fixed_length_field_size;
    unsigned m_element_size;
    unsigned m_log2_element_size;
    std::type_index m_element_type;
    //Function pointers - assigned based on type of data
    //Function pointer that determines validity check
    bool (*m_check_tiledb_valid_element)(const uint8_t* ptr, const size_t num_elements);
    void (*m_print)(std::ostream& fptr, const uint8_t* ptr, const size_t num_elements);
    void (*m_print_csv)(std::ostream& fptr, const uint8_t* ptr, const size_t num_elements,
        const bool is_variable_length_field, const bool is_valid);
    size_t m_buffer_size;
    //Head of list containing free buffers - can be used in invocation of tiledb_array_read()
    GenomicsDBBuffer* m_free_buffer_list_head_ptr;
    //Head of list containing buffers with live data - list is in column-major order
    GenomicsDBBuffer* m_live_buffer_list_head_ptr;
    //Tail of list containing buffers with live data - list is in column-major order
    //tail needed because new data gets appended to end
    GenomicsDBBuffer* m_live_buffer_list_tail_ptr;
    //Index of the element in the live buffer list head that must be read next
    size_t m_curr_index_in_live_buffer_list_tail;
    //Length of live and free lists
    size_t m_live_buffer_list_length;
    size_t m_free_buffer_list_length;
};

#endif
