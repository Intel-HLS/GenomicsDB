/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Intel Corporation
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

#ifndef GENOMICSDB_MULTID_VECTOR_FIELD_H
#define GENOMICSDB_MULTID_VECTOR_FIELD_H

#include "headers.h"
#include "variant_field_data.h"

class GenomicsDBMultiDVectorFieldOperator;
/*
 * Class that holds the binary encoded multiD vector field
 * The internal representation is stored in TileDB and is used for compute operations
 * However, having a class hides the implementation details from clients
 *
 * Example: int F[X][Y] is stored as
 * <size_of_data starting at F> <F[0] vector data> <F[1] vector data> <#elements in F> <offset of F[0], offset of F[1], ...> ...
 * Storing the offsets after the data allows us to determine #elements and offset values and then append them to
 * the buffer.
 *
 * int F[X][Y][Z] is stored as
 * <size_of_data at F>
 *   <size_of_data at F[0]><F[0][0] data><F[0][1] data>.. <#elements in F[0]><offset of F[0][0],F[0][1]..>
 *    ....
 * <#elements in F><offset of F[0]><offset of F[1]>...
 *
 * All sizes and offsets are uint64_t
 *
 */
class GenomicsDBMultiDVectorField
{
  public:
    GenomicsDBMultiDVectorField(const FieldInfo& field_info,
        const uint8_t* ro_data=0, const size_t ro_size=0ull)
    {
      m_field_info_ptr = &field_info;
      m_ro_field_ptr = ro_data;
      m_ro_data_size = ro_size;
      m_rw_data_size = 0ull;
    }
    //Destructor
    ~GenomicsDBMultiDVectorField()
    {
    }
    //Default copy constructor
    GenomicsDBMultiDVectorField(const GenomicsDBMultiDVectorField& other) = default;
    //Default move constructor
    GenomicsDBMultiDVectorField(GenomicsDBMultiDVectorField&& other) = default;

    const uint8_t* get_ro_data_ptr() const { return m_ro_field_ptr; }
    const FieldInfo* get_field_info() const { return m_field_info_ptr; }
    template<class ElementType>
    void parse_and_store_numeric(const char* str, const size_t str_length);
    /*
     * Traverses the multi-d vector and invokes the operator for innermost vector
     * Arguments to the operator include uint8_t* ptr, size of vector, index vector
     */
    void run_operation(GenomicsDBMultiDVectorFieldOperator& multid_vector_field_operator, const uint8_t* data_ptr) const;
    const std::vector<uint8_t>& get_rw_data() const { return m_rw_field_data; }
  private:
    const uint8_t* m_ro_field_ptr; //pointer to the data buffer - read-only used mostly in querying
    std::vector<uint8_t> m_rw_field_data; //useful when constructing a new object - example: deserializing a VCF field encoded as string 
    const FieldInfo* m_field_info_ptr;
    size_t m_ro_data_size;
    size_t m_rw_data_size;
};

//How do you access a multi-D vector?
//Ideally something similar to MultiDVec[i][j][k]
//The first invocation of MultiDVec returns an object of type GenomicsDBMultiDVectorIdx
//Subsequent invocations of the operator modify the returned object to add additional indexes
//
//Indexing could have been achieved by a function similar to:
//index(std::vector<size_t> idx_vec)
//However, this would create many heap operations (std::vector) which I wished to avoid. With
//the GenomicsDBMultiDVectorIdx class, the object is on the stack in most use cases
class GenomicsDBMultiDVectorIdx
{
  public:
    GenomicsDBMultiDVectorIdx(const uint8_t* data_ptr, const FieldInfo* field_info_ptr)
    {
      assert(data_ptr);
      assert(field_info_ptr);
      m_current_dim_idx = -1;
      m_current_index_in_current_dimension = 0u;
      m_num_entries_in_current_dim = 1;
      m_ro_field_ptr = data_ptr;
      m_offsets_ptr = 0;
      m_field_info_ptr = field_info_ptr;
    }
    GenomicsDBMultiDVectorIdx(const uint8_t* data_ptr, const FieldInfo* field_info_ptr, const size_t idx)
      : GenomicsDBMultiDVectorIdx(data_ptr, field_info_ptr)
    {
      advance_to_index_in_next_dimension(idx);
    }
    /*
     * Advances pointer to the index idx in the current dimension
     * Then advances m_current_dim_idx
     */
    void advance_to_index_in_next_dimension(const size_t idx);

    /*
     * Advances index in same dimension
     */
    void advance_index_in_current_dimension();

    /*
     * Should be invoked for the next to last dimension or the last dimension
     * Will provide a pointer to a sequence of elements of type T
     * Example for 2D array int F[X][Y], idx(F,5).get_ptr<int>() will return
     * a pointer to &(F[5][0]) while idx(F,5).advance_index_in_current_dimension(4).get_ptr() will
     * return a pointer to &(F[5][4])
     */
    template<class T>
    inline const T* get_ptr() const
    {
      //Must be the next to lowest dim idx
      assert(static_cast<size_t>(m_current_dim_idx)+2u == m_field_info_ptr->m_length_descriptor.get_num_dimensions());
      return reinterpret_cast<const T*>(m_ro_field_ptr);
    }
    template<class T>
    inline T get_element() const
    {
      //Must be the lowest dim 
      assert(static_cast<size_t>(m_current_dim_idx)+1u == m_field_info_ptr->m_length_descriptor.get_num_dimensions());
      return *(reinterpret_cast<const T*>(m_ro_field_ptr));
    }
    /*
     * Returns current dim idx
     */
    inline int get_current_dim_index() const { return m_current_dim_idx; }
    /*
     * Returns current index in current dimension
     */
    inline uint64_t get_current_index_in_current_dimension() const { return m_current_index_in_current_dimension; }
    /*
     * #bytes in current dimension - will include offsets, total size and the data
     * Except for the innermost dimension - this will be the size of data
     */
    inline uint64_t get_size_of_current_index() const
    {
      assert(static_cast<size_t>(m_current_dim_idx) <
          m_field_info_ptr->m_length_descriptor.get_num_dimensions());
      //the innermost dimension is simply a vector of elements - the size, #elements and offsets are NOT stored on disk,
      //Since the innermost dimension is simply a raw vector, hence the +2u check for dimensions which store 
      //size, #elements, offsets on disk
      if(m_current_dim_idx+1u < m_field_info_ptr->m_length_descriptor.get_num_dimensions())
      {
        assert(m_current_index_in_current_dimension < m_num_entries_in_current_dim);
        return m_offsets_ptr[m_current_index_in_current_dimension+1u]
          - m_offsets_ptr[m_current_index_in_current_dimension];
      }
      else
        return m_field_info_ptr->get_element_size();
    }
    inline uint64_t get_num_entries_in_current_dimension() const { return m_num_entries_in_current_dim; }
  private:
    int m_current_dim_idx;
    uint64_t m_current_index_in_current_dimension;
    uint64_t m_num_entries_in_current_dim;
    const uint8_t* m_ro_field_ptr;
    const uint64_t* m_offsets_ptr;
    const FieldInfo* m_field_info_ptr;
};

class GenomicsDBMultiDVectorFieldOperator
{
  public:
    virtual void operate(const uint8_t* ptr, const size_t size_of_data,
        const std::vector<uint64_t>& idx_vector,
        int outermost_dim_idx_changed_since_last_call_to_operate) = 0;
};

class GenomicsDBMultiDVectorFieldOperatorException : public std::exception {
  public:
    GenomicsDBMultiDVectorFieldOperatorException(const std::string m="Unhandled") : msg_(m) { ; }
    ~GenomicsDBMultiDVectorFieldOperatorException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};


class GenomicsDBMultiDVectorFieldVCFPrinter : public GenomicsDBMultiDVectorFieldOperator
{
  public:
    GenomicsDBMultiDVectorFieldVCFPrinter(std::ostream& fptr, const FieldInfo& field_info);
    void operate(const uint8_t* ptr, const size_t size_of_data,
        const std::vector<uint64_t>& idx_vector,
        int outermost_dim_idx_changed_since_last_call_to_operate);
  private:
    bool m_first_call;
    VariantFieldTypeEnum m_type_enum;
    std::ostream* m_fptr;
    const FieldInfo* m_field_info_ptr;
};

#endif
