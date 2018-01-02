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

#ifndef VID_MAPPER_HD
#define VID_MAPPER_HD

#include "headers.h"
#include "vcf.h"
#include "variant_array_schema.h"
#include "known_field_info.h"
#include "rapidjson/document.h"

typedef std::tuple<std::string, int64_t, int64_t> ContigIntervalTuple;

inline bool contig_offset_idx_pair_cmp(const std::pair<int64_t, int>& first, const std::pair<int64_t, int>& second)
{
  return (first.first < second.first);
}

/*
 * Info object regarding callset
 */
class CallSetInfo
{
  public:
    CallSetInfo()
    {
      m_is_initialized = false;
      m_row_idx = -1;
      m_file_idx = -1;
      m_idx_in_file = 0;
    }
    void set_info(const int64_t row_idx, const std::string& name, const int64_t file_idx=-1, const int64_t idx_in_file=0)
    {
      m_is_initialized = true;
      m_row_idx = row_idx;
      m_file_idx = file_idx;
      m_name = name;
      m_idx_in_file = idx_in_file;
    }
    bool m_is_initialized;
    int64_t m_row_idx;
    int64_t m_file_idx;
    int64_t m_idx_in_file;
    std::string m_name;
};

/*
 * Info object regarding contig
 */
class ContigInfo
{
  public:
    ContigInfo()
    {
      m_contig_idx = -1;
      m_tiledb_column_offset = -1;
    }
    void set_info(const int contig_idx, const std::string& name, const int64_t length, const int64_t offset)
    {
      m_contig_idx = contig_idx;
      m_name = name;
      m_length = length;
      m_tiledb_column_offset = offset;
    }
    int m_contig_idx;
    int64_t m_length;
    int64_t m_tiledb_column_offset;
    std::string m_name;
};

enum VidFileTypeEnum
{
  VCF_FILE_TYPE=0,
  SORTED_CSV_FILE_TYPE,
  UNSORTED_CSV_FILE_TYPE,
  VCF_BUFFER_STREAM_TYPE,
  BCF_BUFFER_STREAM_TYPE
};

class FileInfo
{
  public:
    FileInfo()
    {
      m_file_idx = -1;
      m_owner_idx = -1;
      m_local_file_idx = -1;
      m_local_tiledb_row_idx_pairs.clear();
      m_type = VidFileTypeEnum::VCF_FILE_TYPE;
      //Buffer stream details
      m_buffer_stream_idx = -1;
      m_buffer_capacity = 1024; //1KB
      m_initialization_buffer_num_valid_bytes = 0u;
      //Split files info
      m_single_split_file_path = false;
    }
    void set_info(const int64_t file_idx, const std::string& name)
    {
      m_file_idx = file_idx;
      m_local_file_idx = file_idx;
      m_name = name;
    }
    bool add_local_tiledb_row_idx_pair(const int64_t local, const int64_t global, int64_t& other_row_idx);
    size_t get_num_callsets() const { return m_local_tiledb_row_idx_pairs.size(); }
    size_t get_num_orders() const;
    std::string m_name;
    int64_t m_file_idx;
    //Idx of the entity that handles this file (used when loaders and owners are distinct MPI processes)
    int m_owner_idx;
    //Local idx in the owner's file idx vector list (m_owner_idx_to_file_idx_vec[owner][m_local_file_idx] == m_file_idx)
    int64_t m_local_file_idx;
    //A VCF can contain multiple callsets - each entry in this map contains a vector of pair<local_idx,row_idx>,
    //corresponding to each callset in the VCF
    std::vector<std::pair<int64_t, int64_t>> m_local_tiledb_row_idx_pairs;
    std::unordered_map<int64_t, int64_t> m_local_idx_to_tiledb_row_idx;
    //File/stream type enum
    VidFileTypeEnum m_type;
    //Buffer stream details
    int64_t m_buffer_stream_idx;
    size_t m_buffer_capacity;
    std::vector<uint8_t> m_initialization_buffer;
    size_t m_initialization_buffer_num_valid_bytes;
    //Split files output locations
    bool m_single_split_file_path;
    std::vector<std::string> m_split_files_paths;
};

typedef FileInfo BufferStreamInfo;

enum VCFFieldCombineOperationEnum
{
  VCF_FIELD_COMBINE_OPERATION_SUM=0,
  VCF_FIELD_COMBINE_OPERATION_MEAN,
  VCF_FIELD_COMBINE_OPERATION_MEDIAN,
  VCF_FIELD_COMBINE_OPERATION_DP,      //for the DP INFO field
  VCF_FIELD_COMBINE_OPERATION_MOVE_TO_FORMAT,
  VCF_FIELD_COMBINE_OPERATION_ELEMENT_WISE_SUM,
  VCF_FIELD_COMBINE_OPERATION_CONCATENATE,
  VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION
};

class FieldLengthDescriptorComponent
{
  public:
    FieldLengthDescriptorComponent()
    {
      m_num_elements = 1;
      m_length_descriptor = BCF_VL_FIXED;
    }
    int m_num_elements;
    int m_length_descriptor;
};

class FieldLengthDescriptor
{
  public:
    FieldLengthDescriptor()
    {
      m_length_descriptor_vec.emplace_back();
      m_num_elements = 1u;
      m_is_fixed_length_field = true;
      m_is_length_allele_dependent = false;
      m_is_length_all_alleles_dependent = false;
      m_is_length_genotype_dependent = false;
      m_is_length_ploidy_dependent = false;
    }
    //MultiD vectors
    void resize(const size_t n)
    {
      m_length_descriptor_vec.resize(n);
      m_vcf_delimiter_vec.resize(n);
    }
    size_t get_num_dimensions() const { return m_length_descriptor_vec.size(); }
    //Fixed length components
    void set_num_elements(const int length_dim_idx, const size_t n)
    {
      assert(static_cast<size_t>(length_dim_idx) < m_length_descriptor_vec.size());
      m_length_descriptor_vec[length_dim_idx].m_num_elements = n;
      m_num_elements *= n;
    }
    size_t get_num_elements() const
    {
      assert(is_fixed_length_field());
      return m_num_elements;
    }
    //length descriptor
    void set_length_descriptor(const int length_dim_idx, const int length_descriptor);
    unsigned get_length_descriptor(const int length_dim_idx) const
    {
      assert(static_cast<size_t>(length_dim_idx) < m_length_descriptor_vec.size());
      return m_length_descriptor_vec[length_dim_idx].m_length_descriptor;
    }
    bool is_fixed_size_dimension(const size_t idx) const
    {
      assert(idx < m_length_descriptor_vec.size());
      return (m_length_descriptor_vec[idx].m_length_descriptor == BCF_VL_FIXED);
    }
    size_t get_num_elements_in_dimension(const size_t idx)
    {
      assert(is_fixed_size_dimension(idx));
      return m_length_descriptor_vec[idx].m_num_elements;
    }
    //is fixed length field
    bool is_fixed_length_field() const { return m_is_fixed_length_field; }
    //Allele dependency
    bool is_length_allele_dependent() const { return m_is_length_allele_dependent; }
    bool is_length_only_ALT_alleles_dependent() const
    { return m_is_length_allele_dependent && !m_is_length_all_alleles_dependent; }
    bool is_length_genotype_dependent() const { return m_is_length_genotype_dependent; }
    //Ploidy
    bool is_length_ploidy_dependent() const { return m_is_length_ploidy_dependent; }
    unsigned get_ploidy_step_value() const
    {
      assert(m_length_descriptor_vec.size() == 1u);
      assert(m_is_length_ploidy_dependent);
      return (get_length_descriptor(0u) == BCF_VL_Phased_Ploidy) ? 2u : 1u;
    }
    unsigned get_ploidy(const unsigned n) const
    {
      assert(m_length_descriptor_vec.size() == 1u);
      assert(m_is_length_ploidy_dependent);
      return KnownFieldInfo::get_ploidy(get_length_descriptor(0u), n);
    }
    bool contains_phase_information() const
    {
      assert(m_length_descriptor_vec.size() == 1u);
      assert(m_is_length_ploidy_dependent);
      return (get_length_descriptor(0u) == BCF_VL_Phased_Ploidy);
    }
    size_t get_num_elements(const unsigned num_ALT_alleles, const unsigned ploidy, const unsigned num_elements);
    void set_vcf_delimiter(const size_t dim_idx, const char* vcf_delim)
    {
      assert(vcf_delim);
      assert(dim_idx < m_vcf_delimiter_vec.size());
      m_vcf_delimiter_vec[dim_idx] = vcf_delim[0];
    }
    char get_vcf_delimiter(const size_t dim_idx) const
    {
      assert(dim_idx < m_vcf_delimiter_vec.size());
      return m_vcf_delimiter_vec[dim_idx];
    }
  private:
    //Length descriptors - could be multi-dimensional array: example [ "R", 4 ] 2D vector
    std::vector<FieldLengthDescriptorComponent> m_length_descriptor_vec;
    //Useful only for fixed length fields
    size_t m_num_elements;
    //Flags for fast querying - summarize information in length desc components
    //Is fixed length field
    bool m_is_fixed_length_field;
    //Allele, genotype, ploidy dependency
    bool m_is_length_allele_dependent;
    bool m_is_length_all_alleles_dependent;
    bool m_is_length_genotype_dependent;
    bool m_is_length_ploidy_dependent;
    //VCF delimiter
    std::vector<char> m_vcf_delimiter_vec;
};

/*
 * To take into account tuples
 * Element type
 */
class FieldElementTypeDescriptor
{
  public:
    //Constructors
    FieldElementTypeDescriptor(const unsigned num_entries_in_tuple);
    FieldElementTypeDescriptor(const std::type_index& curr_type, const int ht_type);
    //#elements in tuple
    void resize_num_elements_in_tuple(const unsigned num_entries_in_tuple);
    size_t get_num_elements_in_tuple() const { return m_tuple_element_type_vec.size(); }
    int get_tuple_element_bcf_ht_type(const unsigned idx) const
    {
      assert(idx < m_tuple_element_bcf_ht_type_vec.size());
      return m_tuple_element_bcf_ht_type_vec[idx];
    }
    const std::type_index& get_tuple_element_type_index(const unsigned idx) const
    {
      assert(idx < m_tuple_element_type_vec.size());
      return m_tuple_element_type_vec[idx];
    }
    size_t get_tuple_element_size(const unsigned idx) const
    {
      assert(idx < m_tuple_element_size_vec.size());
      return m_tuple_element_size_vec[idx];
    }
    void set_tuple_element_type(const unsigned idx, const std::type_index& curr_type, const int ht_type);
  private:
    std::vector<std::type_index> m_tuple_element_type_vec;
    std::vector<int> m_tuple_element_bcf_ht_type_vec;
    std::vector<size_t> m_tuple_element_size_vec;
};

class FieldInfo
{
  friend class VidMapper;
  public:
    FieldInfo()
      : m_tiledb_type(1u),
        m_genomicsdb_type(1u),
        m_vcf_type(1u),
        m_length_descriptor()
    {
      m_is_vcf_FILTER_field = false;
      m_is_vcf_INFO_field = false;
      m_is_vcf_FORMAT_field = false;
      m_is_flattened_field = false;
      m_field_idx = -1;
      m_VCF_field_combine_operation = VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION;
      m_element_index_in_tuple = 0u;
    }
    void set_info(const std::string& name, int idx)
    {
      m_name = name;
      m_vcf_name = name;
      m_field_idx = idx;
    }
    //Type information
    /*
     * By default, TileDB, GenomicsDB and VCF types are the same
     */
    void set_type(const FieldElementTypeDescriptor& type)
    {
      m_tiledb_type = type;
      m_genomicsdb_type = type;
      m_vcf_type = type;
      compute_element_size();
    }
    /*
     * Set GenomicsDB type
     */
    void set_genomicsdb_type(const FieldElementTypeDescriptor& type)
    {
      m_genomicsdb_type = type;
      compute_element_size();
    }
    /*
     * Set VCF type
     */
    void set_vcf_type(const FieldElementTypeDescriptor& type)
    {
      m_vcf_type = type;
    }
    /*
     * Set TileDB type
     */
    void set_tiledb_type(const FieldElementTypeDescriptor& type)
    {
      m_tiledb_type = type;
    }
    const FieldElementTypeDescriptor& get_tiledb_type() const { return m_tiledb_type; }
    const FieldElementTypeDescriptor& get_vcf_type() const { return m_vcf_type; }
    const FieldElementTypeDescriptor& get_genomicsdb_type() const { return m_genomicsdb_type; }
    size_t get_element_size() const { return m_element_size; }
    //Composite fields get flattened - in the FieldInfo objects for the
    //flattened field, this contains the index in the composite tuple
    void set_element_index_in_tuple(const unsigned idx)
    {
      m_element_index_in_tuple = idx;
    }
    unsigned get_element_index_in_tuple() const { return m_element_index_in_tuple; }
    void set_is_flattened_field(const bool val)
    {
      m_is_flattened_field = val;
    }
    bool is_flattened_field() const { return m_is_flattened_field; }
    /*
     * For flattened fields, sets the index of the parent composite field
     */
    void set_parent_composite_field_idx(const unsigned idx)
    {
      m_parent_composite_field_idx = idx;
    }
    unsigned get_parent_composite_field_idx() const
    {
      assert(is_flattened_field());
      return m_parent_composite_field_idx;
    }
    //Public members
    std::string m_name;     //Unique per array schema
    std::string m_vcf_name; //VCF naming mess - DP could be FORMAT and INFO - in this case m_name=DP_FORMAT, m_vcf_name = DP
    bool m_is_vcf_FILTER_field;
    bool m_is_vcf_INFO_field;
    bool m_is_vcf_FORMAT_field;
    int m_field_idx;
    //Length descriptor
    FieldLengthDescriptor m_length_descriptor;
    //Combine operation for VCF INFO fields
    int m_VCF_field_combine_operation;
    //Multi-d vector fields - different types in TileDB/VCF/GenomicsDB
    void modify_field_type_if_multi_dim_field();
    void compute_element_size();
  private:
    //Type info
    //TileDB type
    FieldElementTypeDescriptor m_tiledb_type;
    //GenomicsDB type index - could be different from TileDB and VCF
    FieldElementTypeDescriptor m_genomicsdb_type;
    //VCF type info - could be different from TileDB type
    FieldElementTypeDescriptor m_vcf_type;
    //Element size - computed from components of tuple
    size_t m_element_size;
    //Composite fields get flattened - in the FieldInfo objects for the
    //flattened field, this contains the index in the composite tuple
    unsigned m_element_index_in_tuple;
    bool m_is_flattened_field;
    unsigned m_parent_composite_field_idx;
};

/*
 * Base class for mapping callset/contig names to rows/columns
 * Many implementations possible (PostgreSQL, file, SQLite etc)
 */
class VidMapper
{
  public:
    VidMapper()
    {
      clear();
      m_is_initialized = false;
      m_is_callset_mapping_initialized = false;
      m_max_callset_row_idx = -1;
    }
    void clear();
    inline bool is_initialized() const { return m_is_initialized; }
    inline bool is_callset_mapping_initialized() const { return m_is_callset_mapping_initialized; }
    /*
     * Given a position in a flattened 'address' space [TileDB column idx], get the contig_name and location
     * in the contig [0-based]
     * Returns true if valid contig found, false otherwise
     */
    bool get_contig_location(const int64_t position, std::string& contig_name, int64_t& contig_position) const;
    /*
     * Given a position in a flattened 'address' space [TileDB column idx], get the next contig_name and starting
     * location of the contig in the flattened space 
     * Returns true if valid contig found, false otherwise
     */
    bool get_next_contig_location(int64_t query_position, std::string& next_contig_name, int64_t& next_contig_offset) const;
    /*
     * Given a 0-based position in a contig, return the position in the flattened 'address' space [TileDB  column idx]
     * Returns true if valid position found, false otherwise
     */
    bool get_tiledb_position(int64_t& position, const std::string& contig_name, const int64_t contig_position=0ll) const;
    /*
     * Given a TileDB row idx in array, obtain callset name
     * Returns true if valid callset name found, false otherwise
     */
    bool get_callset_name(const int64_t row_idx, std::string& callset_name) const;
    /*
     * Given a callset name, return TileDB row idx in array
     * Returns true if valid callset name found, false otherwise
     */
    bool get_tiledb_row_idx(int64_t& row_idx, const std::string& callset_name) const;
    uint64_t get_num_callsets() const { return m_row_idx_to_info.size(); }
    /*
     * Return callset info
     */
    inline const CallSetInfo& get_callset_info(const int64_t tiledb_row_idx) const
    {
      assert(static_cast<size_t>(tiledb_row_idx) < m_row_idx_to_info.size());
      return m_row_idx_to_info[tiledb_row_idx];
    }
    /*
     * Given a TileDB row idx, return the file_idx 
     */
    bool get_global_file_idx_for_row(int64_t row_idx, int64_t& file_idx) const
    {
      if(row_idx >= static_cast<int64_t>(m_row_idx_to_info.size()))
        return false;
      file_idx = m_row_idx_to_info[row_idx].m_file_idx;
      return (file_idx < 0) ? false : true;
    }
    /*
     * Returns local file idx
     */
    bool get_local_file_idx_for_row(int64_t row_idx, int64_t& file_idx) const
    {
      auto found = get_global_file_idx_for_row(row_idx, file_idx);
      if(!found)
        return false;
      file_idx = m_file_idx_to_info[file_idx].m_local_file_idx;
      return (file_idx < 0) ? false : true;
    }
    /*
     * Given a row idx, returns the owner that is responsible for producing the data
     * local_file_idx is the index in the owner's file idx vector m_owner_idx_to_file_idx_vec
     */
    inline bool get_owner_idx_for_row_idx(const int64_t row_idx, int& owner_idx, int& local_file_idx)
    {
      int64_t file_idx;
      auto found_file = get_global_file_idx_for_row(row_idx, file_idx);
      if(!found_file) return false;
      assert(static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
      owner_idx = m_file_idx_to_info[file_idx].m_owner_idx;
      local_file_idx = m_file_idx_to_info[file_idx].m_local_file_idx;
      return (owner_idx < 0 || local_file_idx < 0) ? false : true;
    }
    /*
     * Total #files
     */
    int64_t get_num_files() const { return m_file_idx_to_info.size(); }
    /*
     * Given a filename, return local-global idx pairs for callsets
     */
    bool get_local_tiledb_row_idx_vec(const std::string& filename, std::vector<int64_t>& row_idx_vec) const
    {
      auto iter = m_filename_to_idx.find(filename);
      if(iter != m_filename_to_idx.end())
      {
        auto file_idx = (*iter).second;
        assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
        return get_local_tiledb_row_idx_vec(file_idx, row_idx_vec);
      }
      else
      {
        for(auto i=0u;i<row_idx_vec.size();++i)
          row_idx_vec[i] = -1;
        return false;
      }
    }
    /*
     * Given a file idx, return local-global idx pairs for callsets
     */
    bool get_local_tiledb_row_idx_vec(const int64_t file_idx, std::vector<int64_t>& row_idx_vec) const
    {
      assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
      for(const auto& pair : m_file_idx_to_info[file_idx].m_local_tiledb_row_idx_pairs)
      {
        if(static_cast<size_t>(pair.first) >= row_idx_vec.size())
          row_idx_vec.resize(static_cast<size_t>(pair.first)+1ull, -1ll);
        row_idx_vec[pair.first] = pair.second;
      }
      return true;
    }
    /*
     * Get idx in file for a given row
     */
    unsigned get_idx_in_file_for_row_idx(const int64_t row_idx) const
    {
      assert(row_idx >= 0 && static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
      return m_row_idx_to_info[row_idx].m_idx_in_file;
    }
    /*
     * Get global file idx for filename, if not exist append and return last index
     */
    int64_t get_or_append_global_file_idx(const std::string& filename)
    {
      auto iter = m_filename_to_idx.find(filename);
      if(iter == m_filename_to_idx.end())
      {
        auto file_idx = m_file_idx_to_info.size();
        iter = m_filename_to_idx.insert(std::make_pair(filename, file_idx)).first;
        m_file_idx_to_info.emplace_back();
        m_file_idx_to_info[file_idx].set_info(file_idx, filename);
      }
      return (*iter).second;
    }
    /*
     * Given a filename, return #callsets within that file being processed 
     */
    bool get_num_callsets_in_file(const std::string& filename, int& num_callsets) const
    {
      auto iter = m_filename_to_idx.find(filename);
      if(iter == m_filename_to_idx.end())
      {
        num_callsets = -1;
        return false;
      }
      else
      {
        auto file_idx = (*iter).second;
        assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
        num_callsets = m_file_idx_to_info[file_idx].m_local_tiledb_row_idx_pairs.size();
        return true;
      }
    }
    /*
     * Return file_idx given filename
     */
    bool get_global_file_idx(const std::string& filename, int64_t& file_idx) const
    {
      auto iter = m_filename_to_idx.find(filename);
      if(iter == m_filename_to_idx.end())
        return false;
      file_idx = (*iter).second;
      return true;
    }
    const FileInfo& get_file_info(const int64_t file_idx) const
    {
      assert(static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
      return m_file_idx_to_info[file_idx];
    }
    /*
     * Get type of file
     */
    inline unsigned get_file_type(const int64_t file_idx) const
    {
      assert(static_cast<size_t>(file_idx) < m_file_idx_to_info.size());
      return m_file_idx_to_info[file_idx].m_type;
    }
    inline bool get_file_type(const std::string& filename, unsigned& file_type) const
    {
      int64_t file_idx = -1;
      auto status = get_global_file_idx(filename, file_idx);
      if(!status)
        return false;
      file_type = get_file_type(file_idx);
      return true;
    }
    /*
     * For a given owner and local file idx, return the buffer stream idx
     * If the "file" is not a buffer stream, the return value will be -1
     */
    inline int64_t get_buffer_stream_idx_for_local_file_idx(const int owner_idx, const int64_t local_file_idx) const
    {
      assert(static_cast<size_t>(owner_idx) < m_owner_idx_to_file_idx_vec.size());
      assert(static_cast<size_t>(local_file_idx) < m_owner_idx_to_file_idx_vec[owner_idx].size());
      auto global_file_idx = m_owner_idx_to_file_idx_vec[owner_idx][local_file_idx];
      assert(static_cast<size_t>(global_file_idx) < m_file_idx_to_info.size());
      return m_file_idx_to_info[global_file_idx].m_buffer_stream_idx;
    }
    /*
     * For a given owner and buffer stream idx, return the local file idx
     * The owner parameter is mostly irrelevant since buffer stream idxs are associated with each process
     * and buffer streams of a process aren't visible to any other process
     */
    inline int64_t get_local_file_idx_for_buffer_stream_idx(const int owner_idx, const int64_t buffer_stream_idx) const
    {
      assert(static_cast<size_t>(buffer_stream_idx) < m_buffer_stream_idx_to_global_file_idx.size());
      auto global_file_idx = m_buffer_stream_idx_to_global_file_idx[buffer_stream_idx];
      assert(static_cast<size_t>(global_file_idx) < m_file_idx_to_info.size());
      return m_file_idx_to_info[global_file_idx].m_local_file_idx;
    }
    /*
     * Return list of files owned by owner_idx
     */
    const std::vector<int64_t>& get_global_file_idxs_owned_by(int owner_idx) const
    {
      assert(static_cast<size_t>(owner_idx) < m_owner_idx_to_file_idx_vec.size());
      return m_owner_idx_to_file_idx_vec[owner_idx];
    }
    inline const std::vector<int64_t>& get_buffer_stream_idx_to_global_file_idx_vec() const { return m_buffer_stream_idx_to_global_file_idx; }
    void build_file_partitioning(const int partition_idx, const RowRange row_partition);
    void verify_file_partitioning() const;
    //Set path of split file
    void set_single_split_file_path(const int64_t global_file_idx, const std::string& split_output_filename)
    {
      assert(static_cast<size_t>(global_file_idx) < m_file_idx_to_info.size());
      auto& file_info = m_file_idx_to_info[global_file_idx];
      file_info.m_single_split_file_path = true;
      file_info.m_split_files_paths.resize(1u);
      file_info.m_split_files_paths[0u] = split_output_filename;
    }
    /*
     * While splitting files, get path of output split file
     */
    std::string get_split_file_path(const std::string& original_filename, const std::string& results_directory,
        std::string& output_type, const int rank) const;
    /*
     * Given a contig name, return global contig idx
     */
    inline bool get_global_contig_idx(const std::string& name, int& contig_idx) const
    {
      auto iter = m_contig_name_to_idx.find(name);
      if(iter != m_contig_name_to_idx.end())
      {
        contig_idx = (*iter).second;
        return true;
      }
      else
      {
        contig_idx = -1;
        return false;
      }
    }
    /*
     * Parse a string length descriptor and fill up the structure
     */
    void parse_string_length_descriptor(const char* field_name,
        const char* length_value_str,
        const size_t length_value_str_length,
        FieldLengthDescriptor& length_descriptor, const size_t length_dim_idx);
    /*
     * Return std::type_index and BCF_HT_* value given a type_string
     */
    std::pair<std::type_index, int> get_type_index_and_bcf_ht_type(const char* type_string);
    /*
     * Get #fields in VidMapper
     */
    inline unsigned get_num_fields() const
    {
      return m_field_idx_to_info.size();
    }
    /*
     * Given a field name, return global field idx
     */
    inline bool get_global_field_idx(const std::string& name, int& field_idx) const
    {
      auto iter = m_field_name_to_idx.find(name);
      if(iter != m_field_name_to_idx.end())
      {
        field_idx = (*iter).second;
        return true;
      }
      else
      {
        field_idx = -1;
        return false;
      }
    }
    /*
     * Given a global field idx, return field info
     */
    inline const FieldInfo& get_field_info(int field_idx) const
    {
      assert(static_cast<size_t>(field_idx) < m_field_idx_to_info.size());
      return m_field_idx_to_info[field_idx];
    }
    /*
     * Given a field name, return FieldInfo ptr
     * If field is not found, return 0
     */
    inline const FieldInfo* get_field_info(const std::string& name) const
    {
      int field_idx = -1;
      auto status = get_global_field_idx(name, field_idx);
      if(!status)
        return 0;
      return &(get_field_info(field_idx));
    }
    const FieldInfo* get_flattened_field_info(const FieldInfo* field_info,
        const unsigned tuple_element_index) const;
    /*
     * Stores the fields, classifying them as FILTER, INFO, FORMAT etc
     */
    void build_vcf_fields_vectors(std::vector<std::vector<std::string>>& vcf_fields) const;
    void build_tiledb_array_schema(VariantArraySchema*& array_schema, const std::string array_name,
        const bool row_based_partitioning, const RowRange& row_range, const bool compress_fields,
	const bool no_mandatory_VCF_fields) const;
    /*
     * Get num contigs
     */
    inline unsigned get_num_contigs() const { return m_contig_idx_to_info.size(); }
    /*
     * Given a global contig idx, return contig info
     */
    inline const ContigInfo& get_contig_info(const int contig_idx) const
    {
      assert(contig_idx >= 0 && static_cast<size_t>(contig_idx) < m_contig_idx_to_info.size());
      return m_contig_idx_to_info[contig_idx];
    }
    /*
     * Given a contig name, find the info
     * Return true if valid contig info found, else false
     */
    inline bool get_contig_info(const std::string& contig_name, ContigInfo& info) const
    {
      auto iter = m_contig_name_to_idx.find(contig_name);
      if(iter == m_contig_name_to_idx.end())
        return false;
      info = get_contig_info((*iter).second);
      return true;
    }
    inline int64_t get_max_callset_row_idx() const { return m_max_callset_row_idx; }
    /*
     * Utility function for obtaining contigs given a column partition
     * is_zero_based - return 0-based or 1-based chromosome intervals
     */
    static std::vector<ContigIntervalTuple> get_contig_intervals_for_column_partition(
        const std::string& loader_filename,
        const int rank, const bool is_zero_based);
    std::vector<ContigIntervalTuple> get_contig_intervals_for_column_partition(
        const int64_t column_partition_begin, const int64_t column_partition_end, const bool is_zero_based) const;
  protected:
    void add_mandatory_fields();
    void flatten_field(int& field_idx, const int original_field_idx);
  protected:
    //Is initialized
    bool m_is_initialized;
    //are callsets initialized
    bool m_is_callset_mapping_initialized;
    //callset mappings
    std::unordered_map<std::string, int64_t> m_callset_name_to_row_idx;
    std::vector<CallSetInfo> m_row_idx_to_info;
    //contig mappings
    std::unordered_map<std::string, int> m_contig_name_to_idx;
    std::vector<ContigInfo> m_contig_idx_to_info;
    //sorted vectors of pair<contig begin/end, idx> 
    std::vector<std::pair<int64_t, int>> m_contig_begin_2_idx;
    std::vector<std::pair<int64_t, int>> m_contig_end_2_idx;
    //field mapping
    std::unordered_map<std::string, int> m_field_name_to_idx;
    std::vector<FieldInfo> m_field_idx_to_info;
    //file mappings
    std::unordered_map<std::string, int64_t> m_filename_to_idx;
    std::vector<FileInfo> m_file_idx_to_info;
    //Buffer stream order index to local file idx
    std::vector<int64_t> m_buffer_stream_idx_to_global_file_idx;
    //owner idx to file_idx vector
    std::vector<std::vector<int64_t>> m_owner_idx_to_file_idx_vec;
    void sort_and_assign_local_file_idxs_for_partition(const int owner_idx);
    //Static members
    static std::unordered_map<std::string, int> m_length_descriptor_string_to_int;
    static std::unordered_map<std::string, std::type_index> m_typename_string_to_type_index;
    static std::unordered_map<std::string, int> m_typename_string_to_bcf_ht_type;
    //INFO field combine operation
    static std::unordered_map<std::string, int> m_INFO_field_operation_name_to_enum;
    //Max row idx in callset idx file
    int64_t m_max_callset_row_idx;
};

//Exceptions thrown 
class FileBasedVidMapperException : public std::exception {
  public:
    FileBasedVidMapperException(const std::string m="") : msg_("FileBasedVidMapperException : "+m) { ; }
    ~FileBasedVidMapperException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class VidMapperException : public std::exception {
  public:
    VidMapperException(const std::string m="") : msg_("VidMapperException : "+m) { ; }
    ~VidMapperException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

/*
 * File based mapper - users provide a file with all the mappings
 */
class FileBasedVidMapper : public VidMapper
{
  public:
    FileBasedVidMapper() : VidMapper()
    {
      m_lb_callset_row_idx = 0;
      m_ub_callset_row_idx = INT64_MAX-1;
    }

    FileBasedVidMapper(
        const std::string& filename,
        const std::string& callset_mapping_file="",
        const int64_t lb_callset_row_idx=0,
        const int64_t ub_callset_row_idx=INT64_MAX-1,
        const bool is_callset_mapping_required=true) : VidMapper() {

      std::vector<BufferStreamInfo> empty_vec;
      common_constructor_initialization(
          filename,
          empty_vec,
          callset_mapping_file,
          "",
          lb_callset_row_idx, ub_callset_row_idx,
          is_callset_mapping_required);
    }

    FileBasedVidMapper(
        const std::string& filename,
        const std::vector<BufferStreamInfo>& buffer_stream_info_vec,
        const std::string& callset_mapping_file="",
        const std::string& buffer_stream_callset_mapping_json_string="",
        const int64_t lb_callset_row_idx=0,
        const int64_t ub_callset_row_idx=INT64_MAX-1,
        const bool is_callset_mapping_required=true) : VidMapper() {

      common_constructor_initialization(
          filename,
          buffer_stream_info_vec,
          callset_mapping_file,
          buffer_stream_callset_mapping_json_string,
          lb_callset_row_idx, ub_callset_row_idx,
          is_callset_mapping_required);
    }

    //Useful when writing partitioned data
    void write_partition_callsets_json_file(
        const std::string& original_callsets_filename,
        const std::string& results_directory,
        const int rank) const;

    void write_partition_loader_json_file(
        const std::string& original_loader_filename,
        const std::string& original_callsets_filename,
        const std::string& results_directory,
        const int num_partition_callset_mapping_files,
        const int rank) const;

  private:
    void common_constructor_initialization(
        const std::string& filename,
        const std::vector<BufferStreamInfo>& buffer_stream_info_vec,
        const std::string& callset_mapping_file="",
        const std::string& buffer_stream_callset_mapping_json_string="",
        const int64_t lb_callset_row_idx=0,
        const int64_t ub_callset_row_idx=INT64_MAX-1,
        const bool is_callset_mapping_required=true);

  private:
    void parse_callsets_json(
        const std::string& filename,
        const std::vector<BufferStreamInfo>& buffer_stream_info_vec,
        const bool is_file);
    void parse_length_descriptor(const char* field_name,
        const rapidjson::Value& length_json_value,
        FieldLengthDescriptor& length_descriptor, const size_t length_dim_idx);

    void parse_type_descriptor(FieldInfo& field_info, const rapidjson::Value& field_info_json_dict);

    int64_t m_lb_callset_row_idx;
    int64_t m_ub_callset_row_idx;
};

#endif  // VID_MAPPER_HD
