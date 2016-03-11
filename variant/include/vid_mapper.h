#ifndef VID_MAPPER_HD
#define VID_MAPPER_HD

#include "headers.h"
#include "vcf.h"
#include "variant_array_schema.h"

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
      m_row_idx = -1;
      m_file_idx = -1;
    }
    void set_info(const int64_t row_idx, const std::string& name, const int64_t file_idx=-1)
    {
      m_row_idx = row_idx;
      m_file_idx = file_idx;
      m_name = name;
    }
    int64_t m_row_idx;
    int64_t m_file_idx;
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

class FileInfo
{
  public:
    FileInfo()
    {
      m_file_idx = -1;
      m_owner_idx = -1;
      m_local_file_idx = -1;
      m_local_tiledb_row_idx_pairs.clear();
    }
    void set_info(const int64_t file_idx, const std::string& name)
    {
      m_file_idx = file_idx;
      m_local_file_idx = file_idx;
      m_name = name;
    }
    void add_local_tiledb_row_idx_pair(int local, int64_t global)
    {
      m_local_tiledb_row_idx_pairs.push_back(std::make_pair(local, global));
    }
    size_t get_num_callsets() const { return m_local_tiledb_row_idx_pairs.size(); }
    std::string m_name;
    int64_t m_file_idx;
    //Idx of the entity that handles this file (used when loaders and owners are distinct MPI processes)
    int m_owner_idx;
    //Local idx in the owner's file idx vector list (m_owner_idx_to_file_idx_vec[owner][m_local_file_idx] == m_file_idx)
    int64_t m_local_file_idx;
    //A VCF can contain multiple callsets - each entry in this map contains a vector of pair<local_idx,row_idx>,
    //corresponding to each callset in the VCF
    std::vector<std::pair<int64_t, int64_t>> m_local_tiledb_row_idx_pairs;
};

class FieldInfo
{
  public:
    FieldInfo()
    {
      m_is_vcf_INFO_field = false;
      m_is_vcf_FORMAT_field = false;
      m_field_idx = -1;
      m_type_info = 0;
      m_bcf_ht_type = BCF_HT_VOID;
      m_length_descriptor = BCF_VL_FIXED;
      m_num_elements = 1;
    }
    void set_info(const std::string& name, int idx)
    {
      m_name = name;
      m_field_idx = idx;
    }
    std::string m_name;
    bool m_is_vcf_INFO_field;
    bool m_is_vcf_FORMAT_field;
    int m_field_idx;
    //Type info
    const std::type_info* m_type_info;
    int m_bcf_ht_type;
    //Length descriptors
    int m_length_descriptor;
    int m_num_elements;
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
    }
    void clear();
    inline bool is_initialized() const { return m_is_initialized; }
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
        for(const auto& pair : m_file_idx_to_info[file_idx].m_local_tiledb_row_idx_pairs)
        {
          assert(static_cast<size_t>(pair.first) < row_idx_vec.size());
          row_idx_vec[pair.first] = pair.second;
        }
        return true;
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
        assert(static_cast<size_t>(pair.first) < row_idx_vec.size());
        row_idx_vec[pair.first] = pair.second;
      }
      return true;
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
     * Return list of files owned by owner_idx
     */
    const std::vector<int64_t>& get_global_file_idxs_owned_by(int owner_idx) const
    {
      assert(static_cast<size_t>(owner_idx) < m_owner_idx_to_file_idx_vec.size());
      return m_owner_idx_to_file_idx_vec[owner_idx];
    }
    void build_file_partitioning(const int partition_idx, const RowRange row_partition);
    void verify_file_partitioning() const;
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
    /*
     * Stores the fields, classifying them as FILTER, INFO, FORMAT etc
     */
    void build_vcf_fields_vectors(std::vector<std::vector<std::string>>& vcf_fields) const;
    void build_tiledb_array_schema(VariantArraySchema*& array_schema, const std::string array_name,
        const bool row_based_partitioning, const RowRange& row_range, const bool compress_fields) const;
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
  protected:
    //Is initialized
    bool m_is_initialized;
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
    //owner idx to file_idx vector
    std::vector<std::vector<int64_t>> m_owner_idx_to_file_idx_vec;
    void sort_and_assign_local_file_idxs_for_partition(const int owner_idx);
    //Static members
    static std::unordered_map<std::string, int> m_length_descriptor_string_to_int;
    static std::unordered_map<std::string, const std::type_info*> m_typename_string_to_typeinfo;
    static std::unordered_map<std::string, int> m_typename_string_to_bcf_ht_type;
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
    FileBasedVidMapper() : VidMapper() { ; }
    FileBasedVidMapper(const std::string& filename, const std::string& callset_mapping_file="",
        const int64_t limit_callset_row_idx=INT64_MAX, const bool callsets_file_required=true);
  private:
    void parse_callsets_file(const std::string& filename);
    int64_t m_limit_callset_row_idx;
};

#endif
