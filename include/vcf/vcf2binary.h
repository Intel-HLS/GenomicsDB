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

#ifndef VCF2BINARY_H
#define VCF2BINARY_H

#ifdef HTSDIR

#include "headers.h"
#include "vid_mapper.h"
#include "htslib/synced_bcf_reader.h"
#include "gt_common.h"
#include "histogram.h"
#include "tiledb_loader_file_base.h" 

//Exceptions thrown 
class VCF2BinaryException : public std::exception {
  public:
    VCF2BinaryException(const std::string m="") : msg_("VCF2BinaryException : "+m) { ; }
    ~VCF2BinaryException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

//Wrapper around VCF's file I/O functions
//Capability of using index only during seek to minimize memory consumption
class VCFReader : public FileReaderBase
{
  public:
    VCFReader();
    virtual ~VCFReader();
    void initialize(const char* filename, const char* regions,
        const std::vector<std::vector<std::string>>& vcf_field_names, const VidMapper* id_mapper, bool open_file);
    //Abstract virtual functions from base class that must be over-ridden
    void add_reader();
    void remove_reader();
    void read_and_advance();
    //Helper functions
    void seek_read_advance(const char* contig, const int pos, bool discard_index);
    bcf_hdr_t* get_header() { return m_hdr; }
    bcf1_t* get_line() { return m_is_line_valid ? m_line : 0; }
  private:
    bcf_srs_t* m_indexed_reader;
    htsFile* m_fptr;
    std::string m_filename;
    bcf_hdr_t* m_hdr;
    bcf1_t* m_line;
    kstring_t m_buffer;
    bool m_is_line_valid;
};

class VCFColumnPartition : public File2TileDBBinaryColumnPartitionBase
{
  friend class VCF2Binary;
  public:
    /*
     * Primary constructor
     */
    VCFColumnPartition()
      : File2TileDBBinaryColumnPartitionBase()
    {
      //buffer for vcf get functions - 16 KB
      m_vcf_get_buffer_size = 16*1024;
      m_vcf_get_buffer = (uint8_t*)malloc(m_vcf_get_buffer_size*sizeof(uint8_t));
    }
    //Delete copy constructor
    VCFColumnPartition(const VCFColumnPartition& other) = delete;
    //Define move constructor
    VCFColumnPartition(VCFColumnPartition&& other);
    ~VCFColumnPartition();
    int64_t get_column_position_in_record() const
    {
      auto vcf_reader_ptr = dynamic_cast<VCFReader*>(m_base_reader_ptr);
      assert(vcf_reader_ptr);
      auto line = vcf_reader_ptr->get_line();
      return (m_contig_tiledb_column_offset + static_cast<int64_t>(line->pos));
    }
  protected:
    //Position in contig from which to fetch next batch of cells
    int m_local_contig_idx;
    int64_t m_contig_position;  //position in contig (0-based)
    int64_t m_contig_tiledb_column_offset;
    //Buffer for obtaining data from htslib 
    uint8_t* m_vcf_get_buffer;
    uint64_t m_vcf_get_buffer_size;
    //Used by move constructor
    void copy_simple_members(const VCFColumnPartition& other);
};

class VCF2Binary : public File2TileDBBinaryBase 
{
  public:
    VCF2Binary(const std::string& vcf_filename, const std::vector<std::vector<std::string>>& vcf_fields,
        unsigned file_idx, VidMapper& vid_mapper, const std::vector<ColumnRange>& partition_bounds,
        size_t max_size_per_callset,
        bool treat_deletions_as_intervals,
        bool parallel_partitions=false, bool noupdates=true, bool close_file=false, bool discard_index=false);
    //Delete default copy constructor as it is incorrect
    VCF2Binary(const VCF2Binary& other) = delete;
    //Define move constructor explicitly
    VCF2Binary(VCF2Binary&& other);
    ~VCF2Binary();
    void clear();
    //Initialization functions
    void initialize(const std::vector<ColumnRange>& partition_bounds);
    void initialize_partition(unsigned idx, const ColumnRange& column_partition);
    //Abstract virtual functions in base class that must be defined 
    bool convert_record_to_binary(std::vector<uint8_t>& buffer, File2TileDBBinaryColumnPartitionBase& partition_info);
    bool seek_and_fetch_position(File2TileDBBinaryColumnPartitionBase& partition_info, bool force_seek, bool advance_reader);
    uint64_t get_num_callsets_in_record(const File2TileDBBinaryColumnPartitionBase& partition_info) const
    { return m_enabled_local_callset_idx_vec.size(); }
    //Helper functions
    void update_local_contig_idx(VCFColumnPartition& vcf_partition, const bcf1_t* line);
    //VCF->TileDB conversion functions
    bool convert_VCF_to_binary_for_callset(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition,
        size_t size_per_callset, uint64_t enabled_callsets_idx);
    /*
     * field_type_idx: BCF_HL_* 
     */
    template<class FieldType>
    bool convert_field_to_tiledb(std::vector<uint8_t>& buffer, VCFColumnPartition& vcf_partition, 
        int64_t& buffer_offset, const int64_t buffer_offset_limit, int local_callset_idx,
        const std::string& field_name, unsigned field_type_idx);
  private:
    bool m_discard_index;
    //Vector of vector of strings, outer vector has 2 elements - 0 for INFO, 1 for FORMAT
    const std::vector<std::vector<std::string>>* m_vcf_fields; 
    std::string m_regions;
    //Local contig idx to global contig idx
    std::vector<int> m_local_contig_idx_to_global_contig_idx;
    //Local field idx to global field idx
    std::vector<int> m_local_field_idx_to_global_field_idx;
    //Used by move constructor
    void copy_simple_members(const VCF2Binary& other);
};

#endif //ifdef HTSDIR

#endif
