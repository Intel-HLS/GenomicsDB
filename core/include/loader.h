/**
 * @file   loader.h
 * @author Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2014 Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * This file defines class Loader. It also defines LoaderException, which is 
 * thrown by Loader.
 */

#ifndef LOADER_H
#define LOADER_H

#include "array_schema.h"
#include "csv_file.h"
#include "storage_manager.h"

/** Special value indicating an invalid tile id. */
#define LD_INVALID_TILE_ID std::numeric_limits<uint64_t>::max()

/**
 * The Loader is the module that creates the array layout from raw data 
 * presented in a CSV format. It can load data into arrays with both 
 * regular and irregular tiles, supporting various cell and tile orders.
 */
class Loader {
 public:
  // CONSTRUCTORS AND DESTRUCTORS
  /** 
   * Simple constructor. The workspace is where the loader will create its
   * data. The storage manager is the module the loader interefaces with.
   */
  Loader(const std::string& workspace, StorageManager& storage_manager);
  /** Empty destructor. */
  ~Loader() {}

  // LOADING FUNCTIONS
  /**
   * Loads a CSV file into an array.
   * \param filename The name of the input CSV file.
   * \param array_schema The schema of the array the CSV file is loaded into.
   */
  void load(const std::string& filename, const ArraySchema& array_schema) const;
  /**
   * Loads a CSV-gVCF file into an array. Note that a CSV-gVCF file is a CSV
   * file derived from a gVCF file. Each row in the CSV-gVCF file has the 
   * following values:
   * 
   * SampleID(int64_t),POS(int64_t),END(int64_t),REF(string),ALT_num(char),  
   * ALT_1(string),...,ALT_{ALT_num}(string),QUAL(float),FILTER_num(char),   
   * FILTER_ID_1(int),...,FILTER_ID_{FILTER_num}(int),BaseQRankSum(float),   
   * ClippingRankSum(float),MQRankSum(float),ReadPosRankSum(float),DP(int),  
   * MQ(float),MQ0(int),DP_FMT(int),MIN_DP(int),GQ(int),SB_1(int),SB_2(int), 
   * SB_3(int),SB_4(int),AD_1(int),...,AD_{ALT_num+1}(int),                   
   * PL_1(int),...,PL_{(ALT_num+1)*(ALT_num+2)/2}(int)
   *
   * Explanation of symbols:
   * 
   * 1. SampleID: The ID of the sample.
   *
   * 2. POS: The starting position in the genome. The position is global, i.e.,
   * derived from concatenating all chromosomes and assuming a single position
   * sequence.
   *
   * 3. END: The ending (global) position.
   *
   * 4. REF: The nuceotide string in the reference sample starting at POS.
   *
   * 5. ALT_num: Number of alleles.
   *
   * 6. ALT_i: The i-th allele.
   *
   * 7. QUAL: 
   *
   * 8. FILTER_num: Number of filters.
   *
   * 9. FILTER_ID_i: The i-th filter ID (we assume that there is somewhere
   * a filter ID index/dictionary.
   *
   * 10. BaseQRankSum:
   *
   * 11. ClippingRankSum: 
   *
   * 12. MQRankSum:
   *
   * 13. ReadPosRankSum:
   *
   * 14. DP:
   *
   * 15. MQ:
   *
   * 16. MQ0:
   * 
   * 17. DP_FMT:
   *
   * 18. MIN_DP:
   *
   * 19. GQ:
   *
   * 20. SB_1:
   * 
   * 21. SB_2:
   * 
   * 22. SB_3:
   * 
   * 23. SB_4:
   * 
   * 24. AD_i:
   *
   * 25. PL_i: The i-th likelihood.
   *
   * Note: A missing value in the CSV-gVCF file is represented by '#'.
   *
   *
   * The function will produce a new array with the following schema:
   *
   * Dimensions:
   *
   * 1. SampleID(int64_t)
   *
   * 2. POS(int64_t)  
   *
   * Attributes:
   *
   * 1. END(int64_t)
   *
   * 2. REF(char): Holds the reference string as a set of characters, ending
   * with '\0'. Permissible characters: 'A', 'T', 'C', 'G'.
   *
   * 3. ALT(char): We store the allele strings as sets of characters, ending 
   * with '\0'. Permissible characters: 'A', 'T', 'C', 'G', '&'. The '&'
   * character stands for \<NON_REF\>. After the last allele, '\0' is NOT stored.
   *
   * 4. QUAL(float) 
   *
   * 5. FILTER_ID(int): First we store the number of filters, then the
   * filter IDs. 
   *
   * 6. BaseQRankSum(float)
   *
   * 7. ClippingRankSum(float)
   *
   * 8. MQRankSum(float)
   *
   * 9. ReadPosRankSum(float)
   *
   * 10. DP(int)
   *
   * 11. MQ(float)
   *
   * 12. MQ0(int)
   *
   * 13. DP_FMT(int) 
   *
   * 14. MIN_DP(int)
   *
   * 15. GQ(int)
   *
   * 16. SB_1(int)
   *
   * 17. SB_2(int)
   *
   * 18. SB_3(int)
   *
   * 19. SB_4(int)
   *
   * 20. AD(int): ALT_num+1 values per cell.
   *
   * 21. PL(int): (ALT_num+1)*(ALT_num+2)/2 values per cell.
   *
   * 22. NULL(int): Stores bitmaps where each bit corresponds to an attribute
   * and 1 indicates that this attribute is NULL (i.e., missing). 
   * The bits correspond to: QUAL, BaseQRankSum, ClippingRankSum, 
   * MQRankSum, ReadPosRankSum, DP, MQ, MQ0, DP_FMT, MIN_DP, GQ, SB_1, SB_2,
   * SB_3, SB_4, AD, PL, from LEFT to RIGHT, where PL is the rightmost bit.
   *
   * 23. OFFSETS(int64_t): We store in this order the offset of this cell in 
   * (i)   the REF tile,
   * (ii)  the ALT tile,
   * (iii) the FILTER_ID tile,
   * (iv)  the AD tile,
   * (v)   the PL tile.
   */
  void load_CSV_gVCF(const std::string& filename, const char* array_name, const uint64_t max_sample_idx) const; 

 private:
  // PRIVATE ATTRIBUTES
  /** The StorageManager object the loader interfaces with. */
  StorageManager& storage_manager_;
  /** A folder in the disk where the loader creates all its data. */
  std::string workspace_;
 
  // PRIVATE METHODS
  /** 
   * Treats the CSV line as a logical cell encompassing coordinates and
   * attribute values, and appends the coordinates to a coordinate tile
   * and the attribute values to the respective attribute tiles. 
   */
  void append_cell(const ArraySchema& array_schema, 
                    CSVLine* csv_line, Tile** tiles) const;
  /** 
   * Treats the CSV line as a logical cell encompassing coordinates and
   * attribute values, and appends the coordinates to a coordinate tile
   * and the attribute values to the respective attribute tiles. It follows
   * the gVCF format (see Loader::load_CSV_gVCF). 
   */
  void append_cell_gVCF(const ArraySchema& array_schema, 
                        CSVLine* csv_line, Tile** tiles) const;
  /** Checks upon invoking the load command. */
  bool check_on_load(const std::string& filename) const;
  /** 
   * Creates an array schema conforming to the gVCF info 
   * (see Loader::load_CSV_gVCF). 
   */
  ArraySchema* create_gVCF_array_schema(const char* array_name, const uint64_t max_sample_idx) const;
  /** Creates the workspace folder. */
  void create_workspace() const;
  /**
   * Injects tile/cell ids to the CSV file prior to loading (applies only to
   * regular tiles with any order, and irregular tiles with Hilbert order). 
   */
  void inject_ids_to_csv_file(const std::string& filename, 
                              const std::string& injected_filename,
                              const ArraySchema& array_schema) const;
  /** Creates the (irregular) tiles and sends them to the storage manager. */
  void make_tiles_irregular(const std::string& filename,
                            const StorageManager::ArrayDescriptor* ad,
                            const ArraySchema& array_schema) const;
  /**  
   * Creates the (irregular) tiles of the genome array and sends them to the 
   * storage manager. 
   */
  void make_tiles_irregular_CSV_gVCF(const std::string& filename,
                                     const StorageManager::ArrayDescriptor* ad,
                                     const ArraySchema& array_schema) const;
  /** Creates the (regular) tiles and sends them to the storage manager. */
  void make_tiles_regular(const std::string& filename,
                          const StorageManager::ArrayDescriptor* ad,
                          const ArraySchema& array_schema) const;
  /** 
   * Creates an array of new tile pointers with the input tile id, 
   * and based on the input array info. 
   */
  void new_tiles(const ArraySchema& array_schema, 
                 uint64_t tile_id, Tile** tiles) const;
  /** Returns true if the input path is an existing directory. */
  bool path_exists(const std::string& path) const;
  /** Simply sets the workspace. */
  void set_workspace(const std::string& path);
  /**  Sorts the csv file depending on the type of tiles and order. */
  void sort_csv_file(const std::string& to_be_sorted_filename,
                     const std::string& sorted_filename,
                     const ArraySchema& array_schema) const;
  /** Sends the tiles to the storage manager. */
  void store_tiles(const StorageManager::ArrayDescriptor* ad,
                   Tile** tiles) const;
};

/** This exception is thrown by Loader. */
class LoaderException {
 public:
  // CONSTRUCTORS & DESTRUCTORS
  /** Takes as input the exception message. */
  LoaderException(const std::string& msg) 
      : msg_(msg) {}
  /** Empty destructor. */
  ~LoaderException() {}

  // ACCESSORS
  /** Returns the exception message. */
  const std::string& what() const { return msg_; }

 private:
  /** The exception message. */
  std::string msg_;
};

#endif
