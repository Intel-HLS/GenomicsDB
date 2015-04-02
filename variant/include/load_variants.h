#ifndef LOAD_VARIANTS_H
#define LOAD_VARIANTS_H

#include "gt_common.h"
#include "loader.h"

class VariantLoader : public Loader
{
  public:
    // CONSTRUCTORS AND DESTRUCTORS
    /** 
     * Simple constructor. The workspace is where the loader will create its
     * data. The storage manager is the module the loader interefaces with.
     */
    VariantLoader(const std::string& workspace, StorageManager& storage_manager)
      : Loader(workspace, storage_manager)
    { ; }

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
   void load_CSV_gVCF(const std::string& filename, const char* array_name, const uint64_t max_sample_idx, const bool is_input_sorted,
       const std::string tmp_dir="") const; 
  private:
     /** 
      * Treats the CSV line as a logical cell encompassing coordinates and
      * attribute values, and appends the coordinates to a coordinate tile
      * and the attribute values to the respective attribute tiles. It follows
      * the gVCF format (see Loader::load_CSV_gVCF). 
      */
     void append_cell_gVCF(const ArraySchema& array_schema, 
         CSVLine* csv_line, Tile** tiles) const;
     /** 
      * Creates an array schema conforming to the gVCF info 
      * (see Loader::load_CSV_gVCF). 
      */
     ArraySchema* create_gVCF_array_schema(const char* array_name, const uint64_t max_sample_idx) const;
     /**  
      * Creates the (irregular) tiles of the genome array and sends them to the 
      * storage manager. 
      */
     void make_tiles_irregular_CSV_gVCF(const std::string& filename,
         const StorageManager::ArrayDescriptor* ad,
         const ArraySchema& array_schema) const;
};

#endif
