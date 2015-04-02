#ifndef VARIANT_H
#define VARIANT_H

#include "gt_common.h"

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

#endif
