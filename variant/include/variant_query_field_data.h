#ifndef VARIANT_QUERY_FIELD_DATA_H
#define VARIANT_QUERY_FIELD_DATA_H

#include "query_field_data.h"

class VariantQueryFieldInfo
{
  public:
    VariantQueryFieldInfo()
    {
      m_NULL_bitidx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_OFFSETS_idx = UNDEFINED_ATTRIBUTE_IDX_VALUE;
      m_length_descriptor = UNDEFINED_ATTRIBUTE_IDX_VALUE;
    }
    unsigned m_NULL_bitidx;
    unsigned m_OFFSETS_idx;
    unsigned m_length_descriptor;
};

class VariantQueryFieldData : public QueryFieldData
{
  public:
    VariantQueryFieldData() 
      : QueryFieldData()
    { ; }
};

#endif
