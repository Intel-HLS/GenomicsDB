#ifndef VARIANT_QUERY_CONFIG_H
#define VARIANT_QUERY_CONFIG_H

#include "query_config.h"
#include "lut.h"

class VariantQueryConfig : public QueryConfig
{
  public:
    VariantQueryConfig()
      : QueryConfig()
    {
      m_query_idx_known_variant_field_enum_LUT.reset_luts();
      m_done_bookkeeping = false;
    }
    //Query idx <--> known fields mapping
    void resize_LUT(unsigned num_known_fields)
    { m_query_idx_known_variant_field_enum_LUT.resize_luts_if_needed(get_num_queried_attributes(), num_known_fields); }
    //Map queryIdx <--> knownEnumIdx
    inline void add_query_idx_known_field_enum_mapping(unsigned queryIdx, unsigned knownEnumIdx)
    {
      assert(queryIdx < get_num_queried_attributes());
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      m_query_idx_known_variant_field_enum_LUT.add_query_idx_known_field_enum_mapping(queryIdx, knownEnumIdx);
    }
    //Get query idx for given knownEnumIdx
    inline unsigned get_query_idx_for_known_field_enum(unsigned knownEnumIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      return m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx);
    }
    //Check whether query contains given knownEnumIdx
    inline bool is_defined_query_idx_for_known_field_enum(unsigned knownEnumIdx) const
    {
      assert(m_query_idx_known_variant_field_enum_LUT.is_defined_value(knownEnumIdx));
      return m_query_idx_known_variant_field_enum_LUT.is_defined_value(
        m_query_idx_known_variant_field_enum_LUT.get_query_idx_for_known_field_enum(knownEnumIdx));
    }
    inline bool is_bookkeeping_done() const { return m_done_bookkeeping; }
    inline void set_done_bookkeeping(bool value) { m_done_bookkeeping = value; }
  private:
    //Mapping between queried idx and known fields enum
    QueryIdxToKnownVariantFieldsEnumLUT m_query_idx_known_variant_field_enum_LUT;
    bool m_done_bookkeeping;
};

#endif
