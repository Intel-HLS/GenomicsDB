#include "known_field_info.h"

std::string g_vcf_SPANNING_DELETION="*";
//Global vector storing info of all known fields
std::vector<KnownFieldInfo> g_known_field_enum_to_info;
//Known field enum to name vector
std::vector<std::string> g_known_variant_field_names = std::vector<std::string>{
    "END",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "BaseQRankSum",
    "ClippingRankSum",
    "MQRankSum",
    "ReadPosRankSum",
    "DP",
    "MQ",
    "MQ0",
    "DP_FORMAT",
    "MIN_DP",
    "GQ",
    "SB",
    "AD",
    "PL",
    "AF",
    "AN",
    "AC",
    "GT",
    "PS",
    "PGT",
    "PID"
};
//Known field name to enum
std::unordered_map<std::string, unsigned> g_known_variant_field_name_to_enum;
//KnownFieldInitializer object
KnownFieldInitializer g_known_field_initializer;

//KnownFieldInfo functions
/*
 * Check whether the known field requires a special creator
 */
bool KnownFieldInfo::requires_special_creator(unsigned enumIdx)
{
  assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
  return (g_known_field_enum_to_info[enumIdx].m_field_creator.get() != 0);
}
const std::shared_ptr<VariantFieldCreatorBase>& KnownFieldInfo::get_field_creator(unsigned enumIdx)
{
  assert(enumIdx < g_known_field_enum_to_info.size());
  return g_known_field_enum_to_info[enumIdx].m_field_creator;
}
/*
 * Function that determines whether length of the field is dependent on the #alleles 
 */
bool KnownFieldInfo::is_length_allele_dependent(unsigned enumIdx)
{
  assert(enumIdx < g_known_field_enum_to_info.size());
  return g_known_field_enum_to_info[enumIdx].is_length_allele_dependent();
}
/*
 * Check whether the known field requires ploidy - e.g. GT, GQ etc
 */
bool KnownFieldInfo::ploidy_required_for_known_field_enum(unsigned enumIdx)
{
  assert(enumIdx >= 0 && enumIdx < GVCF_NUM_KNOWN_FIELDS);
  return g_known_field_enum_to_info[enumIdx].ploidy_required_for_known_field_enum();
}

unsigned KnownFieldInfo::get_num_elements_for_known_field_enum(unsigned known_field_enum,
    unsigned num_ALT_alleles, unsigned ploidy)
{
  assert(known_field_enum < g_known_field_enum_to_info.size());
  return g_known_field_enum_to_info[known_field_enum].get_num_elements_for_known_field_enum(num_ALT_alleles, ploidy);
}

unsigned KnownFieldInfo::get_length_descriptor_for_known_field_enum(unsigned known_field_enum)
{
  assert(known_field_enum < g_known_field_enum_to_info.size());
  return g_known_field_enum_to_info[known_field_enum].get_length_descriptor();
}

bool KnownFieldInfo::get_known_field_enum_for_name(const std::string& field_name, unsigned& known_field_enum)
{
  auto iter = g_known_variant_field_name_to_enum.find(field_name);
  if(iter == g_known_variant_field_name_to_enum.end())
    return false;
  known_field_enum = (*iter).second;
  return true;
}

std::string KnownFieldInfo::get_known_field_name_for_enum(unsigned known_field_enum)
{
  assert(known_field_enum < GVCF_NUM_KNOWN_FIELDS);
  return g_known_variant_field_names[known_field_enum];
}

//Non-static member function
unsigned KnownFieldInfo::get_num_elements_for_known_field_enum(unsigned num_ALT_alleles, unsigned ploidy) const
{
  unsigned length = 0u;
  unsigned num_alleles = num_ALT_alleles + 1u;
  switch(m_length_descriptor)
  {
    case BCF_VL_FIXED:
      length = m_num_elements;
      break;
    case BCF_VL_VAR:
      length = 1u;      //function that reads from tile will know what to do
      break;
    case BCF_VL_A:
      length = num_ALT_alleles;
      break;
    case BCF_VL_R:
      length = num_alleles;
      break;
    case BCF_VL_G:
      length = (num_alleles*(num_alleles+1))/2;
      break;
    case BCF_VL_P:
      length = ploidy;
      break;
    default:
      std::cerr << "Unknown length descriptor "<<m_length_descriptor<<" - ignoring\n";
      break;
  }
  return length;
}

//KnownFieldInitializer constructor
KnownFieldInitializer::KnownFieldInitializer()
{
  for(auto i=0u;i<g_known_variant_field_names.size();++i)
    g_known_variant_field_name_to_enum[g_known_variant_field_names[i]] = i;
  //Mapping from known_field enum
  g_known_field_enum_to_info.resize(GVCF_NUM_KNOWN_FIELDS);
  //set length descriptors and creator objects for special attributes
  for(auto i=0u;i<g_known_field_enum_to_info.size();++i)
    initialize_length_descriptor(i);
}

void KnownFieldInitializer::initialize_length_descriptor(unsigned idx) const
{
  switch(idx)
  {
    case GVCF_REF_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      break;
    case GVCF_ALT_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      g_known_field_enum_to_info[idx].m_field_creator = std::shared_ptr<VariantFieldCreatorBase>(new VariantFieldCreator<VariantFieldALTData>());
      break;
    case GVCF_FILTER_IDX: 
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      break;
    case GVCF_AF_IDX:
    case GVCF_AC_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_A;
      break;
    case GVCF_AD_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_R;
      break;
    case GVCF_PL_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_G;
      break;
    case GVCF_GT_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_P;
      //Ploidy requirements
      g_known_field_enum_to_info[idx].m_ploidy_required = true;
      break;
    case GVCF_SB_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      g_known_field_enum_to_info[idx].m_num_elements = 4u;
      break;
    case GVCF_PGT_IDX:
    case GVCF_PID_IDX:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_VAR;
      break;
    default:
      g_known_field_enum_to_info[idx].m_length_descriptor = BCF_VL_FIXED;
      g_known_field_enum_to_info[idx].m_num_elements = 1u;
      break;
  }
}

bool VariantUtils::contains_deletion(const std::string& REF, const std::vector<std::string>& ALT_vec)
{
  auto REF_length = REF.length();
  if(REF_length <= 1u)
    return false;
  for(auto& alt_allele : ALT_vec)
    if(!(IS_NON_REF_ALLELE(alt_allele)) && alt_allele.length() < REF_length)
      return true;
  return false;
}
