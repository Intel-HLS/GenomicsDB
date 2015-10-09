#ifdef HTSDIR

#include "vcf_adapter.h"

bool contig_offset_idx_pair_cmp(const std::pair<int64_t, int>& first, const std::pair<int64_t, int>& second)
{
  return (first.first < second.first);
}

VCFAdapter::VCFAdapter()
{
  m_sqlite_filename = "";
  m_vcf_header_filename = "";
  m_template_vcf_hdr = 0;
  memset(&m_sqlite_mapping_info, 0, sizeof(sqlite_mappings_struct));
  clear();
}

VCFAdapter::~VCFAdapter()
{
  clear();
  if(m_template_vcf_hdr)
    bcf_hdr_destroy(m_template_vcf_hdr);
  free_sqlite3_data(&m_sqlite_mapping_info);
}

void VCFAdapter::clear()
{
  m_contig_begin_2_idx.clear();
  m_contig_end_2_idx.clear();
}

void VCFAdapter::initialize(const std::string& sqlite_filename, const std::string& vcf_header_filename)
{
  //Sqlite file for mapping ids to names etc
  m_sqlite_filename = sqlite_filename;
  strcpy(m_sqlite_mapping_info.sqlite_file, m_sqlite_filename.c_str());
  open_sqlite3_db(m_sqlite_filename.c_str(), &(m_sqlite_mapping_info.db));
  read_all_from_sqlite(&m_sqlite_mapping_info);
  //Create sorted vector for contig begin, end - useful in querying contig given a position
  m_contig_begin_2_idx.resize(m_sqlite_mapping_info.m_num_contigs);
  m_contig_end_2_idx.resize(m_sqlite_mapping_info.m_num_contigs);
  //Check if sorting needed, or already sorted
  bool sort_needed = false;
  int64_t last_contig_offset = 0;
  for(auto i=0;i<m_sqlite_mapping_info.m_num_contigs;++i)
  {
    auto contig_offset = m_sqlite_mapping_info.input_contig_idx_2_offset[i];
    auto contig_length = m_sqlite_mapping_info.m_contig_lengths[i];
    m_contig_begin_2_idx[i].first = contig_offset;
    m_contig_begin_2_idx[i].second = i;
    m_contig_end_2_idx[i].first = contig_offset + contig_length - 1; //-1 for [begin, end] not [begin, end)
    m_contig_end_2_idx[i].second = i;
    if(contig_offset < last_contig_offset)
      sort_needed = true;
    last_contig_offset = contig_offset;
  }
  if(sort_needed)
  {
    std::sort(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), contig_offset_idx_pair_cmp);
    std::sort(m_contig_end_2_idx.begin(), m_contig_end_2_idx.end(), contig_offset_idx_pair_cmp);
  }
  //Read template header with fields and contigs
  m_vcf_header_filename = vcf_header_filename;
  auto* fptr = bcf_open(vcf_header_filename.c_str(), "r");
  m_template_vcf_hdr = bcf_hdr_read(fptr);
  bcf_close(fptr);
}
    
bool VCFAdapter::get_contig_location(int64_t query_position, std::string& contig_name, int64_t& contig_position) const
{
  int idx = -1;
  std::pair<int64_t, int> query_pair;
  query_pair.first = query_position;
  query_pair.second = 0;
  //find contig with offset >= query_position
  auto iter = std::lower_bound(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), query_pair, contig_offset_idx_pair_cmp);
  if(iter == m_contig_begin_2_idx.end())        //no such contig exists, hence, get last contig in sorted order
  {
    assert(m_contig_begin_2_idx.size() > 0u);
    idx = m_contig_begin_2_idx[m_contig_begin_2_idx.size()-1u].second;
  }
  else
  {
    if((*iter).first == query_position)       //query_position == contig offset at iter, found idx
      idx = (*iter).second;
    else                                //query_position < contig_offset at iter, get idx at previous element
    {
      assert(iter != m_contig_begin_2_idx.begin());     //if iter == begin(), horribly wrong
      auto vector_idx = iter - m_contig_begin_2_idx.begin();
      idx = m_contig_begin_2_idx[vector_idx-1].second;
    }
  }
  if(idx < 0)
    return false;
  assert(idx < m_sqlite_mapping_info.m_num_contigs);
  //query_position is within the contig
  auto contig_offset = m_sqlite_mapping_info.input_contig_idx_2_offset[idx];
  auto contig_length = m_sqlite_mapping_info.m_contig_lengths[idx];
  if((query_position >= contig_offset) && (query_position < contig_offset+contig_length))
  {
    contig_name = m_sqlite_mapping_info.m_contig_names[idx];
    contig_position = query_position - contig_offset;
    return true;
  }
  return false;
}

bool VCFAdapter::get_next_contig_location(int64_t query_position, std::string& next_contig_name, int64_t& next_contig_offset) const
{
  int idx = -1;
  std::pair<int64_t, int> query_pair;
  query_pair.first = query_position;
  query_pair.second = 0;
  //find contig with offset > query_position
  auto iter = std::upper_bound(m_contig_begin_2_idx.begin(), m_contig_begin_2_idx.end(), query_pair, contig_offset_idx_pair_cmp);
  if(iter == m_contig_begin_2_idx.end())        //no such contig exists, hence, set large upper bound
  {
    next_contig_name = "";
    next_contig_offset = INT64_MAX;
  }
  else
  {
    idx = (*iter).second;
    assert(idx >=0 && idx < m_sqlite_mapping_info.m_num_contigs);
    next_contig_name = m_sqlite_mapping_info.m_contig_names[idx];
    next_contig_offset = m_sqlite_mapping_info.input_contig_idx_2_offset[idx];
    assert(next_contig_offset > query_position);
  }
  return true;
}

#endif //ifdef HTSDIR
