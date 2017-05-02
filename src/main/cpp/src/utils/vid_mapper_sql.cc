/*
 * The MIT License (MIT)
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

#ifdef LIBDBI
#include "vid_mapper_sql.h"
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include <sstream>

#define VERIFY_OR_THROW(X) if(!(X)) throw SQLBasedVidMapperException(#X);

SQLBasedVidMapper::SQLBasedVidMapper(const SQLVidMapperRequest& request) : VidMapper() {
  dbi_initialize_r(NULL, &m_instance);
  m_conn = dbi_conn_new_r("pgsql", m_instance);
  VERIFY_OR_THROW(m_conn != NULL);

  dbi_conn_set_option(m_conn, "host", request.host_name.c_str());
  dbi_conn_set_option(m_conn, "username", request.user_name.c_str());
  dbi_conn_set_option(m_conn, "password", request.pass_word.c_str());
  dbi_conn_set_option(m_conn, "dbname", request.db_name.c_str());
  dbi_conn_set_option(m_conn, "encoding", "UTF-8");

  if (dbi_conn_connect(m_conn) < 0) {
    VERIFY_OR_THROW(1 < 0);
  }

  m_work_space = request.work_space;
  m_array_name = request.array_name;
}

SQLBasedVidMapper::~SQLBasedVidMapper() {
  dbi_conn_close(m_conn);
  dbi_shutdown_r(m_instance);
}

int SQLBasedVidMapper::load_contig_info() {
  auto duplicate_contigs_exist = false;
  std::string contig_name;
  int contig_idx = -1;

  std::stringstream ss;
  ss <<"select a.name, a.length, a.tiledb_column_offset ";
  ss <<"from reference a, db_array b, workspace c ";
  ss <<"where ((c.name='" <<m_work_space <<"') and (c.id = b.workspace_id) ";
  ss <<"and (b.name = '" <<m_array_name <<"') and ";
  ss <<"(a.reference_set_id = b.reference_set_id))";

  std::string query = ss.str();
  std::cout <<"QUERY: <" <<query <<">\n";
  dbi_result result = dbi_conn_query(m_conn, query.c_str());
  VERIFY_OR_THROW(result != NULL);

  int num_contigs = (int) dbi_result_get_numrows(result);
  m_contig_idx_to_info.resize(num_contigs);
  m_contig_begin_2_idx.resize(num_contigs);
  m_contig_end_2_idx.resize(num_contigs);

  while (dbi_result_next_row(result)) {
    ++contig_idx;
    contig_name = dbi_result_get_string(result, "name");

    /**
     * There is probably a constraint which prevents the following
     */
    if (m_contig_name_to_idx.find(contig_name) != m_contig_name_to_idx.end()) {
      std::cerr << "Contig/chromosome name "
                << contig_name
                << " appears more than once in reference table\n";
      duplicate_contigs_exist = true;
      continue;
    }

    int64_t tiledb_column_offset = dbi_result_get_longlong(result, "tiledb_column_offset");
    VERIFY_OR_THROW(tiledb_column_offset >= 0LL);
    int64_t length = dbi_result_get_longlong(result, "length");
    VERIFY_OR_THROW(length >= 0LL);

    m_contig_name_to_idx[contig_name] = contig_idx;
    m_contig_idx_to_info[contig_idx].set_info(
                                       contig_idx,
                                       contig_name,
                                       length,
                                       tiledb_column_offset);
    m_contig_begin_2_idx[contig_idx].first = tiledb_column_offset;
    m_contig_begin_2_idx[contig_idx].second = contig_idx;
    m_contig_end_2_idx[contig_idx].first =
        tiledb_column_offset + length - 1; //inclusive
    m_contig_end_2_idx[contig_idx].second = contig_idx;

    if (duplicate_contigs_exist) {
        throw SQLBasedVidMapperException(
          std::string("Duplicate contigs found: ")
          + contig_name);
    }

    std::sort(
        m_contig_begin_2_idx.begin(),
        m_contig_begin_2_idx.end(),
        contig_offset_idx_pair_cmp);
    std::sort(
        m_contig_end_2_idx.begin(),
        m_contig_end_2_idx.end(),
        contig_offset_idx_pair_cmp);
  }

  // Check that there are no spurious overlaps.
  // If found, throw an exception
  auto last_contig_idx = -1;
  auto last_contig_end_column = -1ll;
  auto overlapping_contigs_exist = false;

  for (auto contig_idx = 0UL; contig_idx < m_contig_begin_2_idx.size(); ++contig_idx) {
    const auto& contig_info = m_contig_idx_to_info[contig_idx];

    if (last_contig_idx >= 0) {
      const auto& last_contig_info = m_contig_idx_to_info[last_contig_idx];
      if (contig_info.m_tiledb_column_offset <= last_contig_end_column) {
        std::cerr << "Contig/chromosome "
                  << contig_info.m_name
                  << " begins at TileDB column "
                  << contig_info.m_tiledb_column_offset
                  << " and intersects with contig/chromosome "
                  << last_contig_info.m_name
                  << " that spans columns ["
                  << last_contig_info.m_tiledb_column_offset
                  << ", "
                  << last_contig_info.m_tiledb_column_offset +
                     last_contig_info.m_length-1
                  << "]" << "\n";
        overlapping_contigs_exist = true;
      }
    }

    last_contig_idx = contig_idx;
    last_contig_end_column = (contig_info.m_tiledb_column_offset + contig_info.m_length - 1);
  }

  if (overlapping_contigs_exist) {
    throw SQLBasedVidMapperException(std::string("Overlapping contigs found"));
  }

  dbi_result_free(result);
  return(GENOMICSDB_VID_MAPPER_SUCCESS);
}

int SQLBasedVidMapper::load_mapping_data_from_db() {
  /**
   * ret = parse_contigs_from_vidmap(vid_map_protobuf);
   * assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);
   * ret = parse_infofields_from_vidmap(vid_map_protobuf);
   * assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);
   */
  load_contig_info();
  return(GENOMICSDB_VID_MAPPER_SUCCESS);
}

SQLBasedVidMapperException::~SQLBasedVidMapperException() {}

SQLBasedVidMapperException::SQLBasedVidMapperException(const std::string m) : msg_("SQLBasedVidMapperException : "+m) {}
#endif

