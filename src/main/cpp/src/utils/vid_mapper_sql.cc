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

#define THROW_EXCEPTION(X) throw SQLBasedVidMapperException(X);

#define VERIFY_OR_THROW(X) if(!(X)) throw SQLBasedVidMapperException(#X);

SQLBasedVidMapper::SQLBasedVidMapper(const SQLVidMapperRequest& request) : VidMapper() {
  dbi_initialize_r(NULL, &m_instance);
  m_conn = dbi_conn_new_r(PGSQL_DRIVER.c_str(), m_instance);
  if (NULL == m_conn) {
    THROW_EXCEPTION(DBCONN_NEW_CONNECTION_FAILED);
  }

  dbi_conn_set_option(m_conn, DBCONN_HOST.c_str(), request.host_name.c_str());
  dbi_conn_set_option(m_conn, DBCONN_USERNAME.c_str(), request.user_name.c_str());
  dbi_conn_set_option(m_conn, DBCONN_PASSWORD.c_str(), request.pass_word.c_str());
  dbi_conn_set_option(m_conn, DBCONN_DBNAME.c_str(), request.db_name.c_str());
  dbi_conn_set_option(m_conn, DBCONN_ENCODING.c_str(), UTF8_ENCODING.c_str());

  if (dbi_conn_connect(m_conn) < 0) {
    THROW_EXCEPTION(DBCONN_CONNECT_To_DB_FAILED);
  }

  m_work_space = request.work_space;
  m_array_name = request.array_name;
}

SQLBasedVidMapper::~SQLBasedVidMapper() {
  dbi_conn_close(m_conn);
  dbi_shutdown_r(m_instance);
}

int SQLBasedVidMapper::load_contig_info() {
  std::stringstream ss;
  ss <<"select a.name, a.length, a.tiledb_column_offset ";
  ss <<"from reference a, db_array b, workspace c ";
  ss <<"where ((c.name='" <<m_work_space <<"') and (c.id = b.workspace_id) ";
  ss <<"and (b.name = '" <<m_array_name <<"') and ";
  ss <<"(a.reference_set_id = b.reference_set_id)) ";
  ss <<"order by a.tiledb_column_offset";

  std::string query = ss.str();
  std::cout <<"QUERY: <" <<query <<">\n";
  dbi_result result = dbi_conn_query(m_conn, query.c_str());
  if (NULL == result) {
      const char* errmsg = NULL;
      dbi_conn_error(m_conn, &errmsg);
      std::string error_message = (DBQUERY_FAILED + " : " + errmsg);
      THROW_EXCEPTION(error_message);
  }

  int num_contigs = (int) dbi_result_get_numrows(result);
  m_contig_idx_to_info.resize(num_contigs);
  m_contig_begin_2_idx.resize(num_contigs);
  m_contig_end_2_idx.resize(num_contigs);
  int contig_idx = -1;

  while (dbi_result_next_row(result)) {
    std::string contig_name = dbi_result_get_string(result, DBTABLE_REFERENCE_COLUMN_NAME.c_str());
    ++contig_idx;

    int64_t tiledb_column_offset = dbi_result_get_longlong(result, DBTABLE_REFERENCE_COLUMN_OFFSET.c_str());
    int64_t length = dbi_result_get_longlong(result, DBTABLE_REFERENCE_COLUMN_LENGTH.c_str());

    m_contig_name_to_idx[contig_name] = contig_idx;
    m_contig_idx_to_info[contig_idx].set_info(
                                       contig_idx,
                                       contig_name,
                                       length,
                                       tiledb_column_offset);
    m_contig_begin_2_idx[contig_idx].first = tiledb_column_offset;
    m_contig_begin_2_idx[contig_idx].second = contig_idx;
    m_contig_end_2_idx[contig_idx].first = (tiledb_column_offset + length - 1); 
    m_contig_end_2_idx[contig_idx].second = contig_idx;
  }

  dbi_result_free(result);
  return(GENOMICSDB_VID_MAPPER_SUCCESS);
}

int SQLBasedVidMapper::load_callset_info() {
  std::stringstream ss;
  ss <<"select a.id, a.name from callset a, ";
  ss <<"db_array b, workspace c, callset_to_db_array_association d ";
  ss <<"where ((c.name='" <<m_work_space <<"') and (c.id = b.workspace_id) ";
  ss <<"and (b.name = '" <<m_array_name <<"') and ";
  ss <<"(b.id = d.db_array_id) and (d.callset_id = a.id))";

  std::string query = ss.str();
  std::cout <<"QUERY: <" <<query <<">\n";
  dbi_result result = dbi_conn_query(m_conn, query.c_str());
  if (NULL == result) {
      const char* errmsg = NULL;
      dbi_conn_error(m_conn, &errmsg);
      std::string error_message = (DBQUERY_FAILED + " : " + errmsg);
      THROW_EXCEPTION(error_message);
  }

  int num_callsets = (int) dbi_result_get_numrows(result);
  m_row_idx_to_info.resize(num_callsets);
  m_max_callset_row_idx = -1;

  while (dbi_result_next_row(result)) {
    std::string callset_name = dbi_result_get_string(result, DBTABLE_CALLSET_COLUMN_NAME.c_str());
    int64_t row_idx = dbi_result_get_longlong(result, DBTABLE_CALLSET_COLUMN_ID.c_str());
    int64_t file_idx = row_idx;
    m_max_callset_row_idx = std::max(m_max_callset_row_idx, row_idx);
    m_callset_name_to_row_idx[callset_name] = row_idx;
    m_row_idx_to_info[row_idx].set_info(row_idx, callset_name, file_idx, 0);
  }

  dbi_result_free(result);
  return(GENOMICSDB_VID_MAPPER_SUCCESS);
}

int SQLBasedVidMapper::load_mapping_data_from_db() {
  load_contig_info();
  load_callset_info();
  return(GENOMICSDB_VID_MAPPER_SUCCESS);
}

SQLBasedVidMapperException::~SQLBasedVidMapperException() {}

SQLBasedVidMapperException::SQLBasedVidMapperException(const std::string m) : msg_("SQLBasedVidMapperException : "+m) {}
#endif

