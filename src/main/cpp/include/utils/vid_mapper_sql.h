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
#ifndef GENOMICSDB_VID_MAPPER_SQL_H
#define GENOMICSDB_VID_MAPPER_SQL_H

#include "vid_mapper.h"
#include <string>
#include <vector>
#include <map>
#include <dbi/dbi.h>

const std::string PGSQL_DRIVER = "pgsql";
const std::string DBCONN_HOST = "host";
const std::string DBCONN_USERNAME = "username";
const std::string DBCONN_PASSWORD = "password";
const std::string DBCONN_DBNAME = "dbname";
const std::string DBCONN_ENCODING = "encoding";
const std::string UTF8_ENCODING = "UTF-8";

const std::string DBCONN_NEW_CONNECTION_FAILED = "FAILED to Create new Connection";
const std::string DBCONN_CONNECT_To_DB_FAILED = "FAILED to Connect to DB";
const std::string DBQUERY_FAILED = "DB Query Failed";

const std::string DBTABLE_REFERENCE_COLUMN_NAME = "name";
const std::string DBTABLE_REFERENCE_COLUMN_OFFSET = "tiledb_column_offset";
const std::string DBTABLE_REFERENCE_COLUMN_LENGTH = "length";

const std::string DBTABLE_CALLSET_COLUMN_NAME = "name";
const std::string DBTABLE_CALLSET_COLUMN_ID = "id";

const std::string DBTABLE_FIELD_COLUMN_ID = "id";
const std::string DBTABLE_FIELD_COLUMN_NAME = "field";
const std::string DBTABLE_FIELD_COLUMN_FSID = "field_set_id";
const std::string DBTABLE_FIELD_COLUMN_TYPE = "type";
const std::string DBTABLE_FIELD_COLUMN_FILTER = "is_filter";
const std::string DBTABLE_FIELD_COLUMN_FORMAT = "is_format";
const std::string DBTABLE_FIELD_COLUMN_INFO = "is_info";
const std::string DBTABLE_FIELD_COLUMN_LENTYPE = "length_type";
const std::string DBTABLE_FIELD_COLUMN_LENVAL = "length_intval";
const std::string DBTABLE_FIELD_COLUMN_COMBOP = "vcf_field_combine_operation";

const std::string COMBINE_OPERATION_ERROR_PREFIX = ("VCF field combined operation 'concatenate' can only be used with fields whose length descriptors are 'VAR'; field ");
const std::string COMBINE_OPERATION_ERROR_SUFFIX = (" does not have 'VAR' as its length descriptor");

enum {
  GENOMICSDB_VID_MAPPER_SUCCESS = 0x0,
  GENOMICSDB_VID_MAPPER_FAILURE = 0x1
};

/**
 * This class encapsulates parameters needed to establish a DB connection and
 * also identify the workspace and db_array.
 */
class SQLVidMapperRequest {
  public:
    std::string host_name;
    std::string user_name;
    std::string pass_word;
    std::string db_name;
    std::string work_space;
    std::string array_name;
};

/**
 * This class will be used to load VID Mapper information from DB. In future it
 * may also be used to modify DB info.
 */
class SQLBasedVidMapper : public VidMapper {
  public:
    SQLBasedVidMapper(const SQLVidMapperRequest&);

    /* following two declarations disable copy semantics since not needed */
    SQLBasedVidMapper(const SQLBasedVidMapper&) = delete;

    void operator=(const SQLBasedVidMapper&) = delete;

    /* following two declarations disable move semantics since not needed */
    SQLBasedVidMapper(SQLBasedVidMapper&&) = delete;

    SQLBasedVidMapper&& operator=(SQLBasedVidMapper&&) = delete;

    /* public method to load all mapping data from DB */
    int load_mapping_data_from_db();

    /* following three methods return contig mappings */
    std::vector<ContigInfo>& get_contigs() { return(m_contig_idx_to_info); }

    std::vector<std::pair<int64_t, int>>& get_contig_begin_offsets() {
      return(m_contig_begin_2_idx);
    }

    std::vector<std::pair<int64_t, int>>& get_contig_end_offsets() {
      return(m_contig_end_2_idx);
    }

    /* following two methods return callset mappings */
    std::vector<CallSetInfo>& get_callsets() { return(m_row_idx_to_info); }

    std::unordered_map<std::string, int64_t>& get_callset_name_to_idx_map() {
      return(m_callset_name_to_row_idx);
    }

    /* following two methods return field mappings */
    std::vector<FieldInfo>& get_fields() { return(m_field_idx_to_info); }

    std::unordered_map<std::string, int>& get_field_name_to_idx_map() {
      return(m_field_name_to_idx);
    }

    ~SQLBasedVidMapper();
  protected:
    dbi_conn m_conn;
    dbi_inst m_instance;
    std::string m_work_space;
    std::string m_array_name;

    /* method to retrieve contig info from DB and load into VID mapper */
    int load_contig_info();

    /* method to retrieve callset info from DB and load into VID mapper */
    int load_callset_info();

    /* method to retrieve field info from DB and load into VID mapper */
    int load_field_info();
};

class SQLBasedVidMapperException : public std::exception {
  public:
    /**
     * Default constructor
     */
    SQLBasedVidMapperException(const std::string m="");

    /*
     * Destructor: must be called before end of scope of a
     * VidMapper exception object
     */
    ~SQLBasedVidMapperException();

    /*
     * ACCESSORS
     */

    /**
     * Returns the exception message.
     */
    const char* what() const noexcept { return msg_.c_str(); }

  private:
    std::string msg_;
};

#endif  // GENOMICSDB_VID_MAPPER_SQL_H
#endif
