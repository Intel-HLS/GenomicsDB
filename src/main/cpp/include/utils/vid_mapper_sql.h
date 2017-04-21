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

enum {
  GENOMICSDB_VID_MAPPER_SUCCESS = 0x0,
  GENOMICSDB_VID_MAPPER_FAILURE = 0x1
};

class SQLVidMapperRequest {
  public:
    std::string host_name;
    std::string user_name;
    std::string pass_word;
    std::string db_name;
    std::string work_space;
    std::string array_name;
};

class SQLBasedVidMapper : public VidMapper {
  public:
    SQLBasedVidMapper(const SQLVidMapperRequest&);

    int load_contig_info();

    int load_mapping_data_from_db();

    std::vector<ContigInfo> get_contigs() { return(m_contig_idx_to_info); }

    std::vector<std::pair<int64_t, int>> get_contig_begin() {
      return(m_contig_begin_2_idx);
    }

    std::vector<std::pair<int64_t, int>> get_contig_end() {
      return(m_contig_end_2_idx);
    }

    ~SQLBasedVidMapper();
  protected:
    dbi_conn m_conn;
    dbi_inst m_instance;
    std::string m_work_space;
    std::string m_array_name;
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
