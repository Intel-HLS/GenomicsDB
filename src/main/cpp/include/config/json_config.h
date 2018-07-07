/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
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

#ifndef JSON_CONFIG_H
#define JSON_CONFIG_H

#include "genomicsdb_config_base.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/prettywriter.h"

class JSONConfigBase : public GenomicsDBConfigBase
{
  public:
    JSONConfigBase()
      : GenomicsDBConfigBase()
    {}
    JSONConfigBase(const GenomicsDBConfigBase& x)
      : GenomicsDBConfigBase(x)
    {}
    static void extract_contig_interval_from_object(const rapidjson::Value& curr_json_object,
        const VidMapper* id_mapper, ColumnRange& result);
    static bool extract_interval_from_PB_struct_or_return_false(const rapidjson::Value& curr_json_object,
        const VidMapper* id_mapper,
        ColumnRange& result);
    void read_from_file(const std::string& filename, const int rank=0);
    void read_and_initialize_vid_and_callset_mapping_if_available(const int rank);
    const rapidjson::Document& get_rapidjson_doc() const { return m_json; }
  protected:
    rapidjson::Document m_json;
};

rapidjson::Document parse_json_file(const std::string& s);

#endif
