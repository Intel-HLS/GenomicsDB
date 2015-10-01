#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include <fstream>

class RunConfig
{
  public:
    RunConfig()
    {
      m_workspace = "";
      m_array_name = "";
    }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, int rank);
    rapidjson::Document m_json;
    std::string m_workspace;
    std::string m_array_name;
};

extern RunConfig g_run_config;

#endif
