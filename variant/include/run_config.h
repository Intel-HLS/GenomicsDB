#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"
#include "vcf_adapter.h"

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
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, int rank=0);
    rapidjson::Document m_json;
    std::string m_workspace;
    std::string m_array_name;
};

extern RunConfig g_run_config;

class VCFAdapterRunConfig : public RunConfig
{
  public:
    VCFAdapterRunConfig() : RunConfig()
    {
      m_sqlite_filename = "";
      m_vcf_header_filename = "";
    }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, VCFAdapter& vcf_adapter, int rank=0);
    std::string m_sqlite_filename;
    std::string m_vcf_header_filename;
};

#endif
