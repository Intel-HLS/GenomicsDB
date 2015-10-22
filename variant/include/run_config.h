#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"
#include "vcf_adapter.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"

//Exceptions thrown 
class RunConfigException {
  public:
    RunConfigException(const std::string m="") : msg_("RunConfigException : "+m) { ; }
    ~RunConfigException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

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

#ifdef HTSDIR

class VCFAdapterRunConfig : public RunConfig
{
  public:
    VCFAdapterRunConfig() : RunConfig()
    {
      m_sqlite_filename = "";
      m_vcf_header_filename = "";
    }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, VCFAdapter& vcf_adapter,
        std::string output_format="", int rank=0);
    std::string m_sqlite_filename;
    std::string m_vcf_header_filename;
    std::string m_reference_genome;
    std::string m_vcf_output_filename;
};

#endif

#endif
