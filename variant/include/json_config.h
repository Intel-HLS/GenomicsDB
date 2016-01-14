#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

#include "variant_query_config.h"
#include "vcf_adapter.h"
#include "vid_mapper.h"

#include "rapidjson/document.h"
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"

//Exceptions thrown 
class RunConfigException : public std::exception {
  public:
    RunConfigException(const std::string m="") : msg_("RunConfigException : "+m) { ; }
    ~RunConfigException() { ; }
    // ACCESSORS
    /** Returns the exception message. */
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

class JSONConfigBase
{
  public:
    JSONConfigBase() { ; }
    void read_from_file(const std::string& filename);
  protected:
    rapidjson::Document m_json;
};

class JSONBasicQueryConfig : public JSONConfigBase
{
  public:
    JSONBasicQueryConfig() : JSONConfigBase()
    {
      m_workspace = "";
      m_array_name = "";
    }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config, int rank=0);
  public:
    std::string m_workspace;
    std::string m_array_name;
};

#ifdef HTSDIR

class JSONVCFAdapterConfig : public JSONConfigBase
{
  public:
    JSONVCFAdapterConfig() : JSONConfigBase()
    {
      m_vcf_header_filename = "";
    }
    void read_from_file(const std::string& filename,
        VCFAdapter& vcf_adapter, std::string output_format="", int rank=0);
  protected:
    std::string m_vcf_header_filename;
    std::string m_reference_genome;
    std::string m_vcf_output_filename;
};

class JSONVCFAdapterQueryConfig : public JSONVCFAdapterConfig, public JSONBasicQueryConfig
{
  public:
    JSONVCFAdapterQueryConfig() : JSONVCFAdapterConfig(), JSONBasicQueryConfig() { ; }
    void read_from_file(const std::string& filename, VariantQueryConfig& query_config,
        VCFAdapter& vcf_adapter, FileBasedVidMapper& id_mapper,
        std::string output_format="", int rank=0);
};



#endif

#endif
