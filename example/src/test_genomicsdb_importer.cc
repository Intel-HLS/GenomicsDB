/**
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

#include <getopt.h>
#include "headers.h"
#include "genomicsdb_importer.h"
#include <mpi.h>

class VCFStreamStruct
{
  public:
    VCFStreamStruct(const char* stream_name, const char* filename, const size_t buffer_size)
    {
      m_stream_name = stream_name;
      m_buffer.resize(buffer_size); 
      m_fptr = bcf_open(filename, "r");
      m_hdr = bcf_hdr_read(m_fptr);
      size_t num_bytes_written = 0;
      while(num_bytes_written == 0)
      {
        num_bytes_written = bcf_hdr_serialize(m_hdr, &(m_buffer[0]), 0, m_buffer.size(), 1, 1);
        if(num_bytes_written == 0u)
          m_buffer.resize(2u*m_buffer.size()+1u);
      }
      m_num_valid_bytes_in_buffer = num_bytes_written;
      m_line = bcf_init();
      m_is_line_valid = false;
      memset(&m_tmp_s, 0, sizeof(kstring_t));
    }
    VCFStreamStruct(const VCFStreamStruct& other) = delete;
    VCFStreamStruct(VCFStreamStruct&& other)
    {
      m_is_line_valid = other.m_is_line_valid;
      m_num_valid_bytes_in_buffer = other.m_num_valid_bytes_in_buffer;
      m_stream_name = std::move(other.m_stream_name);
      m_buffer = std::move(other.m_buffer);
      m_line = other.m_line;
      other.m_line = 0;
      m_fptr = other.m_fptr;
      other.m_fptr = 0;
      m_hdr = other.m_hdr;
      other.m_hdr = 0;
      m_tmp_s = other.m_tmp_s;
      memset(&(other.m_tmp_s), 0, sizeof(kstring_t));
    }
    ~VCFStreamStruct()
    {
      if(m_tmp_s.m && m_tmp_s.s)
        free(m_tmp_s.s);
      memset(&m_tmp_s, 0, sizeof(kstring_t));
      if(m_line)
        bcf_destroy(m_line);
      m_line = 0;
      if(m_hdr)
        bcf_hdr_destroy(m_hdr);
      m_hdr = 0;
      if(m_fptr)
        bcf_close(m_fptr);
      m_fptr = 0;
      m_buffer.clear();
      m_stream_name.clear();
    }
    void get_more_data(const size_t val)
    {
      m_num_valid_bytes_in_buffer = val;
      auto is_first_line_in_buffer = true;
      while(1)
      {
        if(m_is_line_valid)
        {
          auto new_offset = bcf_serialize(m_line, &(m_buffer[0]), m_num_valid_bytes_in_buffer, m_buffer.size(), 1, m_hdr, &m_tmp_s);
          if(new_offset == m_num_valid_bytes_in_buffer)
          {
            if(is_first_line_in_buffer)
            {
              m_buffer.resize(2u*m_buffer.size()+1u);
              continue;
            }
            else
              break;
          }
          m_num_valid_bytes_in_buffer = new_offset;
          is_first_line_in_buffer = false;
        }
        auto status = bcf_read(m_fptr, m_hdr, m_line);
        if(status < 0)
        {
          m_is_line_valid = false;
          break;
        }
        m_is_line_valid = true;
      }
    }
  public:
    bool m_is_line_valid;
    size_t m_num_valid_bytes_in_buffer;
    std::string m_stream_name;
    std::vector<uint8_t> m_buffer;
    vcfFile* m_fptr;
    bcf_hdr_t* m_hdr;
    bcf1_t* m_line;
    kstring_t m_tmp_s;
}; 

int main(int argc, char *argv[]) {
  //Initialize MPI environment
  auto rc = MPI_Init(0, 0);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get my world rank
  int my_world_mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_mpi_rank);
  auto num_streams = INT64_MAX;
  // Define long options
  static struct option long_options[] = 
  {
    {"buffer-size",1,0,'b'},
    {"rank",1,0,'r'},
    {"output-format",1,0,'O'},
    {"loader-json-config",1,0,'l'},
    {"num-streams",1,0,'n'},
    {0,0,0,0},
  };
  int c;
  uint64_t buffer_size = 20480u;
  std::string loader_json_config_file = "";
  while((c=getopt_long(argc, argv, "l:b:r:n:", long_options, NULL)) >= 0)
  {
    switch(c)
    {
      case 'b':
        buffer_size = strtoull(optarg, 0, 10);
        break;
      case 'r':
        my_world_mpi_rank = strtoull(optarg, 0 ,10);
        break;
      case 'l':
        loader_json_config_file = std::move(std::string(optarg));
        break;
      case 'n':
        num_streams = strtoll(optarg, 0, 0);
        break;
      default:
        std::cerr << "Unknown command line argument\n";
        exit(-1);
    }
  }
  std::vector<VCFStreamStruct> stream_vector;
  //JSON with stream-filename mappings
  if(optind < argc)
  {
    std::string filename = argv[optind];
    std::ifstream ifs(filename.c_str());
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    rapidjson::Document json_doc;
    json_doc.Parse(str.c_str());
    if(json_doc.HasParseError())
      throw RunConfigException(std::string("Syntax error in JSON file ")+filename);
    std::string stream_name;
    for(auto b=json_doc.MemberBegin(), e=json_doc.MemberEnd();b!=e;++b)
    {
      const auto& curr_obj = *b;
      stream_name = curr_obj.name.GetString();
      filename = curr_obj.value.GetString();
      stream_vector.emplace_back(stream_name.c_str(), filename.c_str(), buffer_size);
      if(stream_vector.size() >= static_cast<size_t>(num_streams))
        break;
    }
  }
  GenomicsDBImporter importer(loader_json_config_file, my_world_mpi_rank);
  for(auto i=0ull;i<stream_vector.size();++i)
    importer.add_buffer_stream(stream_vector[i].m_stream_name, VidFileTypeEnum::BCF_BUFFER_STREAM_TYPE, buffer_size,
        &(stream_vector[i].m_buffer[0]), stream_vector[i].m_num_valid_bytes_in_buffer);
  importer.setup_loader();
  auto& buffer_stream_idx_to_global_file_idx_vec = importer.get_buffer_stream_idx_to_global_file_idx_vec();
  for(auto i=0ull;i<stream_vector.size();++i)
  {
    if(buffer_stream_idx_to_global_file_idx_vec[i] < 0)
      continue;
    stream_vector[i].get_more_data(0u);
    importer.write_data_to_buffer_stream(i, 0, &(stream_vector[i].m_buffer[0]), stream_vector[i].m_num_valid_bytes_in_buffer);
  }
  auto num_iterations = 0ull;
  while(!importer.is_done())
  {
    importer.import_batch();
    const auto& stream_id_vec = importer.get_exhausted_buffer_stream_identifiers();
    for(auto stream_id : stream_id_vec)
    {
      auto stream_idx = stream_id.first;
      stream_vector[stream_idx].get_more_data(0u);
      importer.write_data_to_buffer_stream(stream_idx, 0, &(stream_vector[stream_idx].m_buffer[0]), stream_vector[stream_idx].m_num_valid_bytes_in_buffer);
    }
    ++num_iterations;
  }
  importer.finish();
#ifdef DEBUG
  std::cerr << "Num iterations in importer "<<num_iterations<<"\n";
#endif
  MPI_Finalize();
  return 0;
}
