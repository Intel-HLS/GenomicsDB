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

#include "genomicsdb_importer.h"

#define VERIFY_OR_THROW(X) if(!(X)) throw GenomicsDBImporterException(#X);

void GenomicsDBImporter::add_buffer_stream(const std::string& name, const VidFileTypeEnum buffer_stream_type,
        const size_t capacity, const uint8_t* initialization_buffer, const size_t num_bytes_in_initialization_buffer)
{
  if(m_is_loader_setup)
    throw GenomicsDBImporterException(std::string("Cannot add buffer stream once setup_loader() has been called for a given GenomicsDBImporter object"));
  m_buffer_stream_info_vec.emplace_back();
  auto& curr_buffer_stream_info = m_buffer_stream_info_vec[m_buffer_stream_info_vec.size()-1u];
  curr_buffer_stream_info.m_name = name;
  curr_buffer_stream_info.m_type = buffer_stream_type;
  curr_buffer_stream_info.m_buffer_stream_idx = m_buffer_stream_info_vec.size()-1u;
  curr_buffer_stream_info.m_buffer_capacity = capacity;
  if(initialization_buffer && num_bytes_in_initialization_buffer)
  {
    curr_buffer_stream_info.m_initialization_buffer.resize(num_bytes_in_initialization_buffer);
    memcpy(&(curr_buffer_stream_info.m_initialization_buffer[0]), initialization_buffer, num_bytes_in_initialization_buffer);
    curr_buffer_stream_info.m_initialization_buffer_num_valid_bytes = num_bytes_in_initialization_buffer;
  }
}

void GenomicsDBImporter::copy_simple_members(const GenomicsDBImporter& other)
{
  m_is_loader_setup = other.m_is_loader_setup;
  m_rank = other.m_rank;
  m_lb_callset_row_idx = other.m_lb_callset_row_idx;
  m_ub_callset_row_idx = other.m_ub_callset_row_idx;
}

GenomicsDBImporter::GenomicsDBImporter(GenomicsDBImporter&& other)
{
  copy_simple_members(other);
  //Move-in members
  m_loader_config_file = std::move(other.m_loader_config_file);
  m_buffer_stream_info_vec = std::move(other.m_buffer_stream_info_vec);
  m_loader_ptr = other.m_loader_ptr;
  other.m_loader_ptr = 0;
  m_read_state = other.m_read_state;
  other.m_read_state = 0;
}

GenomicsDBImporter::~GenomicsDBImporter()
{
  m_loader_config_file.clear();
  m_buffer_stream_info_vec.clear();
  if(m_loader_ptr)
    delete m_loader_ptr;
  m_loader_ptr = 0;
  if(m_read_state)
    delete m_read_state;
  m_read_state = 0;
}

void GenomicsDBImporter::setup_loader()
{
  if(m_is_loader_setup) //already setup
    return;
  m_loader_ptr = new VCF2TileDBLoader(m_loader_config_file, m_rank,
      m_buffer_stream_info_vec, m_lb_callset_row_idx, m_ub_callset_row_idx);
  m_read_state = m_loader_ptr->construct_read_state_object();
  m_is_loader_setup = true;
}

void GenomicsDBImporter::import_batch()
{
  if(!m_is_loader_setup)
    throw GenomicsDBImporterException(std::string("Cannot import data till setup_loader() has been called for a given GenomicsDBImporter object"));
  m_loader_ptr->read_all(*m_read_state);
}
