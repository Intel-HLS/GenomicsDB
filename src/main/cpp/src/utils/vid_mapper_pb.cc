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

#include "vid_mapper_pb.h"
#include "known_field_info.h"
#include <ctype.h>
#include <cstdlib>
#include <cstdio>

#define VERIFY_OR_THROW(X) if(!(X)) throw ProtoBufBasedVidMapperException(#X);

GenomicsDBProtoBufInitAndCleanup g_genomicsdb_protobuf_init_and_cleanup;

ProtoBufBasedVidMapper::ProtoBufBasedVidMapper(
  const VidMappingPB* vid_map_protobuf,
  const CallsetMappingPB* callset_map_protobuf,
  const std::vector<BufferStreamInfo>& buffer_stream_info_vec)
  : VidMapper() {

  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  assert (vid_map_protobuf->IsInitialized() &&
      vid_map_protobuf->chromosomes_size()!=0 &&
      vid_map_protobuf->infofields_size()!=0);

  assert (callset_map_protobuf->IsInitialized() &&
      callset_map_protobuf->callset_map_size()!= 0);

  initialize(vid_map_protobuf, callset_map_protobuf, buffer_stream_info_vec);
}

void ProtoBufBasedVidMapper::initialize(
  const VidMappingPB* vid_map_protobuf,
  const CallsetMappingPB* callset_map_protobuf,
  const std::vector<BufferStreamInfo>& buffer_stream_info_vec) {

  int ret = 0;
  ret = parse_callset_protobuf(
          callset_map_protobuf,
          buffer_stream_info_vec);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);
  ret = parse_vidmap_protobuf(vid_map_protobuf);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);

  m_is_initialized = true;
}

int ProtoBufBasedVidMapper::parse_callset_protobuf(
  const CallsetMappingPB* callset_map_protobuf,
  const std::vector<BufferStreamInfo>& buffer_stream_info_vec) {

  auto num_callsets = callset_map_protobuf->callset_map_size();
  m_row_idx_to_info.resize(num_callsets);
  m_max_callset_row_idx = -1;

  google::protobuf::Map<std::string, SampleIDToTileDBIDMap>::const_iterator it =
      callset_map_protobuf->callset_map().cbegin();

  for (; it != callset_map_protobuf->callset_map().cend(); ++it) {
    SampleIDToTileDBIDMap sample_info = it->second;
    std::string callset_name = it->first;
    int64_t row_idx = sample_info.tiledb_row_index();
    int64_t idx_in_file = sample_info.idx_in_file();
    std::string stream_name = sample_info.stream_name();

    if(m_callset_name_to_row_idx.find(callset_name) !=
        m_callset_name_to_row_idx.end()) {

      if(m_callset_name_to_row_idx[callset_name] != row_idx) {
        throw ProtoBufBasedVidMapperException(
            std::string("ERROR: Callset/sample ")
            + callset_name
            + " have two TileDB row indexes: "
            + std::to_string(m_callset_name_to_row_idx[callset_name])
            + ", "
            + std::to_string(row_idx));
      }
    }

    m_max_callset_row_idx = std::max(m_max_callset_row_idx, row_idx);

    auto file_idx = get_or_append_global_file_idx(stream_name);
    const auto& curr_row_info = m_row_idx_to_info[row_idx];

    // If row information for the same callset is found,
    // check whether it conflicts with current information.
    // If so, throw error
    if(curr_row_info.m_is_initialized) {
      if(curr_row_info.m_file_idx >= 0 && file_idx >= 0
          && curr_row_info.m_file_idx != file_idx)
        throw ProtoBufBasedVidMapperException(
          std::string("Callset/sample ")
          + callset_name
          + " has multiple stream names: "
          + m_file_idx_to_info[curr_row_info.m_file_idx].m_name
          + ", "
          + stream_name);
      if(curr_row_info.m_idx_in_file != idx_in_file)
        throw ProtoBufBasedVidMapperException(
          std::string("Conflicting values of \"sample_vcf_index\" ")
          + std::string("specified for Callset/sample ")
          + callset_name
          + ": "
          + std::to_string(curr_row_info.m_idx_in_file)
          + ", "
          + std::to_string(idx_in_file));
    } else {
      assert(file_idx < static_cast<int64_t>(m_file_idx_to_info.size()));
      m_file_idx_to_info[file_idx].add_local_tiledb_row_idx_pair(
          idx_in_file,
          row_idx);
    }

    m_callset_name_to_row_idx[callset_name] = row_idx;
    VERIFY_OR_THROW(static_cast<size_t>(row_idx) < m_row_idx_to_info.size());
    if(m_row_idx_to_info[row_idx].m_is_initialized &&
        m_row_idx_to_info[row_idx].m_name != callset_name)
          throw ProtoBufBasedVidMapperException(
            std::string("Callset/sample ")
            + callset_name
            + " has the same TileDB row index as "
            + m_row_idx_to_info[row_idx].m_name);

    m_row_idx_to_info[row_idx].set_info(
      row_idx,
      callset_name,
      file_idx,
      idx_in_file);
  }  // end of for all callsets

  //For buffer streams
  m_buffer_stream_idx_to_global_file_idx.resize(
    buffer_stream_info_vec.size(), -1);

  auto max_buffer_stream_idx_with_global_file_idx = -1ll;
  for(auto i=0ull;i<buffer_stream_info_vec.size();++i)
  {
    const auto& info = buffer_stream_info_vec[i];
    int64_t global_file_idx;
    auto found = get_global_file_idx(info.m_name, global_file_idx);
    if(found)
    {
      auto& curr_file_info = m_file_idx_to_info[global_file_idx];
      curr_file_info.m_type = info.m_type;
      curr_file_info.m_buffer_stream_idx = i;
      curr_file_info.m_buffer_capacity = info.m_buffer_capacity;
      curr_file_info.m_initialization_buffer = info.m_initialization_buffer;
      curr_file_info.m_initialization_buffer_num_valid_bytes =
        info.m_initialization_buffer_num_valid_bytes;
      m_buffer_stream_idx_to_global_file_idx[i] = global_file_idx;
      max_buffer_stream_idx_with_global_file_idx = i;
    }
  }
  m_buffer_stream_idx_to_global_file_idx.resize(
    max_buffer_stream_idx_with_global_file_idx+1);

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_callset_protobuf

int ProtoBufBasedVidMapper::parse_vidmap_protobuf(
  const VidMappingPB* vid_map_protobuf) {

  int ret = 0;
  ret = parse_contigs_from_vidmap(vid_map_protobuf);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);
  ret = parse_infofields_from_vidmap(vid_map_protobuf);
  assert (ret == GENOMICSDB_VID_MAPPER_SUCCESS);

  return GENOMICSDB_VID_MAPPER_SUCCESS;
}

int ProtoBufBasedVidMapper::parse_contigs_from_vidmap(
  const VidMappingPB* vid_map_protobuf) {

  auto num_contigs = vid_map_protobuf->chromosomes_size();
  m_contig_idx_to_info.resize(num_contigs);
  m_contig_begin_2_idx.resize(num_contigs);
  m_contig_end_2_idx.resize(num_contigs);
  auto duplicate_contigs_exist = false;
  std::string contig_name;

  for (auto contig_idx = 0L; contig_idx < num_contigs; ++contig_idx) {
    contig_name = vid_map_protobuf->chromosomes(contig_idx).name();

    if(m_contig_name_to_idx.find(contig_name) != m_contig_name_to_idx.end())
    {
      std::cerr << "Contig/chromosome name "
                << contig_name
                << " appears more than once in vid map\n"
                << " Please check merge headers method "
                << "in GenomicsDBImporter.java\n";
      duplicate_contigs_exist = true;
      continue;
    }

    auto tiledb_column_offset =
        vid_map_protobuf->chromosomes(contig_idx).tiledb_column_offset();

    VERIFY_OR_THROW(tiledb_column_offset >= 0LL);
    auto length = vid_map_protobuf->chromosomes(contig_idx).length();

    VERIFY_OR_THROW(length >= 0LL);

    m_contig_name_to_idx[contig_name] = contig_idx;
    m_contig_idx_to_info[contig_idx].set_info(
                                       contig_idx,
                                       contig_name,
                                       length,
                                       tiledb_column_offset);
    m_contig_begin_2_idx[contig_idx].first = tiledb_column_offset;
    m_contig_begin_2_idx[contig_idx].second = contig_idx;
    m_contig_end_2_idx[contig_idx].first =
        tiledb_column_offset + length - 1; //inclusive
    m_contig_end_2_idx[contig_idx].second = contig_idx;

    if(duplicate_contigs_exist) {
        throw ProtoBufBasedVidMapperException(
          std::string("Duplicate contigs found: ")
          + contig_name);
    }

    std::sort(
        m_contig_begin_2_idx.begin(),
        m_contig_begin_2_idx.end(),
        contig_offset_idx_pair_cmp);
    std::sort(
        m_contig_end_2_idx.begin(),
        m_contig_end_2_idx.end(),
        contig_offset_idx_pair_cmp);
  }

  // Check that there are no spurious overlaps.
  // If found, throw an exception
  auto last_contig_idx = -1;
  auto last_contig_end_column = -1ll;
  auto overlapping_contigs_exist = false;
  for (auto contig_idx = 0UL; contig_idx < m_contig_begin_2_idx.size();
      ++contig_idx)
//  for (const auto& offset_idx_pair : m_contig_begin_2_idx)
  {
//    auto contig_idx = offset_idx_pair.second;
    const auto& contig_info = m_contig_idx_to_info[contig_idx];
    if(last_contig_idx >= 0)
    {
      const auto& last_contig_info = m_contig_idx_to_info[last_contig_idx];
      if(contig_info.m_tiledb_column_offset <= last_contig_end_column)
      {
        std::cerr << "Contig/chromosome "
                  << contig_info.m_name
                  << " begins at TileDB column "
                  << contig_info.m_tiledb_column_offset
                  << " and intersects with contig/chromosome "
                  << last_contig_info.m_name
                  << " that spans columns ["
                  << last_contig_info.m_tiledb_column_offset
                  << ", "
                  << last_contig_info.m_tiledb_column_offset +
                     last_contig_info.m_length-1
                  << "]" << "\n";
        overlapping_contigs_exist = true;
      }
    }
    last_contig_idx = contig_idx;
    last_contig_end_column =
        contig_info.m_tiledb_column_offset + contig_info.m_length - 1;
  }

  if(overlapping_contigs_exist) {
    throw ProtoBufBasedVidMapperException(
      std::string("Overlapping contigs found"));
  }

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_contigs_from_vidmap

int ProtoBufBasedVidMapper::parse_infofields_from_vidmap(
  const VidMappingPB* vid_map_protobuf) {

  auto num_fields = vid_map_protobuf->infofields_size();
  m_field_idx_to_info.resize(num_fields);
  std::string field_name;
  std::string field_type;
  auto duplicate_fields_exist = false;

  for (auto pb_field_idx = 0, field_idx = 0; pb_field_idx < num_fields;
      ++pb_field_idx, field_idx++) {
    field_name = vid_map_protobuf->infofields(pb_field_idx).name();
    if(m_field_name_to_idx.find(field_name) != m_field_name_to_idx.end()) {
      std::cerr << "Duplicate field name "
                << field_name
                << " found in vid map\n";
      duplicate_fields_exist = true;
      continue;
    }

    // Known fields
    auto known_field_enum = 0u;
    auto is_known_field =
        KnownFieldInfo::get_known_field_enum_for_name(
          field_name,
          known_field_enum);
    // Map
    m_field_name_to_idx[field_name] = field_idx;
    m_field_idx_to_info[field_idx].set_info(field_name, field_idx);
    auto& ref = m_field_idx_to_info[field_idx];

    //Field type - int, char etc
    field_type = vid_map_protobuf->infofields(pb_field_idx).type();
    if (field_type.compare("Integer") == 0) {
      field_type.assign("int");
    } else if (field_type.compare("String") == 0) {
      field_type.assign("char");
    } else if (field_type.compare("Float") == 0) {
      field_type.assign("float");
    } else if (field_type.compare("Flag") == 0) {
      field_type.assign("flag");
    }

    {
      auto iter = VidMapper::m_typename_string_to_type_index.find(field_type);
      VERIFY_OR_THROW(iter != VidMapper::m_typename_string_to_type_index.end()
          && "Field type not handled");
      ref.m_type_index = (*iter).second;
    }
    {
      auto iter = VidMapper::m_typename_string_to_bcf_ht_type.find(field_type);
      VERIFY_OR_THROW(iter != VidMapper::m_typename_string_to_bcf_ht_type.end()
          && "Field type not handled");
      ref.m_bcf_ht_type = (*iter).second;
    }

    // VCF class type can be an array of values: INFO, FORMAT and FILTER
    auto class_type_size =
        vid_map_protobuf->infofields(pb_field_idx).vcf_field_class_size();

    if (class_type_size > 0L) {
      for (int i = 0; i < class_type_size; ++i) {
        std::string class_name =
          vid_map_protobuf->infofields(pb_field_idx).vcf_field_class(i);
        if(class_name == "INFO")
          ref.m_is_vcf_INFO_field = true;
        else if(class_name == "FORMAT") {
          ref.m_is_vcf_FORMAT_field = true;
        } else
          if(class_name == "FILTER")
            ref.m_is_vcf_FILTER_field = true;
      }
    }
    if (vid_map_protobuf->infofields(pb_field_idx).has_length()) {
      std::string length = vid_map_protobuf->infofields(pb_field_idx).length();
      if (isdigit(length.c_str()[0])) {
        ref.m_num_elements =
            strtol(length.c_str(), NULL, 10);
      } else {
        auto iter = VidMapper::m_length_descriptor_string_to_int.find(length);
        if (iter == VidMapper::m_length_descriptor_string_to_int.end()) {
          ref.m_length_descriptor = BCF_VL_VAR;
        } else {
          ref.m_length_descriptor = (*iter).second;
        }
      }
    } else {
      if(is_known_field) {
        auto length_descriptor =
            KnownFieldInfo::get_length_descriptor_for_known_field_enum(
                known_field_enum);
        ref.m_length_descriptor = length_descriptor;
        if(length_descriptor == BCF_VL_FIXED)
          ref.m_num_elements =
              KnownFieldInfo::get_num_elements_for_known_field_enum(
                  known_field_enum,
                  0u,
                  0u);  //don't care about ploidy
      } else {
        ref.m_num_elements = 1;
        ref.m_length_descriptor = BCF_VL_FIXED;
      }
    }

    // Both INFO and FORMAT, throw another entry <field>_FORMAT
    if (ref.m_is_vcf_INFO_field &&
        ref.m_is_vcf_FORMAT_field) {
      auto new_field_idx = field_idx+1u;
      m_field_idx_to_info.resize(m_field_idx_to_info.size()+1u);
      auto& ref = m_field_idx_to_info[field_idx];
      //Copy field information
      m_field_idx_to_info[new_field_idx] = ref;
      auto& new_field_info =  m_field_idx_to_info[new_field_idx];
      //Update name and index - keep the same VCF name
      new_field_info.m_name = field_name+"_FORMAT";
      new_field_info.m_is_vcf_INFO_field = false;
      new_field_info.m_field_idx = new_field_idx;
      new_field_info.m_VCF_field_combine_operation =
          VCFFieldCombineOperationEnum::VCF_FIELD_COMBINE_OPERATION_UNKNOWN_OPERATION;
      //Update map
      m_field_name_to_idx[new_field_info.m_name] = new_field_idx;
      //Set FORMAT to false for original field
      ref.m_is_vcf_FORMAT_field = false;
      field_idx++;
    }
  } // for (auto field_idx = 0; field_idx < num_fields; ++field_idx)

  //Force add END as a field
  auto iter = m_field_name_to_idx.find("END");
  if(iter == m_field_name_to_idx.end()) {
    auto end_idx = m_field_idx_to_info.size();
    m_field_idx_to_info.emplace_back();
    m_field_name_to_idx["END"] = end_idx;
    auto& field_info = m_field_idx_to_info[end_idx];
    field_info.set_info("END", end_idx);
    field_info.m_is_vcf_INFO_field = true;
    field_info.m_type_index = std::move(std::type_index(typeid(int)));
    field_info.m_bcf_ht_type = BCF_HT_INT;
  }
  if(duplicate_fields_exist) {
    throw ProtoBufBasedVidMapperException(
      std::string("Duplicate fields exist in vid map"));
  }

  return GENOMICSDB_VID_MAPPER_SUCCESS;
} // end of parse_infofields_from_vidmap

ProtoBufBasedVidMapperException::~ProtoBufBasedVidMapperException() { ; }

ProtoBufBasedVidMapperException::ProtoBufBasedVidMapperException(
  const std::string m) : msg_("ProtoBufBasedVidMapperException : "+m) { ; }
