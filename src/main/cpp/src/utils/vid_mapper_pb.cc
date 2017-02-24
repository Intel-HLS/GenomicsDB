
#include "vid_mapper_pb.h"

ProtoBufBasedVidMapper::ProtoBufBasedVidMapper()
  : VidMapper() {

  initialize();
}

void ProtoBufBasedVidMapper::initialize() {
//  m_lb_callset_row_idx = 0;
//  m_ub_callset_row_idx = INT64_MAX-1;
}

ProtoBufBasedVidMapperException::~ProtoBufBasedVidMapperException() { ; }

ProtoBufBasedVidMapperException::ProtoBufBasedVidMapperException(
  const std::string m) : msg_("ProtoBufBasedVidMapperException : "+m) { ; }
