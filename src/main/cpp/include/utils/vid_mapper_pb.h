#ifndef GENOMICSDB_VID_MAPPER_PB_H
#define GENOMICSDB_VID_MAPPER_PB_H

#include "vid_mapper.h"

class ProtoBufBasedVidMapper : public VidMapper {
  public:
    ProtoBufBasedVidMapper();
    ~ProtoBufBasedVidMapper();

    void initialize();
  protected:
    std::string m_msg;
};

class ProtoBufBasedVidMapperException : public std::exception {
  public:
    /**
     * Default constructor
     */
    ProtoBufBasedVidMapperException(const std::string m="");

    /*
     * Destructor: must be called before end of scope of a
     * VidMapper exception object
     */
    ~ProtoBufBasedVidMapperException();

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

#endif  // GENOMICSDB_VID_MAPPER_PB_H
