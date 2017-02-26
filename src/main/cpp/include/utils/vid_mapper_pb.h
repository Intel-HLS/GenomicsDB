#ifndef GENOMICSDB_VID_MAPPER_PB_H
#define GENOMICSDB_VID_MAPPER_PB_H

#include "vid_mapper.h"
#include "genomicsdb_vid_mapping.pb.h"
#include "genomicsdb_callsets_mapping.pb.h"

enum {
  GENOMICSDB_VID_MAPPER_SUCCESS = 0x0,
  GENOMICSDB_VID_MAPPER_FAILURE = 0x1
};

class ProtoBufBasedVidMapper : public VidMapper {
  public:

    /**
     * Constructor
     * 1. Callset is required. This is not an option (compared to
     *    FileBasedVidMapper
     * 2. VidMapping protocol buffer structure contains the
     *    merged header from all input GVCFs
     * 3. CallsetMap protocl buffer structure contains the
     *    callset to TileDB row mapping. It also contains
     *    stream names. Stream names cannot be empty
     */
    ProtoBufBasedVidMapper(
      const VidMapping&,
      const CallsetMap&);

    /**
     * Destructor -- must be called once you allocate an object
     */
    ~ProtoBufBasedVidMapper();

    /**
     * Initialization routine to populate this jurassic VidMapper
     * data structure. Super important for maintaining correctness
     * while filling the underlying TileDB array
     */
    void initialize(
      const VidMapping&,
      const CallsetMap&);

    /**
     * Parse the callset map protocol buffer structure and
     * populate data structures of the base jurassic VidMapper
     * class
     */
    int parse_callset_protobuf(const CallsetMap&);

    /**
     * Parse the variant id map protocol buffer structure which
     * contains the merged header. These headers are picked
     * from each of the input GVCF files. Populate the data
     * structures of the base jurassic VidMapper class
     */
    int parse_vidmap_protobuf(const VidMapping& callsetMapProto);
    int parse_contigs_from_vidmap(const VidMapping& vidMapProto);
    int parse_infofields_from_vidmap(const VidMapping& vidMapProto);

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
