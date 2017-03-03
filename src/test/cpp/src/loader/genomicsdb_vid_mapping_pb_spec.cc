
#include "genomicsdb_vid_mapping_pb_spec.h"
#include "vid_mapper_pb.h"
#include "genomicsdb_callsets_mapping.pb.h"

/**
 * Test VidMapper protobuf object read
 */
TEST(VidMappingSpec, test_basic_functions) {

  VidMapping vid_mapping;

  InfoField *pass = vid_mapping.add_infofields();
  pass->set_name("PASS");
  pass->set_type("int");
  pass->add_vcf_field_class_type("INFO");
  pass->set_length("1");

  InfoField *dp = vid_mapping.add_infofields();
  dp->set_name("DP");
  dp->add_vcf_field_class_type("INFO");
  dp->add_vcf_field_class_type("FORMAT");
  dp->set_type("int");
  dp->set_length("G");

  CallsetMap callsetmap;

  SampleIDToTileDBIDMap *sample0 = callsetmap.add_callset_map();
  sample0->set_sample_name("HG00141");
  sample0->set_tiledb_row_index(0);
  sample0->set_sample_vcf_index(0);
  sample0->set_stream_name("HG00141_stream");
  SampleIDToTileDBIDMap *sample1 = callsetmap.add_callset_map();
  sample1->set_sample_name("HG00193");
  sample1->set_tiledb_row_index(1);
  sample1->set_sample_vcf_index(0);
  sample1->set_sample_name("HG00193_stream");

  ProtoBufBasedVidMapper *testVidMap =
      new ProtoBufBasedVidMapper(vid_mapping, callsetmap);

  ASSERT_TRUE(testVidMap!=NULL);
}

