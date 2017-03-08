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
      new ProtoBufBasedVidMapper(&vid_mapping, &callsetmap);

  ASSERT_TRUE(testVidMap!=NULL);
}

