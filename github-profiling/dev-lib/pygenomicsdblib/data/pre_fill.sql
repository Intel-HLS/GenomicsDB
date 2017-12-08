--- pre_fill.sql

--- clean all tables, sqlite does not have truncate
DELETE FROM host;
DELETE FROM template;
DELETE FROM loader_config_tag;

--- fill host TABLE, KD@compute-2-26; docker@compute-2-23 
INSERT INTO host (hostname, availibility) VALUES ("compute-2-22", 0);
INSERT INTO host (hostname, availibility) VALUES ("compute-2-23", 0);
INSERT INTO host (hostname, availibility) VALUES ("compute-2-24", 0);
INSERT INTO host (hostname, availibility) VALUES ("compute-2-25", 0);
INSERT INTO host (hostname) VALUES ("compute-2-26");
INSERT INTO host (hostname) VALUES ("compute-2-27");
INSERT INTO host (hostname) VALUES ("compute-2-28");
INSERT INTO host (hostname) VALUES ("compute-2-29");

---scheme_info
INSERT INTO scheme_info (version) VALUES ("0.2");

--- template 
INSERT INTO template (name, my_type, file_path) VALUES ("broad_vid", 1 '$WS_HOME/templates/vid.json' );
INSERT INTO template (name, my_type, file_path, params) 
VALUES ("callsets_1K", 2, '$WS_HOME/templates/callsets1000.json', '{"histogram": "@histogram_1K@"}');
INSERT INTO template (name, my_type, file_path, params) 
VALUES ("callsets_gen40K", 2, '$WS_HOME/templates/callsets_gen40000.json','{"histogram": "@histogram_1K@"}' );
INSERT INTO template (name, my_type, file_path) 
VALUES ("vcf_header", 3, '$WS_HOME/templates/template_vcf_header.vcf');
INSERT INTO template (name, my_type, file_path) 
VALUES ("ref_genome", 4, '/data/broad/samples/joint_variant_calling/broad_reference/Homo_sapiens_assembly19.fasta' );
INSERT INTO template (name, my_type, file_path) VALUES ("histogram_1K", 5 '$WS_HOME/templates/1000_histogram');
INSERT INTO template (name, my_type, file_path, params) 
VALUES ("loader_tiledb", 6, "$WS_HOME/templates/loader_cfg_tiledb.template", "{'alter': 'partition.templet'}");
INSERT INTO template (name, my_type, file_path, params) 
VALUES ("loader_tiledb", 6, "$WS_HOME/templates/loader_cfg_vcf.template", "{'alter': 'partition.templet'}");
INSERT INTO template (name, my_type, file_path) 
VALUES ("loader_tiledb", 7, "$WS_HOME/templates/query_cfg_vcf.template");

--- loader_config_tag 
INSERT INTO loader_config_tag (name, my_type, default_value) 
VALUES ("row_based_partitioning", 'Boolean', "false");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("produce_tiledb_array", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("treat_deletions_as_intervals", 'Boolean', "false");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("vcf_header_filename", 'Template', "vcf_header");
INSERT INTO loader_config_tag (name, type, default_value, for_combine) 
VALUES ("reference_genome", 'Template', "ref_genome", 1);
INSERT INTO loader_config_tag (name, type, default_value, for_combine) 
VALUES ("offload_vcf_output_processing", 'Boolean', "true", 1);
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("ub_callset_row_idx", 'Number', 999);
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("discard_vcf_index", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("produce_tiledb_array", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("disable_synced_writes", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("delete_and_create_tiledb_array", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("ignore_cells_not_in_partition", 'Boolean', "false");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("tiledb_compression_level", 'Number', "6");

-- for produce_combined_vcf = true ??
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("produce_combined_vcf", 'Boolean', "true");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("vcf_output_format", 'String', "z");
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("produce_GT_field", 'Boolean', "false");

-- unit MiB, default 100 MiB
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("size_per_column_partition", '()', "lambda x: (x.ub_callset_row_idx-x.lb_callset_row_idx+1)*16384", "ul"); 
-- default compress. 
INSERT INTO loader_config_tag (name, type, default_value) 
VALUES ("compress_tiledb_array", 'Boolean', "true");

-- user definable columns 
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("column_partitions", "make_col_partition()", '1', 'bn'); 
-- number parallel read to 1, so far 2 didn't improve' 
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("num_parallel_vcf_files", 'Number', "1", 'pf');
-- number cells per tile
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("num_cells_per_tile", 'Number', 1000, 'nt');
---added after feb 17
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("segment_size", 'MB', "10", 'sg');
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("do_ping_pong_buffering", 'Boolean', "false", 'pb');
---added may 17
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("callset_mapping_file", 'Template', "callsets", 'cm');
INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("vid_mapping_file", 'Template', "vid", "vd");

INSERT INTO loader_config_tag (name, type, default_value, tag_code) 
VALUES ("vcf_output_filename", 'String', "combined.vcf.gz", "vo");

-- for query_config_tag, not in use yet --
INSERT INTO query_config_tag (name, type, default_value, tag_code) 
VALUES ("segment_size", 'MB', 10, "qsg");
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("num_position", "Number", 0);
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("scan_full", "Boolean", 'true');
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("array", "String", 'TEST0');
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("query_attributes", "String", "REF,ALT,BaseQRankSum,MQ,MQ0,ClippingRankSum,MQRankSum,ReadPosRankSum,DP,GT,GQ,SB,AD,PL,DP_FORMAT,MIN_DP");

-- partition --
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("chromosome", "Number", 1);
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("begin", "Number", 0);
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("end", "Number", -1);
INSERT INTO query_config_tag (name, type, default_value) 
VALUES ("workspace", "String", '');
