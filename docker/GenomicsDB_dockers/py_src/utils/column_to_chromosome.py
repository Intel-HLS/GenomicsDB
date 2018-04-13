import json
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--vid", help="Vid JSON file", type=str, required=True)
parser.add_argument('--query-json', '-q', help='Query json', type=str, required=True)
parser.add_argument('--merge-query-column-ranges', help='Merge query column ranges in to single list',
    action='store_true')
args = parser.parse_args()


def tiledb_column_to_contig_position(contigs_dict, val):
  for contig, contig_info in contigs_dict.iteritems():
    if(val >= contig_info['tiledb_column_offset'] and
        val < contig_info['tiledb_column_offset']+contig_info['length']):
      return (contig, val-contig_info['tiledb_column_offset']+1)
  return None

with open(args.vid, 'rb') as vid_fptr:
  vid_dict = json.load(vid_fptr, object_pairs_hook=OrderedDict)
  contigs_dict = vid_dict['contigs']
  with open(args.query_json, 'rb') as query_fptr:
    query_dict = json.load(query_fptr, object_pairs_hook=OrderedDict)
    input_query_column_ranges = query_dict['query_column_ranges']
    query_column_ranges = []
    curr_output_list = []
    for curr_input_list in input_query_column_ranges:
      curr_input_list.sort()
      if(not args.merge_query_column_ranges):
        curr_output_list = []
      for entry in curr_input_list:
        if(type(entry) is int):
          contig_pos_tuple = tiledb_column_to_contig_position(contigs_dict, entry)
          curr_output_list.append({ contig_pos_tuple[0]: contig_pos_tuple[1]})
        else:
          contig_pos_begin_tuple = tiledb_column_to_contig_position(contigs_dict, entry[0])
          contig_pos_end_tuple = tiledb_column_to_contig_position(contigs_dict, entry[1])
          if(contig_pos_end_tuple[0] != contig_pos_begin_tuple[0]): #passed over contigs
            contig_pos_end_tuple[0] = contig_pos_begin_tuple[0]
            contig_pos_end_tuple[1] = contigs_dict[contig_pos_begin_tuple[0]]['length']
          curr_output_list.append({ contig_pos_begin_tuple[0]: [ contig_pos_begin_tuple[1], contig_pos_end_tuple[1] ] })
      if(not args.merge_query_column_ranges):
        query_column_ranges.append(curr_output_list)
    if(args.merge_query_column_ranges):
      query_column_ranges.append(curr_output_list)
    query_dict['query_column_ranges'] = query_column_ranges
    print(json.dumps(query_dict))
