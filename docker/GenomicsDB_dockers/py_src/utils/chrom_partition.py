#! python3

import os
import os.path
import json
import argparse
import itertools
from collections import OrderedDict

def to_json(by_chromosome_dict, fn_prefix):
    for k, v in by_chromosome_dict.items():
        outfn = "%s_%s.json" % (fn_prefix, k)
        with open(outfn, 'w') as ofp:
            json.dump(v, ofp)
        print("created '%s' for chromosome '%s'" % (outfn, k))

def process_be(fn, vid_fn, chrom_list=None):
    ''' "query_column_ranges": [[188920061, 11298738, 32065060, 84081239], [...]] '''
    assert vid_fn
    with open(vid_fn, 'r') as ifp:
        vid_json = json.load(ifp)
    vidmap = sorted([(vals["length"] + vals["tiledb_column_offset"], vals["tiledb_column_offset"]+1, chrname) for chrname, vals in vid_json['contigs'].items()], key=lambda t: t[0])

    with open(fn, 'r') as ifp:
        data = json.load(ifp)
    irange = list(itertools.chain(*data["query_column_ranges"]))
    by_chrom = {k : [v] for k, v in {chrom : [x for x in irange if x >= lval and x <= hval] for hval, lval, chrom in vidmap if chrom in chrom_list}.items() if v}
    to_json(by_chrom, os.path.splitext(fn)[0])

def process_chr_be(fn, chromosome_list=None):
    '''"query_column_ranges": [ [{"1": [849898, 854972]}, {"1": [876528, 885685]},...], [...,{"X": [153746665, 153760508]}] ]'''
    with open(fn, 'r') as ifp:
        data = json.load(ifp)
    by_chromosome_dict = {x: [[]] for x in chromosome_list} if chromosome_list else {}
    for ll1 in data["query_column_ranges"]:
        for di in ll1:
            keychr = list(di.keys())[0]
            if keychr not in by_chromosome_dict and not chromosome_list:
                by_chromosome_dict[keychr] = [[]]
            if keychr in by_chromosome_dict:
                # print('match: ', keychr, di[keychr])
                if not isinstance(di[keychr], list):
                    by_chromosome_dict[keychr][0].append({keychr : di[keychr]})
                else:
                    ddict = OrderedDict()
                    ddict[keychr] = di[keychr]
                    dddict = dict(ddict)
                    by_chromosome_dict[keychr][0].append(dddict)

    to_json(by_chromosome_dict, os.path.splitext(fn)[0])

if __name__ == '__main__':
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-f', type=str, required=True)
    my_parser.add_argument('-c', type=str, required=False)
    my_parser.add_argument('-v', type=str, required=False)
    args = my_parser.parse_args()
    chrlist = args.c.split(',') if args.c else None
    if args.v:
        process_be(args.f, args.v, chrlist)
    else:
        process_chr_be(args.f, chrlist)
    print("Done")
