#! /usr/bin/python3
#
# create 1000g_interesting_positions.var that only contains variant rows from KG's ~/fromKarthik/1000g_interesting_positions
import genopp_process_utils

src = "/path/to/1000g_interesting_positions"
dst = "/path/to/1000g_interesting_positions.var"

df = genopp_process_utils.get_variant_df(src)
df2 = df.drop('num_variant', 1)
df2.to_csv(dst, index=False, header=None, sep=' ')
