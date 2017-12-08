#! /usr/bin/python3
import os
import os.path
import platform
import pandas as pd

class IntPosBuilder:
    DefaultIntPosRoot = "/path/to/interesting_pos"

    def __init__(self, path_root=None):
        self.interesting_pos_toproot = path_root if path_root else self.DefaultIntPosRoot
        self.df_intpos = self.build_intpos_info(self.interesting_pos_toproot)

    @staticmethod
    def build_intpos_info(topdir):
        sample_dirs = [(sd, sd.split('sample')[1]) for sd in os.listdir(topdir) if os.path.isdir(os.path.join(topdir, sd)) and sd.startswith('sample')]
        tlist = []
        for item in sample_dirs:     # (dirname, num_sample)
            n_sample = int(item[1][:-3])
            iposfiles = os.listdir(os.path.join(topdir, item[0]))
            iposfnlist = sorted(iposfiles, key=lambda x: int(x.split('_')[-1]))
            tlist.append([n_sample, item[0], iposfnlist])
        df_intpos = pd.DataFrame(tlist, columns=["n_sample", "dname", 'files'])
        df_intpos['max_partition'] = df_intpos.files.apply(lambda x: len(x))
        df_intpos['combined'] = df_intpos.apply(lambda r: 'merged_%sp_%s.intpos' % (r['max_partition'], r['n_sample']), axis=1)
        print(df_intpos)
        return df_intpos

    def build_intpos(self, num_sample, num_parallel):
        selected = self.df_intpos[(self.df_intpos.n_sample <= num_sample)].sort_values(by='n_sample', ascending=False).head(1)
        assert selected.iloc[0]['max_partition'] >= num_parallel
        ret_fn = os.path.join(self.interesting_pos_toproot, selected.iloc[0]['combined'])
        if os.path.isfile(ret_fn):
            print("INFO found %s for num_sample=%d, num_parallel=%d" % (ret_fn, num_sample, num_parallel))
            return ret_fn

        try:
            with open(ret_fn, "w") as outfile:
                num_line = 0
                for i in range(selected.iloc[0]['max_partition']) :
                    ipfn = os.path.join(self.interesting_pos_toproot, selected.iloc[0]['dname'], "interesting_pos_%d" % i)
                    with open(ipfn, "r") as infile:
                        for aline in infile:
                            try:
                                _ = [int(token) for token in aline.split()]
                                num_line += 1
                                outfile.write(aline)
                            except ValueError:
                                print("INFO skip ", aline)
            print("INFO made %s contains %d lines for num_sample=%d, num_parallel=%d" % (ret_fn, num_line, num_sample, num_parallel))
        except Exception as ex:
            print("Failed: removed file %s, exception=%s", (ret_fn, ex))
            if os.path.isfile(ret_fn):
                os.remove(ret_fn)
            return None
        return ret_fn

def make_files(top):
    if platform.system() == 'Windows':
        sample_dir = ["sample10000_28", "sample7500_28", "sample5000_26", "sample3000_26", "sample1000_26"]
        for ss in sample_dir:
            dname = os.path.join(top, ss)
            if not os.path.isdir(dname):
                os.makedirs(dname)
            for i in range(16):
                ip_fn = os.path.join(dname, 'interesting_pos_%d'%i)
                with open(ip_fn, 'w') as ofd:
                    ofd.write("I am %s" % ip_fn)
                print('make_files: ip_fn =', ip_fn)

if __name__ == "__main__":
    # make_files(IntPosBuilder.DefaultIntPosRoot)
    ipbuiler = IntPosBuilder()
    print('run ipbuiler.build_intpos(1000, 16)...')
    out_fn = ipbuiler.build_intpos(1000, 3)
    print('Done ipbuiler.build_intpos(1000, 3), out_fn=', out_fn)
    print('-------')
    # print('run ipbuiler.build_intpos(5500, 7)...')
    # out_fn = ipbuiler.build_intpos(5500, 7)
    # print('Done ipbuiler.build_intpos(5500, 7), out_fn=', out_fn)

    # print('run ipbuiler.build_intpos(7400, 9)...')
    # out_fn = ipbuiler.build_intpos(7400, 9)
    # print('Done ipbuiler.build_intpos(7400, 9), out_fn=', out_fn)
