#! /usr/bin/python3
# pylint: disable=broad-except
#
#
import os
import os.path
import platform
from subprocess import call
from pprint import pprint

# jps => vm-id
# jstat -gccause vm-id 5s

TEST_SET = "test-set"
WIN_SET = "win-set"

output_root = os.path.join('mnt', 'app_hdd1', 'scratch', 'test_GATK')

cmd_fmr = "/usr/bin/time -v java -Xms256G -Xmx256G -XX:+UseG1GC -XX:+UseStringDeduplication -jar /home/mingrutar/fromKarthik/gatk-package-distribution-3.7.jar -T CombineGVCFs -L 1:1-155200000 -R /data/broad/samples/joint_variant_calling/broad_reference/Homo_sapiens_assembly19.fasta -o %s -V %s"

_src_dst_dirs = {
    TEST_SET : ("/scratch/generated40k/", [
        ("run_1_10000.combined.gz", 10000),
        ("run_2_100.added.gz", ("run_1_10000.combined.gz", 10000, 100)),
        ("run_3_10100.combined.gz", 10100)]),
    WIN_SET : ("\\Users\\mrutarx\\myprojects\\genomicsDB\\tests\\inputs\\vcfs", [
        ("run_1_10000.combined.gz", 3),
        ("run_2_100.added.gz", ("run_1_10000.combined.gz", 3, 2)),
        ("run_3_10100.combined.gz", 5)])}

def get_src_dir():
    return _src_dst_dirs[TEST_SET] if platform.system() == 'Linux' else _src_dst_dirs[WIN_SET]

def concat_list(first, i_start, num):
    ll = [first]
    ll.extend(src_vcf_files[i_start : i_start + num])
    return ll

def combine_vcfs(vcf_list, out_vcf_fn):
    run_cmd = cmd_fmr % (out_vcf_fn, vcf_list)
    run_cmd = run_cmd.split()
    log_path = "%s.log" % (__file__)
    with open(log_path, 'w') as lfd:
        pprint("CMD = %s" % run_cmd)
        call(run_cmd, shell=False, stdout=lfd, stderr=lfd)

if __name__ == '__main__':
    src_vcfs_dir, run_cfg = get_src_dir()
    src_vcf_files = [os.path.join(src_vcfs_dir, fn) for fn in os.listdir(src_vcfs_dir) if fn[-3:] == ".gz" and not fn.startswith(".claimed.")]

    for out_fn, vcf_i in run_cfg:
        vcfs = concat_list(*vcf_i) if isinstance(vcf_i, tuple) else src_vcf_files[:vcf_i]
        assert vcfs
        list_fn = out_fn[:-2] + "list"
        with open(list_fn, 'w') as ofd:
            _ = [ofd.write("%s\n" % line) for line in vcfs]
        full_out_path = os.path.join(output_root, out_fn)
        combine_vcfs(list_fn, full_out_path)
        print("DONE combined to %s" % out_fn)
