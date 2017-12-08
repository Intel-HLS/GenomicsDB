import sys
import os
import os.path
import pickle
from pprint import pprint
import profiling_env as ppenv

DefaultTestName = 'variant-density'

select_by = {'by_density' : lambda qc, criteria: qc.split('.')[0].split(_)[-1] in criteria,
            'by_num_pos' : lambda qc, criteria: qc.split('_')[-2] in criteria }

#'\\Users\\mrutarx\\myprojects\\GenomicsDBPerfTest\\dev_scratch\\run_mini_tests\\variant-density\\run_list'
def alt_qcfg_test(run_cmd):
    ret = {'run_options' : run_cmd['run_options']}
    ret['run_options']['querier']['produce_opt'] = '-s 10000 ' + run_cmd['run_options']['querier']['produce_opt']
    # all sparse ret['configs'] = [line for line in run_cmd['configs'] if select_by['by_density'](line[2], ['1-1'])]
    # 1 or 10 positions
    ret['configs'] = [line for line in run_cmd['configs'] if select_by['by_num_pos'](line[2], ['1', '10'])]
    return ret

def alt_qcfg(run_cmd, *arg):
    ret = run_cmd.copy()
    ret['run_options']['querier']['produce_opt'] = '-s 10000 ' + run_cmd['run_options']['querier']['produce_opt']
    return ret

alt_func = {'variant-density' : alt_qcfg, 'variant-density_test' : alt_qcfg_test}

if __name__ == '__main__':
    test_name = sys.argv[1] if len(sys.argv) > 1 else DefaultTestName
    fname = os.path.join(ppenv.get_cfg_json_ws(), test_name, 'run_list')
    print('-- test_name=%s, fname=%s' % (test_name, fname))

    with open(fname, 'rb') as ifd:
        context = pickle.load(ifd)
    print('-- original -- has %d #cfg' % len(context['configs']))
    pprint(context)
    alted = alt_func[test_name](context)
    print('-- altered -- has %d #cfg' % len(alted['configs']))
    pprint(alted)
    os.rename(fname, "%s.orig" % fname)
    with open(fname, 'wb') as ofd:
        pickle.dump(alted, ofd)
    print('Alter file %s. Saved the original to %s.orig' % (fname, fname))
