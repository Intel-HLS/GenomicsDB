#! /use/python3
# u_result_verifier.py
import os
import os.path
import sys
import json
import pickle
import pandas as pd

TEST_ROOT = '/path/to/logs'
TEST_NAME = 'variant-density'

def __load_json(root, fn):
    with open(os.path.join(root, fn), 'r') as ifd:
        context = json.load(ifd)
    return {fn.split('.')[0] : context}

def __is_query_cfg(fn, fn_type):
    fn_names = fn.split('.')
    return fn.endswith(fn_type) and fn_names[0].split('-')[1].startswith('qc')

def get_query_configs(root_dir):
    ''' return {'1-qc_pos_1_1-1' : qc} '''
    ret = {}
    for root, _, files in os.walk(root_dir, topdown=False):
        _ = [ret.update(__load_json(root, fn)) for fn in files if __is_query_cfg(fn, 'json')]
    return ret

def get_query_outputs(my_root_dir, my_test_name):
    ''' return {'1-qc_pos_1_1-1' : qc} '''
    test_dir = os.path.join(my_root_dir, my_test_name)
    ret = {}
    for root, _, files in os.walk(test_dir, topdown=True):
        for fn in files:
            if __is_query_cfg(fn, 'output'):
                host_name = os.path.basename(os.path.dirname(root))
                if not host_name in ret:
                    ret[host_name] = {}
                ret[host_name].update(__load_json(root, fn))
    return ret

def verify_outputs(outputs):
    ''' outputs is dict of {output_fn : output_json} '''
    keys = list(outputs.keys())
    for k in keys[1:]:
        if outputs[k] != outputs[keys[0]]:
            return  False
    return True

def __get_variant(fvar):
    ''' fvar is variant_calls[i][fields]; returns list of variant values  '''
    if "<NON_REF>" in fvar['ALT']:
        altent = [x for x in  fvar['ALT'] if x != "<NON_REF>"]
        return altent if altent else None
    else:
        raise ValueError('ERROR: invalid field format %s' % fvar)

def build_vcall_df(variant_calls, expand_field=True):
    ''' per output file '''
    cols = ['row', 'int_low', 'int_high']
    col_flds = {'REF', 'ALT', 'GT', 'DP', 'GQ', 'MIN_DP', 'PL'}

    variants = []
    bad_calls = []
    for vcall in variant_calls:
        try:
            altval = __get_variant(vcall["fields"])
            if altval:
                flds1 = {cols[0]: int(vcall['row']), cols[1]: int(vcall['interval'][0]), cols[2]: int(vcall['interval'][1])}
                flds2 = set(vcall['fields'].keys())
                if flds2 and expand_field:
                    extra = flds2 - col_flds
                    if extra:   # remove the extra
                        # print("!! WARN: build_vcall_df got extra fields :", extra)
                        col_flds.update(extra)
                # else:
                #     print("!!!WARN: no fields:", vcall)
                # variants.append({**flds1, **vcall['fields']})
                variants.append(flds1.update(vcall['fields']))
        except Exception as ex:
            print("caught exception,", ex)
            bad_calls.append(vcall)
    df_cols = cols + list(col_flds)
    print("INFO: build_vcall_df, df_col =", df_cols)
    return pd.DataFrame(variants, columns=df_cols), bad_calls

def build_output_df_dict(my_qc_outputs):
    variant_df = {}
    for subtest, out_values in my_qc_outputs.items():
        print("---- %s -----" % subtest)
        df, invalid_calls = build_vcall_df(out_values["variant_calls"])
        if invalid_calls:
            print("WARN: # empty variant calls is ", len(invalid_calls))
        variant_df[subtest] = df        # {'1-qc_1000_1-1': df}
    return variant_df

def __var_range_from_fn(fn):
    num_var_range = fn.split('_')[3].split('-')
    a, b = int(num_var_range[0]), int(num_var_range[1])
    return min(a, b), max(a, b)

def verify_qc_pos2(my_qc_query, my_variant_df):
    ''' not work yet '''
    for subtest, qc in my_qc_query.items():
        if subtest in my_variant_df:
            df = my_variant_df[subtest]
            low_max = df['int_low'].max        #.sort_values(by='int_low') ?
            high_min = df['int_high'].min
            s_cols = pd.Series([int(spos) for spos in qc["query_column_ranges"][0]])
            # s_cols >

def verify_qc_pos(my_qc_query, my_variant_df):
    #TODO:  my_variant_df.sort_values(by='int_low') ?
    for subtest, qc in my_qc_query.items():
        if subtest in my_variant_df:
            df = my_variant_df[subtest]
            expected_size = __var_range_from_fn(subtest)
            counts = [0] * 3        # 0 - pass, 1 - warning, 3 - fail
            for spos in qc["query_column_ranges"][0]:
                pos = int(spos)
                df1 = df[(pos >= df.int_low) & (pos <= df.int_high)]
                # print("!!! pos = %d" % pos)
                # print(df1)
                ret_size = df1.shape[0]
                if ret_size < expected_size[0] or ret_size > expected_size[1]:
                    print("WARN: test %s pos %d, expected %s, got %d" % (subtest, pos, expected_size, ret_size))
                    counts[1] += 1
                else:
                    counts[0] += 1
            print("INFO test %s with #pos=%d, pass %d, pass with glitch %d, failed %d" % (subtest, len(qc["query_column_ranges"][0]), counts[0], counts[1], counts[2]))
        else:
            print('INFO: no output for test %s ' % subtest)

def verify_results(test_root, my_test_name):
    check_point_fn = os.path.join(test_root, my_test_name, "%s_verify.pickle" % my_test_name)
    if not os.path.isfile(check_point_fn):
        qc_outputs = get_query_outputs(test_root, my_test_name)
        # verify outputs on 1 host sincce the others are the same
        assert len(qc_outputs) > 0 and verify_outputs(qc_outputs)
        var_df = build_output_df_dict(qc_outputs[list(qc_outputs.keys())[0]]) # use one host output
        with open(check_point_fn, 'wb') as ofd:
            pickle.dump(var_df, ofd)
        print('INFO saved to ', check_point_fn)
    else:
        print('INFO load from ', check_point_fn)
        with open(check_point_fn, 'rb') as ifd:
            var_df = pickle.load(ifd)

    #var_df.describe()
    verify_qc_pos(get_query_configs(test_root), var_df)

def remove_checkpoint(test_root, my_test_name):
    check_point_fn = os.path.join(test_root, my_test_name, "%s_verify.pickle" % my_test_name)
    if os.path.isfile(check_point_fn):
        os.remove(check_point_fn)
    else:
        print("INFO %s not exist" % check_point_fn)

if __name__ == '__main__':
    test_name = sys.argv[1] if len(sys.argv) > 1 else TEST_NAME
    verify_results(TEST_ROOT, test_name)
