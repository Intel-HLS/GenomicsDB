# /usr/bin/python3
# pylint: disable=broad-except, arguments-differ

import sys
import os
import os.path
import re
import unittest
import unittest.mock as mock
import logging
import traceback
from importlib import import_module
import test_helper as helper

if __name__ == "__main__":
    myloglevel = logging.INFO
    myformat = '%(levelname)s %(module)s:%(lineno)d: %(message)s' if myloglevel == logging.DEBUG else '%(levelname)s %(module)s: %(message)s'
    logging.basicConfig(stream=sys.stdout, format=myformat)
    logging.getLogger().setLevel(os.environ.get('LOGLEVEL', myloglevel))
logger = logging.getLogger(__name__)

def pre_parse4test(theparser):
    theparser.add_argument('--print-AC', type=bool, nargs='?', const=True, default=True)
    theparser.add_argument('--print-calls', type=bool, nargs='?', const=True, default=False)

class TestParams(unittest.TestCase):
    h_passthrough = None
    @staticmethod
    def get_str(key, val):
        def val2str():
            try:
                ret = ','.join(val)
            except TypeError:
                ret = str(val) if val else ""
            else:
                ret = ','.join(val)
            return ret
        v = val if isinstance(val, str) else val2str()
        return (("-%s %s" if len(key) == 1 else "--%s %s") % (key, v)) if val else ("-%s" if len(key) == 1 else "--%s") % key

    def get_commandline(self, arg_list):
        cmd_str = ' '.join([self.get_str(k, v) for k, v in arg_list.items()])
        # cmd_str = [key2fmt % (k, v if isinstance(v, str) else self.val2str(v)) for k, v in arg_list.items()]
        # cmd_str = ' '.join([("-%s %s" if len(k) == 1 else "--%s %s") % (k, v) for k, v in arg_list.items()])
        return cmd_str

    def check_cmd(self, arg_list, msg=""):
        ''' not final '''
        cmd_str = self.get_commandline(arg_list)
        logger.debug(cmd_str)
        try:
            # myparser.add_argument('--dry-run', type=bool, nargs='?', const=True, default=False)
            with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
                self.target.process(None)
                helper.print_good_usage("Correct command: '%s %s'" % (self.m_target, cmd_str))
        except AssertionError as ex:
            info_text = "%s The wrong input was '%s'" % (msg, cmd_str)
            ex.args = ex.args + (info_text,) if ex.args else (info_text,)
            raise
        except Exception as ex:
            traceback.print_exc(file=sys.stdout)
            info_text = "%s The wrong input was '%s'" % (msg if msg else "Invalid options", cmd_str)
            raise AssertionError(info_text) from ex

    # def _load_param_handler(self, all_args, param_hdl_name):
    #     if '.' in param_hdl_name:
    #         pkg_parts = param_hdl_name.split('.')
    #         mod = import_module('.'.join(pkg_parts[:-1]))
    #         logger.debug("mod=%s, cls=%s", '.'.join(pkg_parts[:-1]), pkg_parts[-1])
    #         self.handler = getattr(mod, pkg_parts[-1])(all_args, helper.get_mock_env())
    #     else:
    #         self.handler = getattr(self.target, param_hdl_name)(all_args, helper.get_mock_env())

    def setUp(self, target_name):
        assert target_name
        self.m_target = target_name
        self.target = import_module(self.m_target)
        all_args = self.target.required_args + self.target.optional_args
        self.all_options = helper.my_options(all_args)           # dict {'C': values}

@unittest.skip("temporarily")
class QuickTest(TestParams):
    def setUp(self):
        pass

#@unittest.skip("temporarily")
class TestImporter(TestParams):
    def setUp(self):
        sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'VCFImporter_builder'))
        TestParams.setUp(self, 'vcf_importer')

    def test_0_all(self):
        self.check_cmd(self.all_options)

    def test_0_vid_predefined(self):
        self.all_options['V'] = helper.get_arg('V', 0)
        self.check_cmd(self.all_options)

    @unittest.skip("somehow predefined not work")
    def test_0_vid_customized(self):
        self.all_options['V'] = helper.get_arg('V', 1)
        self.check_cmd(self.all_options)

    def test_0_range_1(self):
        self.all_options['range'] = helper.get_arg('range', 0)
        self.check_cmd(self.all_options)
    def test_0_range_2(self):
        self.all_options['range'] = helper.get_arg('range', 1)
        self.check_cmd(self.all_options)
    def test_0_range_3(self):
        self.all_options['range'] = helper.get_arg('range', 2)
        self.check_cmd(self.all_options)
    def test_0_range_4(self):
        self.all_options['range'] = helper.get_arg('range', 3)
        self.check_cmd(self.all_options)

    def test_4_miss_required(self):
        miss_opt = re.sub(r'\W+', '', self.target.required_args[0])
        self.all_options.pop(miss_opt)
        with self.assertRaises(AssertionError) as context:
            msg = "miss required option '%s'." % miss_opt
            self.check_cmd(self.all_options, msg)
        helper.print_bad_usage(str(context.exception))

    @unittest.expectedFailure
    def test_4_range_extra_space(self):
        bad_range = "CHR4, CHR11:4-3000"
        self.all_options['range'] = bad_range
        msg = "'range %s', extra space." % bad_range
        cmd_str = self.get_commandline(self.all_options)
        helper.print_bad_usage("%s The wrong input is '%s'" % (msg, cmd_str))
        with self.assertRaises(AssertionError) as context:
            self.check_cmd(self.all_options, "'range %s', extra space." % bad_range)
        helper.print_bad_usage(str(context.exception))

    def test_4_range_end_begin(self):
        bad_range = "CHR1:40000000-3"
        self.all_options['range'] = bad_range
        with self.assertRaises(AssertionError) as context:
            self.check_cmd(self.all_options, "'range %s', the end position is less than the start position." % bad_range)
        helper.print_bad_usage(str(context.exception))

    def test_6_ld_cfg_be(self):
        cmd_str = self.get_commandline(self.all_options)
        cfg_fpath = helper.get_config_fpath('test_6_ld_cfg_be.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        helper.print_good_usage('Generated load config @ %s' % cfg_fpath)
        print('Correctly generated load config @ ', cfg_fpath)
        print()

    def test_6_ld_cfg_chrom(self):
        self.all_options['range'] = helper.get_arg('range', 1)
        cmd_str = self.get_commandline(self.all_options)
        cfg_fpath = helper.get_config_fpath('test_6_ld_cfg_chrom.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        helper.print_good_usage('Generated load config @ %s' % cfg_fpath)
        print('Correctly generated load config @ ', cfg_fpath)
        print()

    def test_8_exec_load_p1(self):
        cmd_str = self.get_commandline(self.all_options)
        cfg_fpath = helper.get_config_fpath('test_8_exec_load_p1.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        with mock.patch('subprocess.call', return_value=123) as mock_call_func:
            retval = self.target.exec_cmd(cfg_fpath)
            assert mock_call_func.called
            assert retval == 123, 'expected 123, got %d' % retval
            args, _ = mock_call_func.call_args
            assert args == (['vcf2tiledb', cfg_fpath],), "expecpect vcf2tildb, got %s" % args[0]
        print('Correctly executed vcf2tiledb with loading 1 partition')
        print()

    @mock.patch('subprocess.call')
    def test_8_exec_load_p4(self, mock_call_func):
        cmd_str = self.get_commandline(self.all_options)
        cfg_fpath = helper.get_config_fpath('test_8_exec_load_p4.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        mock_call_func.return_value = 234
        num_partition = helper.get_number_part(0)
        retval = self.target.exec_cmd(cfg_fpath, num_partition)
        assert mock_call_func.called
        assert retval == 234, 'expected 123, got %d' % retval
        args, _ = mock_call_func.call_args
        logger.info(args)
        assert os.path.basename(args[0][0]) == 'mpirun', "expecpect mpirun, got %s" % os.path.basename(args[0][0])
        assert args[0][1:3] == ['-np', str(num_partition)], "expecpect mpirun, got %s" % args[0][1:3]
        assert args[0][3:] == ['vcf2tiledb', cfg_fpath], "expecpect vcf2tildb, got %s" % args[0][3:]
        print('Correctly executed mpirun vcf2tiledb with loading {} partition'.format(num_partition))
        print()

class TestCallsetsGenerator(TestParams):
    def setUp(self):
        sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'VCFImporter_builder'))
        TestParams.setUp(self, 'callsets_generator')

    def test_gen_callsets(self):
        self.check_cmd(helper.get_inputs4callsets())

#@unittest.skip("temporarily")
class TestQuerier(TestParams):
    def setUp(self):
        sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'GDBQuerier_builder'))
        TestParams.setUp(self, 'genomicsdb_querier')
        self.h_passthrough = pre_parse4test

    def test_0_good(self):
        my_args = helper.my_options(self.target.required_args + self.target.passthrough_args[:1])
        # {'L': helper.get_ldcfg_4query(), 'positions' : helper.get_pos_json4query(), : None}
        self.check_cmd(my_args)

    def test_0_good_attr(self):
        my_args = helper.my_options(self.target.required_args + self.target.optional_args + self.target.passthrough_args[:1])
        self.check_cmd(my_args)

    def test_0_good_calls(self):
        my_args = helper.my_options(self.target.required_args + self.target.passthrough_args[1:])
        # my_args = {'L': helper.get_ldcfg_4query(), 'positions' : helper.get_pos_json4query(), 'print-calls': None}
        self.check_cmd(my_args)

    def test_0_good_AC_calls(self):
        my_args = helper.my_options(self.target.required_args + self.target.passthrough_args[:1])
        my_args['print-calls'] = None
        # my_args = {'L': helper.get_ldcfg_4query(), 'positions' : helper.get_pos_json4query(), 'print-calls': None, 'print-AC': None}
        self.check_cmd(my_args)

    def test_6_config(self):
        my_args = helper.my_options(self.target.required_args + self.target.passthrough_args[:1])
        cmd_str = self.get_commandline(my_args)
        cfg_fpath = helper.get_config_fpath('query_cfg.json')
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        helper.print_good_usage('Generated query config @ %s, cmd=%s' % (cfg_fpath, cmd_str))

    def test_6_config_attr(self):
        my_args = helper.my_options(self.target.required_args + self.target.optional_args + self.target.passthrough_args[:1])
        cmd_str = self.get_commandline(my_args)
        cfg_fpath = helper.get_config_fpath('query_cfg_attr.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        helper.print_good_usage('Generated query config @ %s, cmd=%s' % (cfg_fpath, cmd_str))

    def test_8_exec_query(self):
        my_args = helper.my_options(self.target.required_args + self.target.optional_args + self.target.passthrough_args[:1])
        cmd_str = self.get_commandline(my_args)
        cfg_fpath = helper.get_config_fpath('query_cfg_exec.json')
        try:
            os.remove(cfg_fpath)
        except OSError:
            pass
        with mock.patch.object(sys, 'argv', [self.m_target] + cmd_str.split()):
            _, passthrough_cmds = self.target.process(cfg_fpath)
        assert os.path.exists(cfg_fpath)
        with mock.patch('subprocess.call', return_value=123) as mock_call_func:
            retval = self.target.exec_cmd(cfg_fpath, passthrough_cmds)
            assert mock_call_func.called
            assert retval == 123, 'expected 123, got %d' % retval
            args, _ = mock_call_func.call_args
            assert args[0][0] == 'gt_mpi_gather', "expecpect gt_mpi_gather, got %s" % args[0][0][0]
            assert args[0][1:3] == ['-j', cfg_fpath], "expecpect '-j' %s, got %s" % (cfg_fpath, args[0][1:3])
            for pt in passthrough_cmds:
                assert pt in args[0], 'cannot find the expected %s in args %s' % (pt, args[0])
        helper.print_good_usage('Correctly executed gt_mpi_gather, cmd="%s"' % cmd_str)

if __name__ == '__main__':
    try:
        logger.debug('cwd=%s, helper.TEST_DATA_ROOT=%s', repr(os.getcwd()), repr(helper.TEST_DATA_ROOT))
        unittest.main()
    except Exception as ex:
        print('Error, caught exception ', ex)
        TestImporter.tearDownClass()
        traceback.print_exc(file=sys.stdout)
