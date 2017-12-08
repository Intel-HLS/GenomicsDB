#! /usr/bin/python3
import json
import os.path
from utils.common import str2num, get_my_logger
from utils.histogram import HistogramManager

one_KB = 1024
one_MB = 1048576

class LoaderConfigMgr(object):
    __extra_data = {}
    my_templates = {}
    tdb_ws = None
    transformer = {
        'String': lambda x: x if isinstance(x, str) else None,
        'Number': str2num,
        'Boolean': lambda x: x if isinstance(x, bool) else x.lower() == 'true',
        'MB': lambda x: int(x) * one_MB,
        'KB': lambda x: int(x) * one_KB}

    user_loader_conf_def = {}       # { uuid: { tag: val} }

    def __init__(self, h_data, wspace):
        self.logger = get_my_logger(self.__class__.__name__)
        self.__getTemplates(h_data, wspace)
        self.loader_tags, self.overridable_tags = h_data.getConfigTags()
        self.defined_loaders = h_data.getAllUserDefinedConfigItems()
        self.data_handler = h_data

    def __getTemplates(self, h_data, wdir):
        template_generator = h_data.getTemplates()
        for name, file_path, params, extra in template_generator:
            temp_name = str(file_path).replace("'", "\"")
            jstrParams = json.loads(str(params).replace("'", "\"")) if params else None
            jstrMore = json.loads(str(extra).replace("'", "\"")) if extra else None
            input_file = self.__proc_template(wdir, temp_name, jstrParams, jstrMore) 
            self.my_templates[name] = input_file if input_file else temp_name

    def __proc_template(self, working_dir, templateFile, sub_key_val=None, extra_args=None):
        f_path = templateFile.replace("$WS_HOME", working_dir)
        if os.path.isfile(f_path):
            if f_path[-5:] == '.temp':
                with open(f_path, 'r') as fd:
                    context = fd.read()
                jf_path = "%s.json" % f_path[:-5] 
                if sub_key_val:
                    for key, val in sub_key_val.items():
                        context = context.replace(key, val)
                    with open(jf_path, 'w') as ofd:
                        ofd.write(context)
                f_path = jf_path
            if extra_args:
                for key, val in extra_args.items():
                    getattr(self, key)(val, working_dir) 
            return f_path
        else:
            self.logger.warn("template file %s not found" % f_path)
            return None

    def histogram(self, val, working_dir):
        self.__extra_data['histogram_file_path'] = val.replace("$WS_HOME", working_dir)

    def __get_loader_name(self, config): 
        ''' list of tuple, make a short name, return unique name '''
        name = []
        as_defaults = []
        for key in sorted(config.keys()):
            val = config[key]
            
            if key in self.overridable_tags and str(val) != self.overridable_tags[key][2]:
                if self.overridable_tags[key][1] == 'Boolean':
                    name.append("%s%s" % (self.overridable_tags[key][3], 't' if val else 'f'))
                elif isinstance(val, list):
                    name.append("%s%s%s" % (self.overridable_tags[key][3], val[0], val[1][:2])) 
                else:
                    name.append("%s%s" % (self.overridable_tags[key][3], val))
            else:
                as_defaults.append(key)
        for x in as_defaults:
            del config[x]
        return ''.join(name)

    def add_user_configs(self, user_defined):
        if isinstance(user_defined, list):
            ret_val = []
            for config in user_defined: #  config is dict,
                lcname = self.__get_loader_name(config)
                if lcname not in self.defined_loaders:
                    loader_id = self.data_handler.addUserDefinedConfig(lcname, str(config))
                    self.defined_loaders[lcname] = (config, loader_id) 
                self.user_loader_conf_def[lcname] = self.defined_loaders[lcname][0]
                ret_val.append(lcname)
            return ret_val
        else:
            print("WARN add requies a list object") 
            return None

    def __get_value(self, itemType, itemVal):
        ''' String, Number, Boolean, Template, MB, KB, func() '''
        if itemType in self.transformer:
            return self.transformer[itemType](itemVal)
        elif itemType == 'Template':
            return self.my_templates[itemVal]
        elif itemType[-2:] == '()':
            return getattr(self, itemType[:-2])(itemVal)
        else:
            return None

    def gen_load_config(self, lc_id, tdb_ws):
        ''' return json load file '''
        self.tdb_ws = tdb_ws
        lc_items = self.user_loader_conf_def[lc_id]
        load_conf = {}
        partition_num = 1
        for key, val in self.loader_tags.items():
            load_conf[key] = self.__get_value(val[1], val[2])
        for key, val in self.overridable_tags.items():
            if key in lc_items:
                uval = lc_items[key]
                if isinstance(uval, list):    # user override the type 
                    load_conf[key] = self.__get_value(uval[1], uval[0])
                else:
                    load_conf[key] = self.__get_value(val[1], uval)
                if key == "column_partitions":
                    partition_num = int(lc_items[key])
            else:                      # use default
                load_conf[key] = self.__get_value(val[1], val[2])
        return load_conf, partition_num

    def get_user_loader_cfg(self):
        return self.user_loader_conf_def
    
    def make_col_partition(self, bin_num):
        bin_num = int(bin_num)
        partitions = []        
        if 'histogram_file_path' in self.__extra_data:
            histogram_fn = self.__extra_data['histogram_file_path']
            hm = HistogramManager(histogram_fn)
            begin_list = hm.calc_bin_begin_pos(bin_num)
            for parnum, begin in enumerate(begin_list):
                partitions.append({"array":"TEST%d" % parnum, "begin": begin, "workspace": self.tdb_ws})
        else:
            # TODO: run histogram utility to generate ?
            print("WARN: no histogram file. 1 partition only")      
            partitions.append({"array":"TEST0", "begin": 0, "workspace": self.tdb_ws})
        return partitions
    