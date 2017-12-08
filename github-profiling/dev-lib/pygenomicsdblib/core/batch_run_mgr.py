#! /usr/bin/python3
class BatchManager:
    type_specifics = {
        'loader' : {'f_get_all':'getAllLoaderBatches', 'f_insert':'addBatchLoaderConfig'},
        'query' : {'f_get_all':'getAllQueryBatches', 'f_insert':'addBatchQueryConfig'}} 
    class __BatchItem(object):
        def __init__(self, *args, **kwargs):
            '''lcdef_list, target_cmd, tdb_ws, libs, {'0': tdb_prefix, '1': mpiruns} '''
            self.my_args = args
            self.my_kwargs = kwargs

        def __eq__(self, right):
            if isinstance(right, self.__class__):
                ll = min(len(self.my_args), len(right.my_args))
                return self.my_args[:ll] == right.my_args[:ll] and self.my_kwargs == right.my_kwargs
            else:
                return False
    
    batch_list = {}
    def __init__(self, h_data, test_type):
        if test_type not in self.type_specifics.keys():
            raise RuntimeError("Invalid test type %s" % test_type)
        mysettings = self.type_specifics[test_type]
        assert hasattr(h_data, mysettings['f_get_all']), \
            'invalid method %s' % mysettings['f_get_all']
        assert hasattr(h_data, mysettings['f_insert']), \
            'invalid method %s' % mysettings['f_insert']
        self.__test_type = test_type
        env_libs = h_data.getAdditionalLibs()
        self.__a_libs = {v : k for k, v in env_libs.items()} if env_libs else {} # we want reverse 

        list_generator = getattr(h_data, mysettings['f_get_all'])()
        for batch_id, name, cmd_path, loader_configs, libs, description, varargs in list_generator:
            self.__add2list(batch_id, name, cmd_path, loader_configs, libs, description,  \
             **{'0':varargs[0], '1':varargs[1]})
        self.data_handler = h_data

    def __add2list(self, bid, name, *args, **kwargs):
        batch_item = self.__BatchItem(*args, **kwargs)
        self.batch_list[(bid, name)] = batch_item
        return batch_item
    
    def get_env_id(self, alib_path):
        if alib_path not in self.__a_libs:
            new_id = self.data_handler.addAdditionalLibs(alib_path)
            self.__a_libs[alib_path] = new_id
        else:
            new_id = self.__a_libs[alib_path]
        return new_id

    def __find_batch(self, id_or_name):
        ''' x is id or name '''
        i = 0 if isinstance(id_or_name, int) else 1
        def func(k, n):
            return k[i] == n
        for key, val in self.batch_list.items():
            if func(key, id_or_name):
                return key, val 
        return None, None

    def add_config(self, name, *args, **kwargs):
        ''' check the name
         from loader (name, target_cmd, lcdef_list, tdb_ws, libs, {'0': tdb_prefix, '1': mpiruns})
        '''
        assert len(args) >= 2 and len(kwargs) >= 2, 'invalid arguments: #args=%d, #kwargs=%d' % \
        (len(args), len(kwargs)) 
        batch_id_name = None
        if self.batch_list:
            batch_id_name, known_batch = self.__find_batch(name)
            if batch_id_name:
                new_item = self.__BatchItem(*args, **kwargs)
                if known_batch != new_item:         
                    print('INFO: batch does not match the found one! Replace with current')
        if not batch_id_name:        # new batch
            rt_env_id = self.get_env_id(args[2]) if len(args) == 3 and args[2] else None
            f_add_batch = self.type_specifics[self.__test_type]['f_insert']
            batch_id = getattr(self.data_handler, f_add_batch)(name, args[0], args[1], kwargs['0'], \
             kwargs['1'], rt_env_id)
            self.__add2list(batch_id, name, *args, **kwargs)
            return batch_id, True
        else:
            return batch_id_name[0], False
    
    def __find_my_batch(self, name):
        b_key, b_val = self.__find_batch(name)
        return (b_key, b_val.my_args), b_val.my_kwargs if b_key else None
        
    def find_batch_by_name(self, name):
        return self.__find_my_batch(name)
        
    def find_batch_by_id(self, batch_id):
        return self.__find_my_batch(batch_id)
