
class PreloaderInterface(BaseInterface):
    
def preloader(config_fn):
    
    assert os.path.isdir(from_vcf_path), 'Invalid vcf directory: %s' % from_vcf_path
    vcf_list = [os.path.join(from_vcf_path, fn) for fn in os.listdir(from_vcf_path) if fn[-3:] == ".gz" and not fn.startswith('.claimed')]
    assert vcf_list, "no vcf file found in %s" % from_vcf_path
    file_map = generate_vcfs(dest_vcf_dir, num_gen_files, vcf_list)
    
    gen_vcf_list = [fitem[0] for fitem in file_map]
    print("gen_vcf_list=", gen_vcf_list[:5])
    callsets_fname = os.path.join(ws_dir, "callsets_%s.json" % datetime.now().strftime("%y%m%d%H%M"))
    num_part_units = generate_callsets_json(gen_vcf_list, callsets_fname)
