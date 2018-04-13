from setuptools import setup, find_packages
# python3 setup.py install --record install.log

setup(name='PyGenomicsDBPerf',
      version='0.2.0.7',
      description='python library for GenomicsDB https://github.com/Intel-HLS/GenomicsDB profiling and docker scripts',
      url='https://github.com/mingrutar/pyPerfProject/tree/dev/py_src',
      download_url = 'https://github.com/mingrutar/pyPerfProject/tree/dev/py_src/archive/0.2.tar.gz',
      author='Ming Rutar',
      author_email='mingx.rutar@intel.com',
      license='MIT',
      packages=find_packages(exclude=['docs', 'tests', '__pycache__']),
      keywords="genomicsdb performance profiling docker",
      # package_data={
      #   'config_jsons': ['*']  },
      # data_files=[('json_cfg', ['config_jsons/16-bn32nt1000000.json'])], not work
      zip_safe=False)
