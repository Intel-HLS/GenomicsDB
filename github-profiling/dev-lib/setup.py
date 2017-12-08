from setuptools import setup, find_packages
# python3 setup.py install --record install.log

setup(name='PyGenomicsdbLib',
      version='0.6.0.13',
      description='python library for GenomicsDB https://github.com/Intel-HLS/GenomicsDB profiling',
      url='https://github.com/mingrutar/GenomicsLoaderConfigCreator',
      author='Ming Rutar',
      author_email='mingx.rutar@intel.com',
      license='MIT',
      packages=find_packages(exclude=['docs', 'tests', 'run_tests', '__pycache__', 'cfg_jsons']),
      keywords="genomicsdb performance profiling",
      # package_data={
      #   'config_jsons': ['*']  },
      # data_files=[('json_cfg', ['config_jsons/16-bn32nt1000000.json'])], not work
      zip_safe=False)