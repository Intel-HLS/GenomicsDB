### GenomicsDB docker builder and  python performance test library

Two dockers are available: importer and querier.

##### How to build dockers
To build importer docker, run:
<code>$ bash docker_src/VCFImporter_builder/build_docker.bash </code>

To build querier docker, run:
<code>$ bash docker_src/GDBQuerier_builder/build_docker.bash </code>

#### How to use the dockers

Assum we have the following directory and files layout:
<code>
ROOT
  |_ mapping_data
        |_ customer_vid.json            # vid mapper file
      [ |_ callsets.json  ]             # callsets file. if not preset, code will generates one
  |_ Homo_sapiens_assembly19.fasta      # reference file
  |_ vcfs
        |_ xxx.gz                       # block compressed VCF file
        |_ xxx.gz.tbi                   # index file
        |_ ....
  |_ genomicsdb_ws                      # where genomicsDB will be

See https://github.com/Intel-HLS/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB for more details on how to make xxx.gz and xxx.gz.tbi from xxx.vcf file.

</code>

For syntax,
<code>
$ docker run vcf_importer:0.6
or
$ docker run genomicsdb_querier:0.9.2-93da4b0-0.5
<code>

<b> example of importing VCF files into GenomicsDB</b>
<code>
$ docker run -v $ROOT:/myenv vcf_importer:0.6 vcf_importer.py -R /myenv/Homo_sapiens_assembly19.fasta -V /myenv/mapping_data/customer_vid.json  -C /myenv/mapping_data/ -o /myenv/genomicsdb_ws/ -i /myenv/vcfs --range 1:1-249250621
</code>

<b> example of querying positions with GenomicsDB</b>
<code>
$ docker run -v $ROOT:/myenv genomicsdb_querier:0.9.2-93da4b0-0.5 genomicsdb_querier.py -R /myenv/Homo_sapiens_assembly19.fasta -V /myenv/mapping_data/customer_vid.json -C /myenv/mapping_data/callsets.json -o /myenv/genomicsdb_ws/ --print-AC --positions /examples/query_1_10_f1.json
</code>

/examples/query_1_10_f1.json is an example query position file. You can find more by

<code>
$ docker run genomicsdb_querier:0.9.2-93da4b0-0.5 ls /examples
</code>

