## How to use generate combined VCF file using docker image vcf_combiner
#### Create vcf_combiner docker image

1. clone the project from github.
2. create a directory name tmp
3. copy your built vcf2tiledb to tmp
4. Then RUN:

<pre><code>sudo docker -t vcf_combiner --build-arg GENOMICSDB_BUILD_PATH=/path/to/GenomicsDB</code></pre>

adjust `network setting` as `needed`.

## Command combine_vcf
<p>The command performs VCF files combining is **combine_vcf**.</p>

<b>The usage</b>
<pre><code>combine_vcf -R /container_path/to/ref_file -o /container_path/output/combined.vcf.gz -i /container_path/tn.list</code></pre>

*OR*
 <pre><code>combine_vcf -R /container_path/to/ref_file -o /container_path/output/combined.vcf.gz -i /container_path/vcf_file_1,...,/container_path/vcf_file_n</code></pre>

_where_

     -o : path to outpur file name
     -i : either a list of vcf files or vaf files separated by ','
     -R : reference file
     -c : optional, callset file
*We can combine VCF files at host command line or at docker container command line.*

#### Examples of running at the HOST
>Run with -i a_list_of_vcf_files
<pre><code>sudo docker run -it -v /host_docker_data/:/data/ -v /host_docker_reference/:/ref/ vcf_combiner combine_vcf -R /ref/Homo_sapiens_assembly19.fasta -o /data/output/combined.vcf.gz -i /data/tn.list</code></pre>

>Run with -i vcf_files ...
<pre><code>sudo docker run -it /host_docker_data/:/data/ -v /host_docker_reference/:/ref/ vcf_combiner combine_vcf -R /ref/Homo_sapiens_assembly19.fasta -o /data/output/combined.vcf.gz -i /data/vcfs/info_op1.vcf.gz,/data/vcfs/t7.vcf.gz,/data/vcfs/t1.vcf.gz</code></pre>

#### Examples of running inside the docker container

<p>Create a docker container with bash command
<pre><code>sudo docker run -it -v /host_docker_data/:/data/ -v /host_docker_reference/:/ref/ vcf_combiner bash</code></pre>
</p>
<p>Run with -i a_list_of_vcf_files
<pre><code>combine_vcf -R /ref/Homo_sapiens_assembly19.fasta -o /data/output/combined.vcf.gz -i /data/tn.list</code></pre>
*OR run with -i vcf_files ...*
<pre><code>combine_vcf -R /ref/Homo_sapiens_assembly19.fasta -o /data/output/combined.vcf.2.gz -i /data/vcfs/info_op1.vcf.gz,/data/vcfs/t7.vcf.gz,/data/vcfs/t1.vcf.gz</code></pre>
</p>

##### The look of tn.list of above examples

    /data/vcfs/info_op1.vcf.gz
    /data/vcfs/t7.vcf.gz
    /data/vcfs/t1.vcf.gz
