#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

yaml_parameter = new FileInputStream(new File(params.yamlConfig))
new Yaml().load(yaml_parameter).each { k, v -> params[k] = v }

  reads1 = Channel
      .fromPath( params.pair1 )
      .map { path -> tuple(sample(path), path) }


  reads2 = Channel
      .fromPath( params.pair2 )
      .map { path -> tuple(sample(path), path) }

  alignmentReadPairs = Channel.create()
  read_pairs = reads1
          .phase(reads2)
          .map { pair1, pair2 -> [pair1[0], pair1[1], pair2[1] ] }.tap(alignmentReadPairs)

process bwa_map_sort {
    cpus 4     
    echo true

input:
    set s, r1, r2 from read_pairs
output:
	file "${s}.bwa_map_sort.bam" into contigsBam

	"""
	bwa mem    -t 4  -M ${params.human_ref_bwa}  ${r1} ${r2} | sambamba view --sam-input  --compression-level 0   --nthreads 1  --format bam   --output-filename ${s}.bwa_map.bam /dev/stdin
		
	sambamba sort --nthreads=1 --memory-limit=8G  --tmpdir= tempdir   --out=${s}.bwa_map_sort.bam   --compression-level=0    ${s}.bwa_map.bam
	
	"""
}

process merging {
cpus 4
echo true

input:
  file bamfiles from contigsBam.collect()
  file  'PADDED_TARGET_BED_FILE' from params.PADDED_TARGET_BED_FILE

output:
  file "${params.merge_bam}" into mapped_reads
  file 'dedup.bam' into dedup
  file 'sample.vcf' into outputVcf
  
"""	
    sambamba merge --nthreads 8 --compression-level 0 ${params.merge_bam} ${bamfiles} 
	
    sambamba markdup --nthreads 8 --compression-level 9   --tmpdir temp_dir ${params.merge_bam} dedup.bam 
	
	eval "module load gcc/5.2.0; module load vardictjava/20160521; module load perl/5.20.0; module load R/3.0.2; module load vcftools/0.1.14; module load vt/20160129;"
	
	unset _JAVA_OPTIONS; export _JAVA_OPTIONS="-Xms2g -Xmx20g";
	
	/apps/gcc/5.2.0/vardictjava/20160521/build/install/VarDict/bin/VarDict \
            -G ${params.human_ref_fasta} \
            -f ${params.AF_THR} -N sample_name -b dedup.bam \
            -c 1 -S 2 -E 3 -g 4 \
            -th {threads} -t -r 4 -B 2 -m 6 -P 5 -O 25 -q 25 ${PADDED_TARGET_BED_FILE} | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/teststrandbias.R | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/var2vcf_valid.pl -N sample_name -E -f ${params.AF_THR} \
        > sample_variant.vcf
	
	     unset _JAVA_OPTIONS;

    bgzip -f sample_variant.vcf
	
	tabix -f sample_variant.vcf.gz 
	
	vt decompose -s sample_variant.vcf.gz |	vt normalize -r ${params.human_ref_fasta} - | vcf-sort >sample.vcf 
"""
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("_R"))
}