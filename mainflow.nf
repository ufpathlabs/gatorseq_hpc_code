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

echo true

input:
    set s, r1, r2 from read_pairs
	
output:
	file "${s}.bwa_map_sort.bam" into contigsBam

	"""
	
	bwa mem    -t ${params.bwa_map_sort.threads}  -M ${params.human_ref_bwa}  ${r1} ${r2} 2> ${params.bwa_map_sort.log} | sambamba view --sam-input  --compression-level 0   --nthreads 1  --format bam   --output-filename ${s}.bwa_map.bam /dev/stdin 2>> ${params.bwa_map_sort.log}
		
	sambamba sort --nthreads=${params.bwa_map_sort.threads} --memory-limit=${params.bwa_map_sort.memory}  --tmpdir= tempdir   --out=${s}.bwa_map_sort.bam   --compression-level=0    ${s}.bwa_map.bam 2>> ${params.bwa_map_sort.log}
	
	"""
}

process merging {

echo true

PADDED_TARGET_BED_FILE=file(params.PADDED_TARGET_BED_FILE)

input:
  file bamfiles from contigsBam.collect()
  file bedpath from PADDED_TARGET_BED_FILE

output:
  file "${params.merging.merge_bam}" into mapped_reads
  file "${params.merging.dedup_bam}" into dedup
  file "${params.merging.dedup_bam_bai}" into dedupbai
  file "${params.merging.vcf_file}" into outputVcf
  
"""	
    sambamba merge --nthreads ${params.merging.threads} --compression-level 0 ${params.merging.merge_bam} ${bamfiles} &>${params.merging.log} 
	
    sambamba markdup --nthreads ${params.merging.threads} --compression-level 9   --tmpdir temp_dir ${params.merging.merge_bam} ${params.merging.dedup_bam} &>>${params.merging.log} 
	
	eval "module load gcc/5.2.0; module load vardictjava/20160521; module load perl/5.20.0; module load R/3.0.2; module load vcftools/0.1.14; module load vt/20160129;"
	
	unset _JAVA_OPTIONS; export _JAVA_OPTIONS="-Xms2g -Xmx20g";
	
	/apps/gcc/5.2.0/vardictjava/20160521/build/install/VarDict/bin/VarDict \
            -G ${params.human_ref_fasta} \
            -f ${params.AF_THR} -N sample_name -b ${params.merging.dedup_bam} \
            -c 1 -S 2 -E 3 -g 4 \
            -th 8 -t -r 4 -B 2 -m 6 -P 5 -O 25 -q 25 ${bedpath} 2>> ${params.merging.log}  | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/teststrandbias.R | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/var2vcf_valid.pl -N sample_name -E -f ${params.AF_THR} \
        > ${params.merging.vardict_vcf_file}  2>> ${params.merging.log} 
	
	unset _JAVA_OPTIONS;

    bgzip -f sample_variant.vcf &>> ${params.merging.log} 
	
	tabix -f sample_variant.vcf.gz &>> ${params.merging.log} 
	
	vt decompose -s sample_variant.vcf.gz  2>> ${params.merging.log}  |	vt normalize -r ${params.human_ref_fasta} - 2>> ${params.merging.log}  | vcf-sort >${params.merging.vcf_file}  2>> ${params.merging.log} 
	
"""

}

process generate_metrics_file {

echo true

MERGED_TARGET_BED_FILE= file(params.MERGED_TARGET_BED_FILE) 
PADDED_TARGET_BED_FILE= file(params.PADDED_TARGET_BED_FILE )
PROBE_TARGET_BED_FILE= file(params.TARGET_BED_FILE) 
METRICS_SCRIPT=file("metrics.sh")

input:
  file bamfile from dedup
  file vcffile from outputVcf
  file baibamfile from dedupbai
  file merged_target_bed_file from MERGED_TARGET_BED_FILE
  file padded_target_bed_file from PADDED_TARGET_BED_FILE
  file probe_target_bed_file from PROBE_TARGET_BED_FILE
  file metrics_script from METRICS_SCRIPT

output:
  file "${params.generate_metrics_file.metrics_file}" into summary
  file "${params.generate_metrics_file.coverage_file}" into coverage

"""	
    bam stats --in  ${bamfile} --bufferSize ${params.generate_metrics_file.buffer_size} --basic &>${params.generate_metrics_file.total_bamutils_out}
	
	bam stats --in  ${bamfile} --bufferSize ${params.generate_metrics_file.buffer_size}\
	--basic --regionList ${merged_target_bed_file} &> ${params.generate_metrics_file.target_bamutils_out}

	bam stats --in  ${bamfile} --bufferSize ${params.generate_metrics_file.buffer_size} \
		--basic --regionList ${padded_target_bed_file} &>${params.generate_metrics_file.paddedtarget_bamutils_out}
		
	sambamba view -H --nthreads=1  ${bamfile} |\
	grep -P  "^@SQ\tSN:" |cut -f2,3 |sed 's/SN://g' |\
	sed 's/LN://g' >${params.generate_metrics_file.sort_order_chr_names}
	
	bedtools sort -i ${merged_target_bed_file} \
	-faidx ${params.generate_metrics_file.sort_order_chr_names} >${params.generate_metrics_file.resort_merged_target_bed_file}

	bedtools coverage -sorted -d -a ${params.generate_metrics_file.resort_merged_target_bed_file} \
		-b ${bamfile} \
		-g ${params.generate_metrics_file.sort_order_chr_names} >${params.generate_metrics_file.coverage_out_merged_target_bed_file} 2> ${params.generate_metrics_file.log}
		
    bedtools sort -i ${probe_target_bed_file} \
	   -faidx ${params.generate_metrics_file.sort_order_chr_names} >${params.generate_metrics_file.resort_probe_target_bed_file}

    bedtools coverage -sorted -d -a ${params.generate_metrics_file.resort_probe_target_bed_file} \
		-b ${bamfile} \
		-g ${params.generate_metrics_file.sort_order_chr_names} >${params.generate_metrics_file.coverage_out_probe_target_bed_file} 2> ${params.generate_metrics_file.log}

    bedtools groupby -g 1,2,3,4 -c 6 -o mean,stdev -i ${params.generate_metrics_file.coverage_out_probe_target_bed_file} >${params.generate_metrics_file.mean_stdev_coverage_out_probe_target_bed_file} 
	
	bamtools  stats -insert -in ${bamfile} >${params.generate_metrics_file.bamtools_stats_out} 2>> ${params.generate_metrics_file.log}
	
	rtg vcfstats ${vcffile} >${params.generate_metrics_file.rtgtools_stats_out} 2>> ${params.generate_metrics_file.log}
	
	sh ${metrics_script} \
            ${bamfile} \
            ${vcffile} \
            ${params.generate_metrics_file.metrics_file} \
            ${params.generate_metrics_file.total_bamutils_out} \
            ${params.generate_metrics_file.target_bamutils_out} \
            ${params.generate_metrics_file.paddedtarget_bamutils_out} \
            ${params.generate_metrics_file.coverage_out_merged_target_bed_file} \
            ${params.generate_metrics_file.bamtools_stats_out} \
            ${params.generate_metrics_file.rtgtools_stats_out} \
            ${params.generate_metrics_file.padding_size} \
            ${params.generate_metrics_file.merged_target_bed_file} \
            ${params.generate_metrics_file.padded_target_bed_file} \
            ${params.generate_metrics_file.mean_stdev_coverage_out_probe_target_bed_file} \
            ${params.generate_metrics_file.sample_id} \
            ${params.generate_metrics_file.sample_name} \
            ${params.generate_metrics_file.run_folder} \
            ${params.human_ref_fasta} \
            ${params.generate_metrics_file.coverage_file} 2>> ${params.generate_metrics_file.log}
"""
}


def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("_R"))
}