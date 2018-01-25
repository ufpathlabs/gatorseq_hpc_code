#!/usr/bin/env nextflow


ON_SUCCESS_FILE = params.SAMPLE_DIR + '/' + 'SUCCESS.txt'
ON_FAILURE_FILE = params.SAMPLE_DIR + '/' + 'FAILED.txt'
ON_START_FILE = params.SAMPLE_DIR + '/' + 'START.txt'

onstart()

params.pairs = params.FASTQ_FILES_DIR+"/*_R{1,2}_*.fastq.gz"

Channel
    .fromFilePairs( params.pairs, flat: true)                                     
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .splitFastq(by: 5000000, pe:true, file:true, decompress:true)
    .set {read_pairs}

print read_pairs 

process bwa_map_sort {
lane_num=1
params.dummy_rg_id="'"+'@RG\\tID:'+ params.SAMPLE_NAME + "_" + lane_num + '\\tLB:' + params.RUN_NAME + '\\tSM:' + params.SAMPLE_NAME + '\\tPL:ILLUMINA\\tPU:' +lane_num+"'"
println params.dummy_rg_id


echo true

input:
    set s, r1, r2 from read_pairs
	
output:
	file "${s}.bwa_map_sort.bam" into contigsBam
	file "${s}${params.bwa_map_sort.log}" into bwa_map_sort_log

	"""
    #touch ${s}${params.bwa_map_sort.log}
    #touch ${s}.bwa_map_sort.bam

	bwa mem  -R ${params.dummy_rg_id} \
        -t ${params.bwa_map_sort.threads} \
        -M ${params.human_ref_bwa} \
        ${r1} ${r2} 2> ${s}${params.bwa_map_sort.log} |\
    sambamba view --sam-input \
        --compression-level ${params.bwa_map_sort.compression_level}   \
        --nthreads ${params.bwa_map_sort.sambamba_threads}  \
        --format bam   \
        --output-filename ${s}.bwa_map.bam \
        /dev/stdin 2>> ${s}${params.bwa_map_sort.log}
		
	sambamba sort \
        --nthreads=${params.bwa_map_sort.threads} \
        --memory-limit=${params.bwa_map_sort.memory}  \
        --tmpdir=temp_dir  \
        --out=${s}.bwa_map_sort.bam   \
        --compression-level=${params.bwa_map_sort.compression_level}    \
        ${s}.bwa_map.bam 2>> ${s}${params.bwa_map_sort.log}

    #if [ -d ${params.bwa_map_sort.temp_dir} ]; then rmdir ${params.bwa_map_sort.temp_dir}; fi
	
	"""
}

process merging {
publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.bam*"

input:
  file bamfiles from contigsBam.collect()

output:
  file "${params.merging.dedup_bam}" into dedup
  file "${params.merging.dedup_bam_bai}" into dedupbai
  file "${params.merging.log}" into merging_log
  
"""	
    sambamba merge --nthreads ${params.merging.threads} --compression-level 0 ${params.merging.merge_bam} ${bamfiles} &>${params.merging.log} 
    sambamba markdup --nthreads ${params.merging.threads} --compression-level 9   --tmpdir temp_dir ${params.merging.merge_bam} ${params.merging.dedup_bam} &>>${params.merging.log} 
	
"""

}

process vardict {

publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.vcf"
PADDED_TARGET_BED_FILE=file(params.PADDED_TARGET_BED_FILE)

input:
  file bamfile from dedup
  file baibamfile from dedupbai
  file bedpath from PADDED_TARGET_BED_FILE

output:
  file "${params.vardict.vcf_file}" into outputVcf
  file "${params.vardict.log}" into vardict_log
  
"""	
	
	eval "module load gcc/5.2.0; module load vardictjava/20160521; module load perl/5.20.0; module load R/3.0.2; module load vcftools/0.1.14; module load vt/20160129;"
	
	unset _JAVA_OPTIONS; export _JAVA_OPTIONS="-Xms2g -Xmx20g";
	
	/apps/gcc/5.2.0/vardictjava/20160521/build/install/VarDict/bin/VarDict \
            -G ${params.human_ref_fasta} \
            -f ${params.vardict.AF_THR} -N ${params.SAMPLE_NAME} -b ${bamfile} \
            -c 1 -S 2 -E 3 -g 4 \
            -th 8 -t -r 4 -B 2 -m 6 -P 5 -O 25 -q 25 ${bedpath} 2>> ${params.vardict.log}  | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/teststrandbias.R | \
        /apps/gcc/5.2.0/vardictjava/20160521/VarDict/var2vcf_valid.pl -N ${params.SAMPLE_NAME} -E -f ${params.vardict.AF_THR} \
        > ${params.vardict.vardict_vcf_file}  2>> ${params.vardict.log} 
	
	unset _JAVA_OPTIONS;

    bgzip -f ${params.vardict.vardict_vcf_file} &>> ${params.vardict.log} 
	
	tabix -f ${params.vardict.vardict_vcf_file}.gz &>> ${params.vardict.log} 
	
	vt decompose -s ${params.vardict.vardict_vcf_file}.gz  2>> ${params.vardict.log}  |	vt normalize -r ${params.human_ref_fasta} - 2>> ${params.vardict.log}  | vcf-sort >${params.vardict.vcf_file}  2>> ${params.vardict.log} 
	
"""

}

process generate_metrics_file {

publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.csv"

MERGED_TARGET_BED_FILE= file(params.MERGED_TARGET_BED_FILE) 
PADDED_TARGET_BED_FILE= file(params.PADDED_TARGET_BED_FILE )
PROBE_TARGET_BED_FILE= file(params.TARGET_BED_FILE) 
METRICS_SCRIPT=file(params.generate_metrics_file.METRICS_SCRIPT)

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
  file "${params.generate_metrics_file.log}" into generate_metrics_file_log

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
            ${params.PADDING_SIZE} \
            ${params.MERGED_TARGET_BED_FILE} \
            ${params.PADDED_TARGET_BED_FILE} \
            ${params.generate_metrics_file.mean_stdev_coverage_out_probe_target_bed_file} \
            ${params.generate_metrics_file.sample_id} \
            ${params.generate_metrics_file.sample_name} \
            ${params.generate_metrics_file.run_folder} \
            ${params.human_ref_fasta} \
            ${params.generate_metrics_file.coverage_file} 2>> ${params.generate_metrics_file.log}
"""
}

process merge_log_benchmark_files {

publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.txt"

echo true

input:
  file metrics_file from summary
  file bwa_sort_log from bwa_map_sort_log.collect()
  file merge_log from merging_log.collect()
  file vardict_log from vardict_log.collect()
  file metrics_log from generate_metrics_file_log.collect()


output:
  file "${params.merge_log_benchmark_files.final_logs}" into final_log
  
  
script:
"""	
   echo "MERGED LOG FILES" >${params.merge_log_benchmark_files.final_logs}
        for f in ${bwa_sort_log};
        do
            echo "" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##-------------" \$f "-------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            cat \$f >>${params.merge_log_benchmark_files.final_logs}
        done

        for f in ${merge_log};
        do
            echo "" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##-------------" \$f "-------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            cat \$f >>${params.merge_log_benchmark_files.final_logs}
        done

        for f in ${vardict_log};
        do
            echo "" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##-------------" \$f "-------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            cat \$f >>${params.merge_log_benchmark_files.final_logs}
        done


"""
}

workflow.onComplete { 
	if( workflow.success)
		{
			onsuccess()
		}
	else 
		{
			onfailure()			
		}
   
}

workflow.onError {
    onfailure()
}


def onstart(){
	createFile(ON_START_FILE,"Workflow started at\t")
	println "STARTED"
	}

def onsuccess(){
	createFile(ON_SUCCESS_FILE,"Workflow finished successfully at\t")
	cleanup()
	println "SUCCESS"
}

def onfailure(){
	createFile(ON_FAILURE_FILE,"Workflow failed at\t")
	println "FAILED"
}

def createFile(name,message){
	file = new File(name)
	file.append( "GSBW_VERSION:\t")
	file.append(params.GSBW_VERSION)
	file.append("\n")
	file.append(message)
	file.append("yyyy:mm:dd:hh:mm:ss\t" )
	file.append(new Date().format( 'yyyy:MM:dd:hh:mm:ss' ))
	file.append("\n")
}

def cleanup(){

	println "Deleting Work Directory :" + params.SAMPLE_DIR + "/work"
	result= new File(params.SAMPLE_DIR+"/work").deleteDir()
	if ( result )
		println "Work Directory deleted Successfully"
	else 
		println "Work Directory deletion Failed"
	
	
	println "Deleting Nextflow Cache Directory :" + params.SAMPLE_DIR + "/.nextflow"
	result= new File(params.SAMPLE_DIR+"/.nextflow").deleteDir()
	if ( result )
		println "Nextflow Cache Directory deleted Successfully"
	else 
		println "Nextflow Cache Directory deletion Failed"
		
	println "Deleting Nextflow Log file :" + params.SAMPLE_DIR + "/.nextflow.log"
	
	
	result= new File(params.SAMPLE_DIR+"/.nextflow.log").delete()
	if ( result )
		println "Nextflow Log file deleted Successfully"
	else 
		println "Nextflow Log file deletion Failed"
}