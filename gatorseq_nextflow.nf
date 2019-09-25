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
    .map{b->[sample(b[1]),b[1],b[2]]}
    .set {read_pairs}

def sample(Path path){
    def name = path.getFileName().toString()
    int start = Math.max(0, name.lastIndexOf('/'))
    return name.substring(start)
}

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
	
    #https://gitter.im/nextflow-io/nextflow/archives/2018/01/30
    #https://gitter.im/nextflow-io/nextflow/archives/2017/06/02

    which VarDict 2> ${params.vardict.log}
    VarDict \
        -G ${params.human_ref_fasta} \
        -f ${params.vardict.AF_THR} -N ${params.SAMPLE_NAME} -b ${bamfile} \
        -c 1 -S 2 -E 3 -g 4 \
        -th 8 -t -r 4 -B 2 -m 6 -P 5 -O 25 -q 25 ${bedpath} \
        > ${params.vardict.VarDict_out_file} 2>> ${params.vardict.log}  

    which teststrandbias.R 2>> ${params.vardict.log}

    cat ${params.vardict.VarDict_out_file} | \
        teststrandbias.R \
        > ${params.vardict.teststrandbias_out_file} 2>> ${params.vardict.log} 

    which var2vcf_valid.pl 2>> ${params.vardict.log}
    cat ${params.vardict.teststrandbias_out_file} | \
            var2vcf_valid.pl \
            -N ${params.SAMPLE_NAME} \
            -E -f ${params.vardict.AF_THR} \
        > ${params.vardict.vardict_vcf_file}  2>> ${params.vardict.log} 
	

    bgzip -f ${params.vardict.vardict_vcf_file} &>> ${params.vardict.log} 
	tabix -f ${params.vardict.vardict_vcf_file}.gz &>> ${params.vardict.log} 

	
    which vt 2>> ${params.vardict.log}
	vt decompose -s ${params.vardict.vardict_vcf_file}.gz  2>> ${params.vardict.log}  |	\
    vt normalize -r ${params.human_ref_fasta} - 2>> ${params.vardict.log}  | \
    vcf-sort >${params.vardict.vcf_file}  2>> ${params.vardict.log} 

    sed -i '1 a ##reference=NCBIb37' ${params.vardict.vcf_file}

    
"""

}


process vepAnnot {
publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.vep.vcf.txt"

input:
  file varDictVcf from outputVcf 

output:
  file "${params.vepAnnot.vep_vcf_file}" into outputVepVcf
  file "${params.vepAnnot.log}" into vepAnnot_log
  
"""	
    echo "TEMP PLACEHOLDER FOR VEP" &>> ${params.vepAnnot.log}
    echo "TEMP PLACEHOLDER FOR VEP" &>> ${params.vepAnnot.vep_vcf_file}
    
	##COMMENT #### vt decompose -s ${params.vardict.vardict_vcf_file}.gz -o ${params.vardict.decompose_vcf_file} &>> ${params.vardict.log}
    #normVCF -o ${params.vepAnnot.norm_vcf_file} --reference ${params.human_ref_fasta} --sample ${varDictVcf} &>> ${params.vepAnnot.log}

    ##transvar
    #transvar ganno --vcf ${params.vepAnnot.norm_vcf_file} --reference ${params.TRANSVAR_REF_GENOME} > ${params.vepAnnot.transvar_temp_file} 2>> ${params.vepAnnot.log}

    #cat ${params.vepAnnot.transvar_temp_file} |grep  "^##" > ${params.vepAnnot.transvar_vcf_file} 2>> ${params.vepAnnot.log}
    #cat ${params.vepAnnot.transvar_temp_file} |grep  "^#CHR" | cut -f1-10 >> ${params.vepAnnot.transvar_vcf_file} 2>> ${params.vepAnnot.log}
    #cat ${params.vepAnnot.transvar_temp_file} |grep -v "^#" | cut -f1-10,14,16 | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 "\t" \$6 "\t" \$7 "\t" \$8 ";COORDINATES=" \$11 ";" \$12 "\t" \$9 "\t" \$10}' 2>> ${params.vepAnnot.log} |sort -k1,1 -k1,1n >> ${params.vepAnnot.transvar_vcf_file} 2>> ${params.vepAnnot.log}

    ##module load vep/88.8
    #${params.vepAnnot.vep_program} \
    #    --species ${params.VEP_SPECIES} \
    #    --assembly ${params.VEP_ASSEMBLY} \
    #    --input_file ${params.vepAnnot.transvar_vcf_file} \
    #    --format vcf \
    #    --output_file ${params.vepAnnot.vep_vcf_file} \
    #    --force_overwrite \
    #    --stats_file ${params.vepAnnot.vep_stat_file} \
    #    --cache \
    #    --dir ${params.VEP_FOLDER} \
    #    --dir_cache ${params.VEP_FOLDER} \
    #    --dir_plugins ${params.VEP_FOLDER} \
    #    --offline \
    #    --fasta ${params.VEP_FASTA} \
    #    --merged  \
    #    --exclude_predicted \
    #    --cache_version ${params.VEP_CACHE_VERSION} \
    #    --everything \
    #    --hgvsg \
    #    --flag_pick \
    #    --flag_pick_allele \
    #    --flag_pick_allele_gene \
    #    --xref_refseq \
    #    --check_existing \
    #    --fork ${params.VEP_FORK} \
    #    --vcf --vcf_info_field ANN 2>> ${params.vepAnnot.log}


"""	
}


process iCallSV {

//publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}*"

input:
  file bamfile from dedup
  file baibamfile from dedupbai
  file dummyOutputVepVcf from outputVepVcf 

output:
  file "${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_only_final.txt" into outputICallSV
  file "${params.iCallSV.log}" into iCallSV_log
  
"""	
mkdir ${params.iCallSV.outDir} 
cat  ${params.iCallSV.svConfig} | sed "s|PATH_FOR_REFRENCE_GENOME|${params.human_ref_fasta}|g" | sed "s|ICALLSV_RESOURCES|${params.icallsv_resources_folder}|g" >${params.iCallSV.iCallSV_Config}

python ${params.iCallSV.icallsv_program} \
    -sc ${params.iCallSV.iCallSV_Config} \
    --caseBam ${bamfile} \
    --controlBam ${params.iCallSV.controlBAMFile} \
    --caseId ${params.SAMPLE_NAME} \
    --controlId ${params.iCallSV.controlID} \
    --outDir ${params.iCallSV.outDir} \
    --outPrefix ${params.iCallSV.outPrefix} 2>> ${params.iCallSV.log}

cat ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_final.txt | head -n 1 >${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_only_final.txt 2>> ${params.iCallSV.log}
cat ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_final.txt | grep -P '\tTRA\t' |cat >>${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_only_final.txt 2>> ${params.iCallSV.log}
cat ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_final.txt | grep -P '\tINV\t' |cat >>${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_only_final.txt 2>> ${params.iCallSV.log}

cp ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_merged.txt  ${params.SAMPLE_DIR}/${params.iCallSV.outPrefix}_merged.xls 
cp ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_final.txt  ${params.SAMPLE_DIR}/${params.iCallSV.outPrefix}_final.xls
cp ${params.iCallSV.DellyDir}/${params.iCallSV.outPrefix}_only_final.txt  ${params.SAMPLE_DIR}/${params.iCallSV.outPrefix}_only_final.xls

"""	
}


process mSINGS {

publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.mSINGS.txt"

input:
  file bamfile from dedup
  file baibamfile from dedupbai

output:
  file "${params.merging.msi_analyzer_output}" into outputMSINGS 
  
"""	
	
    current_dir=$(pwd)
    echo "Starting MSI Analysis of ${params.SAMPLE_NAME}" > ${params.mSINGS.log};

    #echo "sorting bam" >> ${params.mSINGS.log};
    #date +"%D %H:%M" >> ${params.mSINGS.log};
    #samtools sort ${bamfile} ${params.mSINGS.sorted_bam} && samtools index ${params.mSINGS.sorted_bam}.bam &>> ${params.mSINGS.log};

    #echo "Making mpileups" >> ${params.mSINGS.log};
    #date +"%D %H:%M" >> ${params.mSINGS.log};
    #samtools mpileup -f ${params.human_ref_fasta} -d 100000 -A -E  ${params.mSINGS.sorted_bam}.bam -l ${params.mSINGS.intervals_file} | awk '{if($4 >= 6) print $0}' > ${params.mSINGS.mpileup_file} &>> ${params.mSINGS.log};

    
    #echo "Varscan Readcounts start" >> ${params.mSINGS.log};
    #date +"%D %H:%M" >> ${params.mSINGS.log};
    #java -Xmx4g -jar VarScan.v2.3.7.jar readcounts ${params.mSINGS.mpileup_file} --variants-file ${params.mSINGS.intervals_file} --min-base-qual 10 --output-file ${params.mSINGS.varscan_readcounts_file}  &>> ${params.mSINGS.log}; &
    #wait

    #echo "MSI Analyzer start" >> ${params.mSINGS.log};
    #date +"%D %H:%M" >> ${params.mSINGS.log};
    #msi analyzer ${params.mSINGS.varscan_readcounts_file} ${params.mSINGS.bed_file} -o ${params.mSINGS.msi_analyzer_output} &>> ${params.mSINGS.log};

    #echo "MSI calls start" >> ${params.mSINGS.log};
    #date +"%D %H:%M" >> ${params.mSINGS.log};
    #msi count_msi_samples ${params.mSINGS.msi_baseline} $current_dir -m ${params.mSINGS.multiplier} -t ${params.mSINGS.msi_min_threshold} ${params.mSINGS.msi_max_threshold} -o ${params.SAMPLE_NAME}.MSI_Analysis.txt &>> ${params.mSINGS.log};

    touch ${params.mSINGS.msi_analyzer_output}
    echo "Completed MSI Analysis" >> ${params.mSINGS.log};
    date +"%D %H:%M" >> ${params.mSINGS.log};

"""

}


process generate_metrics_file {

publishDir params.SAMPLE_DIR, mode: 'copy', overwrite: true,pattern: "*.csv"

MERGED_TARGET_BED_FILE= file(params.MERGED_TARGET_BED_FILE) 
PADDED_TARGET_BED_FILE= file(params.PADDED_TARGET_BED_FILE )
PROBE_TARGET_BED_FILE= file(params.TARGET_BED_FILE) 
METRICS_SCRIPT=file(params.generate_metrics_file.METRICS_SCRIPT)

input:
  file icalFile from outputICallSV 
  file bamfile from dedup
  file vcffile from outputVcf
  file vcfvepfile from outputVepVcf
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
	
	bash ${metrics_script} \
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
  file vepAnnot_log from vepAnnot_log.collect()
  file metrics_log from generate_metrics_file_log.collect()
  file iCallSV_log from iCallSV_log.collect()


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

        for f in ${vepAnnot_log};
        do
            echo "" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##-------------" \$f "-------------##" >>${params.merge_log_benchmark_files.final_logs}
            echo "##---------------------------------------------------------##" >>${params.merge_log_benchmark_files.final_logs}
            cat \$f >>${params.merge_log_benchmark_files.final_logs}
        done

        for f in ${iCallSV_log};
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
	//createFile(ON_START_FILE,"Workflow started at\t")
	file = new File(ON_START_FILE)
	file.append( "STARTED\n")

    file.append("Start Timestamp:\t$workflow.start\n")
    file.append("Project:\t$workflow.projectDir\n")
    file.append("Git info:\t$workflow.repository - $workflow.revision [$workflow.commitId]\n")
    file.append("Cmd line:\t$workflow.commandLine\n")

	}

def onsuccess(){
	createFile(ON_SUCCESS_FILE,"Workflow finished successfully at\t")
	println "SUCCESS"
}

def onfailure(){
	createFile(ON_FAILURE_FILE,"Workflow failed at\t")
	println "FAILED"
}

def createFile(name,message){
	file = new File(name)
	file.append(message)
	file.append("yyyy:mm:dd:hh:mm:ss\t" )
	file.append(new Date().format( 'yyyy:MM:dd:hh:mm:ss' ))
	file.append("\n")
}
