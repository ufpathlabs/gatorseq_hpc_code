#!/usr/bin/env nextflow

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
	
	if [ -d {params.temp_dir} ]; then rmdir {params.temp_dir}; fi
        mkdir -p {params.temp_dir};
		
	sambamba sort --nthreads=1 --memory-limit=8G  --tmpdir= ${params.temp_dir}   --out=${s}.bwa_map_sort.bam   --compression-level=0    ${s}.bwa_map.bam
	
	"""
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("_R"))
}

process merging {
cpus 4
echo true

input:
  file bamfiles from contigsBam.collect()

output:
  file 'merged.bam' into mapped_reads

"""
 samtools merge merged.bam $bamfiles
"""
}

