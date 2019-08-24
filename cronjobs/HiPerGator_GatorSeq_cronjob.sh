#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. /etc/profile.d/modules.sh
. /etc/profile.d/slurm.sh

export PATH="/usr/sbin/mksquashfs:$PATH"

#export SLURM_TMPDIR=/ufrc/chamala/share/GatorSeq_Share/GatorSeq_Analysis
#export $SCRIPT_DIR
echo $(date)
HPC_ANALYSIS_FOLDER=$SCRIPT_DIR"/../analysis"
#SBATCH_CMD='sbatch --export=SLURM_TMPDIR="/ufrc/chamala/share/GatorSeq_Share/GatorSeq_Analysis" '

gatorseq_job_scripts=$HPC_ANALYSIS_FOLDER'/*.cronjob.sh'
#echo $gatorseq_job_scripts

for gatorseq_job in $( ls $gatorseq_job_scripts ); do
    log_file=$gatorseq_job".log"
    job_name=$( echo $gatorseq_job | rev |cut -d '/' -f1 |rev )
    slurm_output_prefix=$gatorseq_job".slurm.%A_%a.log"
    #SBATCH_CMD="sbatch --job-name=$job_name --output=$slurm_output_prefix --ntasks=1 --mem=6gb --time=24:00:00 --qos=chamala-b --account=chamala --export=SCRIPT_DIR=$SCRIPT_DIR"
    SBATCH_CMD="sbatch --job-name=$job_name --output=$slurm_output_prefix --ntasks=1 --mem=6gb --time=24:00:00 --qos=chamala-b --account=chamala "
    #echo $log_file
    #echo $job_name
    #echo $slurm_output_prefix
    if [ ! -f "$log_file" ]
    then
        #bashCommand="nohup sh "$gatorseq_job" &>"$log_file" &"
        #bashCommand="sh "$gatorseq_job" &>"$log_file
        #bashCommand="echo; ulimit -c -l; echo; ulimit -c unlimited; echo;ulimit -c -l;sh "$gatorseq_job" &>"$log_file
        bashCommand=$SBATCH_CMD" "$gatorseq_job" &>"$log_file
        echo $bashCommand
        eval $bashCommand
    fi
done

echo "" 
