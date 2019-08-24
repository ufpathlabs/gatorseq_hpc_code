#!/bin/bash

set -e

BAM_FILE=$1
VCF_FILE=$2
METRICS_OUT_FILE=$3
TOTAL_BAMUTILS_OUT=$4
TARGET_BAMUTILS_OUT=$5
PADDEDTARGET_BAMUTILS_OUT=$6
MERGED_TARGET_COVERAGE_FILE=$7
BAMTOOLS_STATS_OUT=$8
RTGTOOLS_STATS_OUT=$9
PADDING_SIZE=${10}
MERGED_TARGET_BED_FILE=${11}
PADDED_TARGET_BED_FILE=${12}
PROBE_MEAN_COVERAGE_OUT=${13}
SAMPLE_ID=${14}
SAMPLE_NAME=${15}
RUN_FOLDER=${16}
HUMAN_REF=${17}
PROBE_COVERAGE_OUT_FILE=${18}
TARGET_MANIFEST="designed-probe-coords IDT"


######################################################
#   Enrichment Summary
######################################################


Total_Target_Bases=$(eval cat $MERGED_TARGET_COVERAGE_FILE | wc -l )
Padded_Target_Bases=$(eval cat $PADDED_TARGET_BED_FILE | awk '{print $3-$2}' | awk '{s+=$1} END {print s}' )

######################################################
#   Read Level Enrichment
######################################################

Total_Reads=$(eval cat $TOTAL_BAMUTILS_OUT  |grep "^TotalReads" |cut -f2)
Total_Aligned_Reads=$(eval cat $TOTAL_BAMUTILS_OUT  |grep "^MappedReads" |cut -f2)
Percent_Total_Aligned_Reads=$(cat $TOTAL_BAMUTILS_OUT  |grep "^MappingRate" |cut -f2)
Target_Aligned_Reads=$(cat $TARGET_BAMUTILS_OUT  |grep "^MappedReads" |cut -f2)
Padded_Target_Aligned_Reads=$(cat $PADDEDTARGET_BAMUTILS_OUT |grep "^MappedReads" |cut -f2)
Read_Enrichment=$( echo -e "scale =2; 100 * $Target_Aligned_Reads / $Total_Aligned_Reads" | bc)
Padded_Read_Enrichment=$( echo -e "scale =2; 100 * $Padded_Target_Aligned_Reads / $Total_Aligned_Reads" | bc)
Median_Insert_Size=$(cat $BAMTOOLS_STATS_OUT |grep "Median insert size" |rev|cut -d ' ' -f1 |rev)
Percent_Proper_Pair=$(cat $TARGET_BAMUTILS_OUT  |grep "^ProperPair(%)" |cut -f2)

######################################################
#   Base Level Enrichment
######################################################

Total_Bases=$(eval cat $TOTAL_BAMUTILS_OUT  |grep "^TotalBases" |cut -f2)
Total_Aligned_Bases=$(eval cat $TOTAL_BAMUTILS_OUT  |grep "^BasesInMappedReads" |cut -f2)
Percent_Aligned_Bases=$( echo -e "scale =2; 100 * $Total_Aligned_Bases / $Total_Bases" | bc)
Target_Aligned_Bases=$(eval cat $TARGET_BAMUTILS_OUT  |grep "^BasesInMappedReads" |cut -f2)
Padded_Target_Aligned_Bases=$(eval cat $PADDEDTARGET_BAMUTILS_OUT  |grep "^BasesInMappedReads" |cut -f2)
Bases_Enrichment=$( echo -e "scale =2; 100 * $Target_Aligned_Bases / $Total_Aligned_Bases" | bc)
Padded_Bases_Enrichment=$( echo -e "scale =2; 100 * $Padded_Target_Aligned_Bases / $Total_Aligned_Bases" | bc)

######################################################
#   Coverage Summary for Target Regions
######################################################

#eval $DATAMASH_LOAD_CMD
#cat bam1.cov |cut -f6 |datamash min 1 q1 1 median 1 mean 1 q3 1 max 1 sstdev 1
Mean_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash mean 1 |xargs printf "%.*f\n" 2 )
Median_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash  median 1 )
Max_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash max 1)
Min_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash min 1)
q1_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash q1 1)
q3_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash q3 1)
sstdev_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |datamash sstdev 1 |xargs printf "%.*f\n" 2 )

Times_Mean=$( echo -e "scale =2; 0.2 * $Mean_Coverage" | bc)
#Uniformity_of_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk '$1 > $Times_Mean'|wc -l )
Uniformity_of_Coverage=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk  -v Times_Mean="$Times_Mean" '$1 > Times_Mean'|wc -l )
Percent_Uniformity_of_Coverage=$( echo -e "scale =2; 100 * $Uniformity_of_Coverage / $Total_Target_Bases" | bc)

Target_Coverage_at_50X=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk '$1 > 49'|wc -l )
Target_Coverage_at_100X=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk '$1 > 99'|wc -l )
Target_Coverage_at_300X=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk '$1 > 299'|wc -l )
Target_Coverage_at_500X=$(eval cat $MERGED_TARGET_COVERAGE_FILE |cut -f6 |awk '$1 > 499'|wc -l )
Percent_Target_Coverage_at_50X=$( echo -e "scale =2; 100 * $Target_Coverage_at_50X / $Total_Target_Bases" | bc)
Percent_Target_Coverage_at_100X=$( echo -e "scale =2; 100 * $Target_Coverage_at_100X / $Total_Target_Bases" | bc)
Percent_Target_Coverage_at_300X=$( echo -e "scale =2; 100 * $Target_Coverage_at_300X / $Total_Target_Bases" | bc)
Percent_Target_Coverage_at_500X=$( echo -e "scale =2; 100 * $Target_Coverage_at_500X / $Total_Target_Bases" | bc)

######################################################
#   Duplicate Summary
######################################################
Percent_Duplicate_Reads=$(cat $TOTAL_BAMUTILS_OUT  |grep "^DupRate" |cut -f2)

######################################################
#   Small Variants Summary
######################################################
Total_Passing=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Passed Filters" |cut -d ':' -f2|sed 's/^ //g')
Total_Failing=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Failed Filters" |cut -d ':' -f2|sed 's/^ //g')
SNVs=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^SNPs" |cut -d ':' -f2|sed 's/^ //g')
MNVs=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^MNPs" |cut -d ':' -f2|sed 's/^ //g')
Insertions=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Insertions" |cut -d ':' -f2|sed 's/^ //g')
Deletions=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Deletions" |cut -d ':' -f2|sed 's/^ //g')
Indels=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Indels" |cut -d ':' -f2|sed 's/^ //g')
SNV_Transitions_Transversions=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^SNP Transitions" |cut -d ':' -f2|sed 's/^ //g')
Total_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Total Het" |cut -d ':' -f2|sed 's/^ //g')
SNV_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^SNP Het" |cut -d ':' -f2|sed 's/^ //g')
MNV_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^MNP Het" |cut -d ':' -f2|sed 's/^ //g')
Insertion_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Insertion Het" |cut -d ':' -f2|sed 's/^ //g')
Deletion_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Deletion Het" |cut -d ':' -f2|sed 's/^ //g')
Indel_Het_Hom_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Indel Het" |cut -d ':' -f2|sed 's/^ //g')
Insertion_Deletion_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Insertion\/Deletion" |cut -d ':' -f2|sed 's/^ //g')
Indel_SNV_MNV_ratio=$(eval cat $RTGTOOLS_STATS_OUT  | grep "^Indel\/SNP" |cut -d ':' -f2|sed 's/^ //g')


echo -e "Enrichment Summary Report," >$METRICS_OUT_FILE
echo -e '"Note: All enrichment values are calculated without padding (sequence immediately upstream and downstream) unless otherwise stated. If any targeted region overlaps another region, the region positions will be adjusted to remove",' >>$METRICS_OUT_FILE
echo -e '"Note: PCR duplicate reads are removed from statistics.",' >>$METRICS_OUT_FILE
echo -e "," >>$METRICS_OUT_FILE
echo -e "Sample ID,$SAMPLE_ID" >>$METRICS_OUT_FILE
echo -e "Sample Name,$SAMPLE_NAME" >>$METRICS_OUT_FILE
echo -e "Runfolder,$RUN_FOLDER" >>$METRICS_OUT_FILE
echo -e "Reference Genome,$HUMAN_REF" >>$METRICS_OUT_FILE
echo -e "Target manifest,$TARGET_MANIFEST" >>$METRICS_OUT_FILE

echo -e "Total length of targeted reference,$Total_Target_Bases" >>$METRICS_OUT_FILE
echo -e "Padding size,$PADDING_SIZE" >>$METRICS_OUT_FILE

echo -e "Total PF reads(million),$Total_Reads" >>$METRICS_OUT_FILE
echo -e "Total aligned reads(million),$Total_Aligned_Reads" >>$METRICS_OUT_FILE
echo -e "Percent aligned reads,$Percent_Total_Aligned_Reads" >>$METRICS_OUT_FILE

echo -e "Percent duplicate paired reads,$Percent_Duplicate_Reads" >>$METRICS_OUT_FILE

echo -e "Targeted aligned reads(million),$Target_Aligned_Reads" >>$METRICS_OUT_FILE
echo -e "Padded target aligned reads(million),$Padded_Target_Aligned_Reads" >>$METRICS_OUT_FILE
echo -e "Read enrichment,$Read_Enrichment" >>$METRICS_OUT_FILE
echo -e "Padded read enrichment,$Padded_Read_Enrichment" >>$METRICS_OUT_FILE

echo -e "Total PF bases(million),$Total_Bases" >>$METRICS_OUT_FILE
echo -e "Percent Q30,NA" >>$METRICS_OUT_FILE
echo -e "Percent Q30 Aligned,NA" >>$METRICS_OUT_FILE
echo -e "Total aligned bases(million),$Total_Aligned_Bases" >>$METRICS_OUT_FILE
echo -e "Percent aligned bases,$Percent_Aligned_Bases" >>$METRICS_OUT_FILE
echo -e "Targeted aligned bases(million),$Target_Aligned_Bases" >>$METRICS_OUT_FILE
echo -e "Padded target aligned bases(million),$Padded_Target_Aligned_Bases" >>$METRICS_OUT_FILE
echo -e "Base enrichment,$Bases_Enrichment" >>$METRICS_OUT_FILE
echo -e "Padded base enrichment,$Padded_Bases_Enrichment" >>$METRICS_OUT_FILE
echo -e "Mean region coverage depth,$Mean_Coverage" >>$METRICS_OUT_FILE
echo -e "Uniformity of coverage (Pct > 0.2*mean),$Percent_Uniformity_of_Coverage" >>$METRICS_OUT_FILE

echo -e "Percent Target Coverage at 50X,$Percent_Target_Coverage_at_50X" >>$METRICS_OUT_FILE
echo -e "Percent Target Coverage at 100X,$Percent_Target_Coverage_at_100X" >>$METRICS_OUT_FILE
echo -e "Percent Target Coverage at 300X,$Percent_Target_Coverage_at_300X" >>$METRICS_OUT_FILE
echo -e "Percent Target Coverage at 500X,$Percent_Target_Coverage_at_500X" >>$METRICS_OUT_FILE

echo -e "Fragment length median,NA" >>$METRICS_OUT_FILE
echo -e "Fragment length min,NA" >>$METRICS_OUT_FILE
echo -e "Fragment length max,NA" >>$METRICS_OUT_FILE
echo -e "Fragment length SD,NA" >>$METRICS_OUT_FILE

echo -e "SNVs (All),$SNVs" >>$METRICS_OUT_FILE
echo -e "SNVs,$SNVs" >>$METRICS_OUT_FILE
echo -e "SNV Het/Hom ratio,$SNV_Het_Hom_ratio" >>$METRICS_OUT_FILE
echo -e "SNV Ts/Tv ratio,$SNV_Transitions_Transversions" >>$METRICS_OUT_FILE
echo -e "Indels (All),$Indels" >>$METRICS_OUT_FILE
echo -e "Indels,$Indels" >>$METRICS_OUT_FILE
echo -e "Indel Het/Hom ratio,$Indel_Het_Hom_ratio" >>$METRICS_OUT_FILE
echo -e "Insertions (All),$Insertions" >>$METRICS_OUT_FILE
echo -e "Insertions,$Insertions" >>$METRICS_OUT_FILE
echo -e "Insertion Het/Hom ratio,$Insertion_Het_Hom_ratio" >>$METRICS_OUT_FILE
echo -e "Deletions (All),$Deletions" >>$METRICS_OUT_FILE
echo -e "Deletions,$Deletions" >>$METRICS_OUT_FILE
echo -e "Deletion Het/Hom ratio,$Deletion_Het_Hom_ratio" >>$METRICS_OUT_FILE


echo -e "#Enrichment: NA%," >$PROBE_COVERAGE_OUT_FILE
echo -e "#Reads: NA,NA,NA" >>$PROBE_COVERAGE_OUT_FILE
echo -e "#Chromosome,Start,Stop,RegionID,MeanCoverage,StdDevCoverage" >>$PROBE_COVERAGE_OUT_FILE

cat $PROBE_MEAN_COVERAGE_OUT |sed 's/\t/,/g' >>$PROBE_COVERAGE_OUT_FILE 




#echo -e "Total Length of Padded Targeted Reference:\t$Padded_Target_Bases" >>$METRICS_OUT_FILE
#echo -e "###   Read Level Enrichment\t###" >>$METRICS_OUT_FILE

#echo -e "Median Insert Size:\t$Median_Insert_Size" >>$METRICS_OUT_FILE
#echo -e "Percent Proper Pair:\t$Percent_Proper_Pair" >>$METRICS_OUT_FILE
#echo -e " \t " >>$METRICS_OUT_FILE

#echo -e "###   Base Level Enrichment\t###" >>$METRICS_OUT_FILE
#echo -e " \t " >>$METRICS_OUT_FILE


#echo -e "###   Coverage Summary for Target Regions\t###" >>$METRICS_OUT_FILE
#echo -e "Mean Coverage:\t$Mean_Coverage" >>$METRICS_OUT_FILE
#echo -e "Median Coverage:\t$Median_Coverage" >>$METRICS_OUT_FILE
#echo -e "Max Coverage:\t$Max_Coverage" >>$METRICS_OUT_FILE
#echo -e "Min Coverage:\t$Min_Coverage" >>$METRICS_OUT_FILE
#echo -e "q1 Coverage:\t$q1_Coverage" >>$METRICS_OUT_FILE
#echo -e "q3 Coverage:\t$q3_Coverage" >>$METRICS_OUT_FILE
#echo -e "sstdev Coverage:\t$sstdev_Coverage" >>$METRICS_OUT_FILE

#echo -e "Target Coverage at 100X:\t$Target_Coverage_at_100X" >>$METRICS_OUT_FILE
#echo -e "Target Coverage at 300X:\t$Target_Coverage_at_300X" >>$METRICS_OUT_FILE
#echo -e "Target Coverage at 500X:\t$Target_Coverage_at_500X" >>$METRICS_OUT_FILE

#echo -e "Percent Target Coverage at 100X:\t$Percent_Target_Coverage_at_100X" >>$METRICS_OUT_FILE
#echo -e "Percent Target Coverage at 300X:\t$Percent_Target_Coverage_at_300X" >>$METRICS_OUT_FILE
#echo -e "Percent Target Coverage at 500X:\t$Percent_Target_Coverage_at_500X" >>$METRICS_OUT_FILE
#echo -e " \t " >>$METRICS_OUT_FILE

#echo -e "###   Duplicate Summary\t###" >>$METRICS_OUT_FILE

#echo -e " \t " >>$METRICS_OUT_FILE


#echo -e "###   Small Variants Summary\t###" >>$METRICS_OUT_FILE

#echo -e "Total Passed Variants:\t$Total_Passing" >>$METRICS_OUT_FILE
#echo -e "Total Failed Variants:\t$Total_Failing" >>$METRICS_OUT_FILE
#echo -e "MNVs:\t$MNVs" >>$METRICS_OUT_FILE
#echo -e "Indel/SNV+MNV ratio:\t$Indel_SNV_MNV_ratio" >>$METRICS_OUT_FILE
#echo -e "Total Het/Hom ratio:\t$Total_Het_Hom_ratio" >>$METRICS_OUT_FILE
#echo -e "MNV Het/Hom ratio:\t$MNV_Het_Hom_ratio" >>$METRICS_OUT_FILE
#echo -e "Insertion/Deletion ratio:\t$Insertion_Deletion_ratio" >>$METRICS_OUT_FILE

#echo -e " \t " >>$METRICS_OUT_FILE

#eval sed -i 's/\(/[/g' $METRICS_OUT_FILE
#eval sed -i 's/\)/]/g' $METRICS_OUT_FILE


