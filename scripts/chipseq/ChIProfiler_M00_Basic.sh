#!/bin/bash
# ChIP-seq analysis pipeline - Basic workflow
# Original author: gxu@lawlab
# Adapted for container environment by Salk HPI
set -euo pipefail

#################################
################################# set parameters
# Reference genome - mount your bowtie2 index to /references
reference="${REFERENCE:-/references/bowtie2_index}"
cpu_core="${CPU_CORES:-12}" # how many cpu cores used in parallel
SEorPE="${READ_TYPE:-PE}" # PE or SE
callPeaksSample="${IP_SAMPLE:-IP_col0}"
callPeaksCK="${INPUT_SAMPLE:-input_col0}"
callPeakSuffix="${TAG_SUFFIX:-_sorted_rmdup_uniq_TagDir}" # can be "_sorted_rmdup_uniq_TagDir" or "_sorted_uniq_TagDir"

# Script directory for calling other modules
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

################################# trim
################################# output: SE: *trimmed.fq.gz, *trimming_report.txt; PE: *_R1_val_1.fq.gz, _R2_val_2.fq.gz, *trimming_report.txt
if [ $SEorPE = "PE" ]
then
	ls *R1.fastq.gz | sed 's/_R1.fastq.gz//g' > fastqList.txt
	parallel -j $cpu_core 'trim_galore --paired {}_R1.fastq.gz {}_R2.fastq.gz' ::: `cat fastqList.txt`
else
	ls *.fastq.gz | sed 's/.fastq.gz//g' > fastqList.txt
	parallel -j $cpu_core 'trim_galore {}.fastq.gz' ::: `cat fastqList.txt`
fi

################################# bowtie2 mapping
################################# output: *_sorted.bam
if [ $SEorPE = "PE" ]
then
	parallel -j $cpu_core 'bowtie2 -p 1 -q --local -x '$reference' -1 {}_R1_val_1.fq.gz -2 {}_R2_val_2.fq.gz | samtools view -bS | samtools sort -o {}_sorted.bam' ::: `cat fastqList.txt`
else
	parallel -j $cpu_core 'bowtie2 -p 1 -q --local -x '$reference' -U {}_trimmed.fq.gz | samtools view -bS | samtools sort -o {}_sorted.bam' ::: `cat fastqList.txt`
fi


################################# bam --> uniq bam
#################################
#parallel -j $cpu_core 'samtools rmdup -s {}_sorted.bam {}_sorted_rmdup.bam' ::: `cat fastqList.txt` # This is used for SE
parallel -j $cpu_core '
  samtools sort -n {}_sorted.bam -o {}_namesorted.bam &&
  samtools fixmate -m {}_namesorted.bam {}_fixmate.bam &&
  samtools sort {}_fixmate.bam -o {}_positionsorted.bam &&
  samtools markdup {}_positionsorted.bam {}_sorted_markdup.bam &&
  samtools view -b -F 1024 {}_sorted_markdup.bam > {}_sorted_rmdup.bam
' ::: `cat fastqList.txt`

parallel -j $cpu_core 'samtools view -b -q 20 -o {}_sorted_rmdup_uniq.bam {}_sorted_rmdup.bam' ::: `cat fastqList.txt`
parallel -j $cpu_core 'samtools view -b -q 20 -o {}_sorted_uniq.bam {}_sorted.bam' ::: `cat fastqList.txt`

#parallel -j $cpu_core 'samtools index {}' ::: *sorted_rmdup.bam
parallel -j $cpu_core 'samtools index {}' ::: *sorted_rmdup_uniq.bam
parallel -j $cpu_core 'samtools index {}' ::: *sorted_uniq.bam


################################# make tag directory
#################################
parallel -j $cpu_core 'makeTagDirectory {}_sorted_rmdup_uniq_TagDir -format sam -mis 2 -unique {}_sorted_rmdup_uniq.bam' ::: `cat fastqList.txt`
parallel -j $cpu_core 'makeTagDirectory {}_sorted_uniq_TagDir -format sam -mis 2 -unique {}_sorted_uniq.bam' ::: `cat fastqList.txt`

################################# make UCSC tracks for genome browser/IGV
#################################
parallel -j $cpu_core 'makeUCSCfile {}_sorted_rmdup_uniq_TagDir none -fragLength given -norm 10000000 -o auto' ::: `cat fastqList.txt`
parallel -j $cpu_core 'makeUCSCfile {}_sorted_uniq_TagDir none -fragLength given -norm 10000000 -o auto' ::: `cat fastqList.txt`
cp ./*/*.bedGraph.gz ./
rename _TagDir.ucsc.bedGraph.gz .bedGraph.gz *_TagDir.ucsc.bedGraph.gz
parallel -j $cpu_core 'igvtools toTDF {}_sorted_rmdup_uniq.bedGraph.gz {}_sorted_rmdup_uniq.tdf tair10' ::: `cat fastqList.txt`
parallel -j $cpu_core 'igvtools toTDF {}_sorted_uniq.bedGraph.gz {}_sorted_uniq.tdf tair10' ::: `cat fastqList.txt`

################################# M01 count read numbers
################################# output 1_res_read_mapping_stat.txt
"${SCRIPT_DIR}/ChIProfiler_M01_Stat_PE.sh" "$cpu_core"

################################# Module03 call peaks
#################################
"${SCRIPT_DIR}/ChIProfiler_M03_CallPeaks.sh" "$callPeaksSample" "$callPeaksCK" "$callPeakSuffix" 3

################################# delete temp files manually
#################################
# if the *_sorted_rmdup.bam files looks good, delete the following
rm ./*fixmate.bam
rm ./*namesorted.bam
rm ./*positionsorted.bam
rm ./*sorted_markdup.bam










