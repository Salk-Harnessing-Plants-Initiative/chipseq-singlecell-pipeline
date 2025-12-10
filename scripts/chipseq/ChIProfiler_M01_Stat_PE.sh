#!/bin/bash
# ChIP-seq analysis pipeline - Statistics module
# Original author: gxu@lawlab
# v4 changes compared to v2: corrected small errors
# Adapted for container environment by Salk HPI
set -euo pipefail

cpu_core=${1:-12} # how many cpu cores used in parallel

################################# count read number from the bam file
################################# output 1_res_read_mapping_stat.txt

# trimming info
: > trim_info_R1.txt
for file in `cat fastqList.txt`; do echo $file\_R1.fastq.gz_trimming_report.txt; cat $file\_R1.fastq.gz_trimming_report.txt ; done > trim_info_R1.txt
awk '/Total reads processed:/ {print $4}' trim_info_R1.txt > trim_info_R1_raw_reads.txt
awk '/Reads with adapters:/ {print $4,$5}' trim_info_R1.txt > trim_info_R1_reads_wAdapters.txt
awk '/Reads written/ {print $5,$6}' trim_info_R1.txt > trim_info_R1_reads_writen.txt
: > trim_info_R2.txt
for file in `cat fastqList.txt`; do echo $file\_R2.fastq.gz_trimming_report.txt; cat $file\_R2.fastq.gz_trimming_report.txt ; done > trim_info_R2.txt
awk '/Total reads processed:/ {print $4}' trim_info_R2.txt > trim_info_R2_raw_reads.txt
awk '/Reads with adapters:/ {print $4,$5}' trim_info_R2.txt > trim_info_R2_reads_wAdapters.txt
awk '/Reads written/ {print $5,$6}' trim_info_R2.txt > trim_info_R2_reads_writen.txt
for file in trim_info_R1_*.txt; do sed -i 's/,//g' $file; done
for file in trim_info_R2_*.txt; do sed -i 's/,//g' $file; done

#sorted bam file
rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c {}_sorted.bam > tmp_counts/{}.txt' ::: $(cat fastqList.txt)
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_total_read.txt

rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c -F 260 {}_sorted.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_mapped_read.txt

# sorted_rmdup bam file
rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c {}_sorted_rmdup.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_rmdup_total_read.txt

rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c -F 260 {}_sorted_rmdup.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_rmdup_mapped_read.txt

# sorted_rmdup_uniq bam file
rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c {}_sorted_rmdup_uniq.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_rmdup_uniq_total_read.txt

rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c -F 260 {}_sorted_rmdup_uniq.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_rmdup_uniq_mapped_read.txt

# sort_uniq bam file
rm -rf tmp_counts && mkdir tmp_counts
parallel -N1 -j $cpu_core 'samtools view -c {}_sorted_uniq.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_uniq_total_read.txt

rm -rf tmp_counts && mkdir tmp_counts
: > numbers_sorted_uniq_mapped_read.txt
parallel -N1 -j $cpu_core 'samtools view -c -F 260 {}_sorted_uniq.bam > tmp_counts/{}.txt' ::: `cat fastqList.txt`
while read sample; do
    cat tmp_counts/${sample}.txt
done < fastqList.txt > numbers_sorted_uniq_mapped_read.txt

paste fastqList.txt \
	trim_info_R1_raw_reads.txt \
	trim_info_R1_reads_wAdapters.txt \
	trim_info_R1_reads_writen.txt \
	trim_info_R2_raw_reads.txt \
	trim_info_R2_reads_wAdapters.txt \
	trim_info_R2_reads_writen.txt \
	numbers_sorted_total_read.txt \
	numbers_sorted_mapped_read.txt \
	numbers_sorted_rmdup_total_read.txt \
	numbers_sorted_rmdup_mapped_read.txt \
	numbers_sorted_rmdup_uniq_total_read.txt \
	numbers_sorted_rmdup_uniq_mapped_read.txt \
	numbers_sorted_uniq_total_read.txt \
	numbers_sorted_uniq_mapped_read.txt > 1_res_read_mapping_stat.txt
	
sed -i '1isample_name\tR1_read\tR1_w_adapter\tR1_writen\tR2_read\tR2_adapter\tR2_writen\tsorted_total_read\tsorted_mapped_read\tsorted_rmdup_total_read\tsorted_rmdup_mapped_read\tsorted_rmdup_uniq_total_read\tsorted_rmdup_uniq_mapped_read\tsorted_uniq_total_read\tsorted_uniq_mapped_read' 1_res_read_mapping_stat.txt

rm -rf tmp_counts 
rm ./trim_info_*.txt


