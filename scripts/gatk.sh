#!/bin/bash

echo "------------ Running the script ------------"
cp /home/tedwon/Sources/QA/data/silico/silico.fa .

rm *.fastq
rm *.dict
rm *.vcf
rm A.fq
rm B.fq
rm *.bam
rm *.intervals

echo '--------------- Generating Simulations ---------------'
wgsim silico.fa f1.fq f2.fq
echo 'Simulated reads generated'

echo '--------------- Building SAMTools Index ---------------'
samtools faidx silico.fa

echo '--------------- Generate Sequence Dictionary ---------------'
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.95/CreateSequenceDictionary.jar REFERENCE=silico.fa OUTPUT=silico.dict

echo '--------------- Generate Index ---------------'
bwa index silico.fa

echo '--------------- Alignment ---------------'
#bwa mem -t 10 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' silico.fa r1_1.fq r1_2.fq > aligned.sam
bwa mem -t 10 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' silico.fa f1.fq f2.fq > aligned.sam

echo '--------------- Sorting SAM to BAM ---------------'
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.95/SortSam.jar INPUT=aligned.sam OUTPUT=sorted.bam SORT_ORDER=coordinate

echo '--------------- Marking Duplicates ---------------'
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.95/MarkDuplicates.jar INPUT=sorted.bam OUTPUT=marked.bam METRICS_FILE=metrics.txt

echo '--------------- Build BAM Index ---------------'
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.95/BuildBamIndex.jar INPUT=sorted.bam

echo '--------------- Running RealignerTargetCreator ---------------'
java -jar /share/ClusterShare/software/contrib/gi/gatk/3.1/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R silico.fa -I marked.bam -o realigner.intervals

echo '--------------- Running IndelRealigner ---------------'
java -jar /share/ClusterShare/software/contrib/gi/gatk/3.1/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T IndelRealigner -R silico.fa -I marked.bam -targetIntervals realigner.intervals -o realigned.bam

echo '--------------- Haplotype Analysis ---------------'
#java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R silico.fa -I realigned.bam --genotyping_mode DISCOVERY --heterozygosity 0.01 --defaultBaseQualities 30 -o haplotype.vcf
java -jar /share/ClusterShare/software/contrib/gi/gatk/3.1/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T HaplotypeCaller -R silico.fa -I sorted.bam --genotyping_mode DISCOVERY -o haplotype.vcf







