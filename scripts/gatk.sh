#!/bin/bash

echo "Running the script..."

echo 'Copying the reference genome...'
cp /home/tedwon/Sources/QA/data/silico/silico.fa .
rm *.fastq
rm *.dict
rm *.vcf
rm A.fq
rm B.fq
rm *.bam
rm *.intervals

echo '--------------- Generating Simulations ---------------'
wgsim -N 50000 -e 0 -r 0 r1.fa s1.fq s2.fq
#wgsim -r 0.10 silico_m.fa A.fq B.fq
echo 'Simulated reads generated'

echo '--------------- Building SAMTools Index ---------------'
samtools faidx silico.fa

echo '--------------- Generate Sequence Dictionary ---------------'
java -jar ../../Tools/Picard/picard.jar CreateSequenceDictionary REFERENCE=silico.fa OUTPUT=silico.dict

echo '--------------- Generate Index ---------------'
bwa index -a bwtsw silico.fa

echo '--------------- Alignment ---------------'
#bwa mem -t 10 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' silico.fa r1_1.fq r1_2.fq > aligned.sam
bwa mem -t 10 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' silico.fa s1.fq s2.fq > aligned.sam

echo '--------------- Sorting SAM to BAM ---------------'
java -jar ../../Tools/Picard/picard.jar SortSam INPUT=aligned.sam OUTPUT=sorted.bam SORT_ORDER=coordinate

echo '--------------- Marking Duplicates ---------------'
java -jar ../../Tools/Picard/picard.jar MarkDuplicates INPUT=sorted.bam OUTPUT=marked.bam METRICS_FILE=metrics.txt

echo '--------------- Build BAM Index ---------------'
java -jar ../../Tools/Picard/picard.jar BuildBamIndex INPUT=marked.bam

echo '--------------- Running RealignerTargetCreator ---------------'
java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R silico.fa -I marked.bam -o realigner.intervals

echo '--------------- Running IndelRealigner ---------------'
java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R silico.fa -I marked.bam -targetIntervals realigner.intervals -o realigned.bam

echo '--------------- Haplotype Analysis ---------------'
#java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R silico.fa -I realigned.bam --genotyping_mode DISCOVERY --heterozygosity 0.01 --defaultBaseQualities 30 -o haplotype.vcf
java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R silico.fa -I realigned.bam --genotyping_mode DISCOVERY -o haplotype.vcf







