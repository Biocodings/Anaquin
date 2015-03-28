#!/usr/bin/python

import os
import sys
import math

os.system('rm -rf /home/tedwon/Projects/Temp')
os.system('mkdir /home/tedwon/Projects/Temp')

os.system('cp /home/tedwong/Sources/QA/data/silco/silico.fa /home/tedwon/Projects/Temp')

os.system('cd /home/tedwon/Projects/Temp')

os.system('wgsim -e0 -r0 -R0 silico.fa A.fastq B.fastq')

os.system('samtools faidx silico.fa')

os.system('java -jar ../../Tools/Picard/picard.jar CreateSequenceDictionary REFERENCE=silico.fa OUTPUT=silico.dict')

os.system('bwa index -a bwtsw silico.fa')

os.system("bwa mem -t 20 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' silico.fa A.fastq B.fastq > aligned.sam")

os.system('java -jar ../../Tools/Picard/picard.jar SortSam INPUT=aligned.sam OUTPUT=sorted.bam SORT_ORDER=coordinate')

os.system('java -jar ../../Tools/Picard/picard.jar MarkDuplicates INPUT=sorted.bam OUTPUT=marked.bam METRICS_FILE=metrics.txt')

os.system('java -jar ../../Tools/Picard/picard.jar BuildBamIndex INPUT=marked.bam')

os.system('java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R silico.fa -I marked.bam -o realigner.intervals')

os.system('java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R silico.fa -I marked.bam -targetIntervals realigner.intervals -o realigned.bam')

os.system('java -jar ../../Tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R silico.fa -I realigned.bam --genotyping_mode DISCOVERY --heterozygosity 0.01 --defaultBaseQualities 30 -o haplotype.vcf')








