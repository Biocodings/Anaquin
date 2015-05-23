python simulate.py RNA

tophat -o RNA_A_1_1/aligned  combined RNA_A_1_1/reads/simulated_1.fastq  RNA_A_1_1/reads/simulated_2.fastq
tophat -o RNA_A_1_2/aligned  combined RNA_A_1_2/reads/simulated_1.fastq  RNA_A_1_2/reads/simulated_2.fastq
tophat -o RNA_A_1_3/aligned  combined RNA_A_1_3/reads/simulated_1.fastq  RNA_A_1_3/reads/simulated_2.fastq
tophat -o RNA_B_20_1/aligned combined RNA_B_20_1/reads/simulated_1.fastq RNA_B_20_1/reads/simulated_2.fastq
tophat -o RNA_B_20_2/aligned combined RNA_B_20_2/reads/simulated_1.fastq RNA_B_20_2/reads/simulated_2.fastq
tophat -o RNA_B_20_3/aligned combined RNA_B_20_3/reads/simulated_1.fastq RNA_B_20_3/reads/simulated_2.fastq

cufflinks -p 4 -o RNA_A_1_1/assembly  RNA_A_1_1/aligned/accepted_hits.bam
cufflinks -p 4 -o RNA_A_1_2/assembly  RNA_A_1_2/aligned/accepted_hits.bam
cufflinks -p 4 -o RNA_A_1_3/assembly  RNA_A_1_3/aligned/accepted_hits.bam
cufflinks -p 4 -o RNA_B_20_1/assembly RNA_B_20_1/aligned/accepted_hits.bam
cufflinks -p 4 -o RNA_B_20_2/assembly RNA_B_20_2/aligned/accepted_hits.bam
cufflinks -p 4 -o RNA_B_20_3/assembly RNA_B_20_3/aligned/accepted_hits.bam

rm assemblies.txt
echo RNA_A_1_1/assembly/transcripts.gtf  >> assemblies.txt
echo RNA_A_1_2/assembly/transcripts.gtf  >> assemblies.txt
echo RNA_A_1_3/assembly/transcripts.gtf  >> assemblies.txt
echo RNA_B_20_1/assembly/transcripts.gtf >> assemblies.txt
echo RNA_B_20_2/assembly/transcripts.gtf >> assemblies.txt
echo RNA_B_20_3/assembly/transcripts.gtf >> assemblies.txt

cuffmerge -g ../data/rna/standards.gtf -s combined.fa -p 4 assemblies.txt
