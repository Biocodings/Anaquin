rm -rf RNA_Simulation
python simulate.py RNA

python ~/submit.py A1.pbs
python ~/submit.py A2.pbs
python ~/submit.py A3.pbs
python ~/submit.py B1.pbs
python ~/submit.py B2.pbs
python ~/submit.py B3.pbs

#rm assemblies.txt
#echo RNA_A_1_1/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_A_1_2/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_A_1_3/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_B_10_1/assembly/transcripts.gtf >> assemblies.txt
#echo RNA_B_10_2/assembly/transcripts.gtf >> assemblies.txt
#echo RNA_B_10_3/assembly/transcripts.gtf >> assemblies.txt

#cuffmerge -g ../data/rna/standards.gtf -s combined.fa -p 4 assemblies.txt
