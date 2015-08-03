rm -rf RNA_Simulation

python simulate.py RNA

qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh A1.pbs
qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh A2.pbs
qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh A3.pbs
qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh B1.pbs
qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh B3.pbs
qsub -M twon@garvan.org.au -m beas -V -cwd -pe smp 1 -j y -S /bin/sh B3.pbs

#rm assemblies.txt
#echo RNA_A_1_1/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_A_1_2/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_A_1_3/assembly/transcripts.gtf  >> assemblies.txt
#echo RNA_B_10_1/assembly/transcripts.gtf >> assemblies.txt
#echo RNA_B_10_2/assembly/transcripts.gtf >> assemblies.txt
#echo RNA_B_10_3/assembly/transcripts.gtf >> assemblies.txt

#cuffmerge -g ../data/rna/standards.gtf -s combined.fa -p 4 assemblies.txt
