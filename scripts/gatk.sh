echo "Running the script..."

rm -rf /home/tedwon/Projects/Temp


mkdir /home/tedwon/Projects/Temp
cd /home/tedwon/Projects/Temp

cp /home/tedwon/Sources/QA/data/silico/silico.fa /home/tedwon/Projects/Temp


wgsim -e0 -r0 -R0 silico.fa A.fastq B.fastq



