#!/bin.bash

echo "Running the script..."

echo "Generating the output folder..."
rm -rf /home/tedwon/Projects/Temp
mkdir /home/tedwon/Projects/Temp
cd /home/tedwon/Projects/Temp

echo 'Copying the reference genome...'
cp /home/tedwon/Sources/QA/data/silico/silico.fa /home/tedwon/Projects/Temp


wgsim -e0 -r0 -R0 silico.fa A.fastq B.fastq



