#!/bin/bash

while read line; do python /home/shardwick/DEXSeq/inst/python_scripts/dexseq_count.py -p yes -s reverse -r name RNA.v1.DEXSeq.gff ../topout_$line/$line.accepted_hits.sortedbyName.sam $line.DEXSeq.counts.txt; done < ../samples_list.txt
