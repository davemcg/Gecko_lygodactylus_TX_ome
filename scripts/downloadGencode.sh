#!/bin/bash
refGenome_url=$1
refGFF_url=$2
refGenome=$3
refGFF=$4
gentmp={refGenome}.tmp
GFFtmp={refGFF}.tmp
wget -O - $refGenome_url | gunzip -c > $gentmp
wget -O - $refGFF_url | gunzip -c > $GFFtmp
python3 scripts/clean_genome_names.py $gentmp $refGenome
gtf=$GFFtmp
fasta=$refGenome
new_gtf=$refGFF
gtf_recs=/tmp/gtf_recs
fasta_recs=/tmp/fasta_recs
badrec=/tmp/tx_to_remove
grep '>' $fasta | tr -d '>' > $fasta_recs
awk '{print $1 }' $gtf | sort -u | grep -v '#' - >$gtf_recs
grep -v -Ff $fasta_recs $gtf_recs > $badrec
grep -v -Ff $badrec $gtf >  $new_gtf
rm *.tmp
