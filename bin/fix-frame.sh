#!/usr/bin/env bash
reads=$1
#reads=n2
#reads=emales-fl-d3
prots=$2
#prots=a1
#prots=emales
set -x;

# seqkit subseq --bed <(echo -e "EMALE01_Cflag_c017B\t3800\t6530\n") ../emales.fna > n2.ffn         

diamond makedb --in $prots.faa --db $prots
diamond blastx -q $reads.fna -d $prots --range-culling -F 15 --more-sensitive -k 1 --query-gencode 4 \
 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen cigar sseq btop \
 --out ${reads}_${prots}.tsv 2>/dev/null
./fix-frame.pl -D $reads.fna ${reads}_${prots}.tsv > $reads-f1.fna

# mycoplasma has genetic code 4
rm -r mmseqs-tmp

prodigal -g 4 -t mbovis-prdgl.train -i ${reads}.fna -f gff -o ${reads}.gff0
gff-clean ${reads}.gff0 > ${reads}.gff
gff2cds --aa --fna $reads.fna --type CDS --source Prodigal_v2.6.3 $reads.gff > $reads.faa
mmseqs easy-search --format-mode 2 $reads.faa $prots.faa $reads-faa.tsv mmseqs-tmp

prodigal -g 4 -t mbovis-prdgl.train -i ${reads}-f1.fna -f gff -o ${reads}-f1.gff0
gff-clean ${reads}-f1.gff0 > ${reads}-f1.gff
gff2cds --aa --fna $reads-f1.fna --type CDS --source Prodigal_v2.6.3 $reads-f1.gff > $reads-f1.faa
mmseqs easy-search --format-mode 2 $reads-f1.faa $prots.faa $reads-f1-faa.tsv mmseqs-tmp



