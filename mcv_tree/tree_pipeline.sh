#!/bin/bash

#download and unzip mm10 SNP and indel VCF files
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
gunzip mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
gunzip mgp.v5.merged.snps_all.dbSNP142.vcf.gz

#subset to chromosome 7
vcftools --vcf mgp.v5.merged.snps_all.dbSNP142.vcf --chr 7 --recode --out chr7_snps --keep founders.txt
vcftools --vcf mgp.v5.merged.indels.dbSNP142.normed.vcf --chr 7 --recode --out chr7_indels --keep founders.txt

#identify genomic ranges without recombination using 4-gamete test on high-quality, biallelic SNPs
Rscript four_gamete_test.R 7

#construct pseudogenome for genomic range without recombination around causal gene Hbb-b1
#RefSeq location for Hbb-bs: chr7:103826532-103827929 (https://genome.ucsc.edu/cgi-bin/hgc?hgsid=751258183_AbtpzDad45UaAie5k4AaxbEN5xlH&c=chr7&l=103826531&r=103827929&o=103826531&t=103827929&g=refGene&i=NM_001278161)
#this corresponds to line 3776 of chr7_ranges.txt (3777 including header)
sed -n '3777p' chr7_ranges.txt

#subset SNP and indel VCF files using range from pervious step
vcftools --vcf chr7_snps.recode.vcf --chr 7 --from-bp 103807679 --to-bp 103831178 --recode --out chr7_3776_snps --keep founders.txt
vcftools --vcf chr7_indels.recode.vcf --chr 7 --from-bp 103807679 --to-bp 103831178 --recode --out chr7_3776_indels --keep founders.txt

#download and unzip mm10 Chr7 fasta
wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.7.fa.gz
gunzip Mus_musculus.GRCm38.dna.chromosome.7.fa.gz

#create pseudogenomes
Rscript construct_pseudogenome.R chr7_3776 7 103807679 103831178

#create BEAST input file
sh /nas02/home/w/c/wcrouse/BEASTGenv1.0.2/bin/beastgen -D "chain_length=1100000,log_every=1000" beast1_custom_template.xml chr7_3776.fa chr7_3776.xml

#run BEAST
sh /nas02/home/w/c/wcrouse/BEASTv1.8.3/bin/beast -beagle_off -overwrite chr7_3776.xml
rm $JOBID".ops"

#convert trees to coalescent units
Rscript convert_to_coalescent.R