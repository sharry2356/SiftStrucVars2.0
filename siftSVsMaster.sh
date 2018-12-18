#!/bin/bash


windowSize=5000 # put this analysis in first

Rscript siftSVs1.R

###LINUX
grep -F -f GOIs_and_Correlatedgenes_list.txt annotatedgenome.gff > GOIs_and_Cs_annotatedgenome.bed
# subset mRNA and CDS lines from GOIs_and_Cs_annotatedgenome.bed
grep ID=mRNA GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_mRNA.bed 
grep ID=CDS GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_CDS.bed
# **** Find the CDS that intersect SVs and then concatonate with full set of CDS of those genes/mRNAs
bedtools window -a GOIs_and_Cs_CDS.bed -b bedformat_SVs.bed -w 0 -u  > GOIsandCs_CDS_that_intersect_SVs.bed
grep -E -e "Solyc[0-9]+g[0-9]+" GOIsandCs_CDS_that_intersect_SVs.bed -o | uniq > genes_with_CDS_overlapped_by_SVs.txt
grep -F -f genes_with_CDS_overlapped_by_SVs.txt GOIs_and_Cs_CDS.bed > Full_CDS_sets.bed
# which SVs intersect which CDS (and size of overlap)
bedtools intersect -a  bedformat_SVs.bed -b GOIs_and_Cs_CDS.bed -wo  > overlapping_SVs_and_GOIsandCs_CDS.bed
# Process for finding SVs that are 2000 bp upstream promoter or intersect 5'UTR (and don't intersect CDS of any GOIs/Cs 
grep ID=five_prime_UTR GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_five_prime_UTR.bed
bedtools intersect -a bedformat_SVs.bed -b GOIs_and_Cs_CDS.bed -v > SVs_no_intersect_CDS.bed
bedtools intersect -a SVs_no_intersect_CDS.bed -b GOIs_and_Cs_five_prime_UTR.bed -v > SVs_no_interesct_CDS_or_5_prime.bed
bedtools window -a SVs_no_interesct_CDS_or_5_prime.bed -b GOIs_and_Cs_five_prime_UTR.bed -sw -r 0 -l $windowSize > SVs_upstream_of_genes.bed
bedtools intersect -a SVs_no_intersect_CDS.bed -b GOIs_and_Cs_five_prime_UTR.bed -wa -wb > SVs_only_intersect_fives.bed


Rscript siftSVs2.R 

head -n 1 Anthesis.txt > Anthesis_CDS.txt 
head -n 1 Meristems.txt > Meristems_CDS.txt 
head -n 1 Anthesis.txt > Anthesis_RE.txt 
head -n 1 Meristems.txt > Meristems_RE.txt 

grep -F -f SVs_intersect_CDS_genelist.txt Anthesis.txt >> Anthesis_CDS.txt
grep -F -f SVs_intersect_CDS_genelist.txt Meristems.txt >> Meristems_CDS.txt
grep -F -f SVs_intersect_RE_genelist.txt Anthesis.txt >> Anthesis_RE.txt
grep -F -f SVs_intersect_RE_genelist.txt Meristems.txt >> Meristems_RE.txt 

Rscript siftSVs3.R

mkdir -p Intermediates Outputs
mv Anthesis_CDS.txt Anthesis_RE.txt Meristems_CDS.txt Meristems_RE.txt SVs_intersect_RE_genelist.txt SVs_intersect_CDS_genelist.txt GOIs_and_Cs_annotatedgenome.bed GOIs_and_Cs_mRNA.bed GOIs_and_Cs_CDS.bed GOIsandCs_CDS_that_intersect_SVs.bed genes_with_CDS_overlapped_by_SVs.txt Full_CDS_sets.bed overlapping_SVs_and_GOIsandCs_CDS.bed GOIs_and_Cs_five_prime_UTR.bed SVs_no_intersect_CDS.bed SVs_no_interesct_CDS_or_5_prime.bed SVs_upstream_of_genes.bed SVs_only_intersect_fives.bed bedformat_SVs.bed GOIs_and_Correlatedgenes_list.txt annotatedgenomeinter.txt annotatedgenome.gff -t Intermediates
mv Correlation_Table SVs_intersect_CDS.txt SVs_intersect_RE.txt Relevent_Expression_SVs_intersect_CDS.txt Relevent_Expression_SVs_intersect_RE.txt -t Outputs

wc -l *.* > wc_l_Inputs_Intermediates_Outputs 
wc -l Intermediates/*.* >> wc_l_Inputs_Intermediates_Outputs
wc -l Outputs/*.* >> wc_l_Inputs_Intermediates_Outputs 

mv wc_l_Inputs_Intermediates_Outputs -t Outputs
