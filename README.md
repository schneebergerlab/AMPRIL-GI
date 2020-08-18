
The scripts used in the project "The evolutionary dynamics of genetic incompatibilities introduced by duplicated genes in Arabidopsis thaliana"


## step 1: get the SNP markers based on the whole genome resequencing data from the parents

###  map reads and call SNPs using shore pipeline 
http://shore.sourceforge.net/wiki/

shore import -v Fastq -Q illumina -x -y -o --rplot

shore mapflowcell -f -i -n 10 -g 7 -c 30 -p  --rplot

shore correct4pe -l  -x 150 -e 1 -p

shore merge -m -o  -p

shore consensus -n -f -o -i -g 5 -a -r

#Example:

shore import -v Fastq -Q illumina -x 357_F_CTTGTA_L002_R1.fastq -y 357_F_CTTGTA_L002_R2.fastq -o /flowc/ --rplot

shore mapflowcell -f flowc/ -i TAIR10_chr_all.fas.shore -n 10 -g 7 -c 30 -p  --rplot

shore correct4pe -l flowc/1 -x 260 -e 1 -p

shore merge -m map.list.1.gz,map.list.2.gz -o A -p

shore consensus -n A -f TAIR10_chr_all.fas.shore -o A/consensus -i A/map.list.gz -g 5 -a bin/shore/shore_current/shore/Analysis/scoring_matrices/scoring_matrix_hom.txt -v -r


## step 2: map reads of AMPRIL population RAD-seq. 
### import reads
shore import -v Fastq -a genomic -Q sanger -r barcodes.txt -x run2.lane4.reads_1.fastq.gz run2.lane4.reads_2fastq.gz -o flowcell --rplot -h 1

### map reads from whole sequencing run of pooled samples
shore mapflowcell -f flowcell -i TAIR10_chr_all.fas.shore -n 7 -g 4 -c 40 -p --rplot

#read number statistics: perl read_number_per_indi.pl -d dir -o output

#### Correcting for paired-end information for each sample
shore correct4pe -l sample_RAD-F198/ -x 250 -p -d

#alternatively, run perl correct_pe.pl -b barcode.txt

#### shore merge 
shore merge -m sample_RAD-F150/ -o ../../AlignmentFolder_sample_RAD-F150 -p

#alternatively, run perl merge_map_list.pl -b barcode.txt


## step 3: genotying using HMM
### 1) get allele count at maker positons from map.list files generated by shore, using the script allele_count_from_map_list2.pl

#e.g:

perl allele_count_from_map_list2.pl -n name -m map.list -p marker -r restrictions_sites -o output -g beginning and end cut

#alternatively, using the script allele_count_from_map_list_batch_submit.pl to submit batch jobs

#output files, e.g: allele_count.population_corrected.txt

#plant_id:189

1   83  M   0   0   0   0   0   T C T T

1   92  M   0   0   0   0   0   C C A A

1   110 M   0   0   0   0   0   T G G G

1   123 M   0   0   0   0   0   T C C C


### 2) genotyping
run basecall_hmm_batch.pl which runs the base_caller.jar, hmm.jar, hmm_output_to_marker.pl, create_break_file.pl

#Example:

java -Xmx1G -jar hmm/base_caller.jar -r allele_count.population_corrected.TE_filtered.txt -o allele_count.population_corrected.TE_filtered.txt.base_call.txt -n multi -m;

#hmm

java -jar hmm/hmm.jar -r allele_count.population_corrected.TE_filtered.txt.base_call.txt -o allele_count.population_corrected.TE_filtered.txt.base_call.hmm.txt -t multi -z hmm/hmm_profile

#conversion

perl hmm_output_to_marker.pl -n barcode.txt -m allele_count.population_corrected.TE_filtered.txt -p allele_count.population_corrected.TE_filtered.txt.base_call.hmm.txt -o allele_count.population_corrected.TE_filtered.txt.base_call.hmm_marker.txt

#create break file

perl create_break_file.pl -m allele_count.population_corrected.TE_filtered.txt.base_call.hmm_marker.txt -c tair10.chr.size.txt -o allele_count.population_corrected.TE_filtered.txt.base_call.hmm_marker_break.txt
 


# GWAS
1) run gwas, e.g.: with the phenotype of presence of non-functionalized HPA copy 

Rscript gwas.LM.multipro.r gwas.data.Rdata HPA.phenotype.01.txt HPA 20

2) draw the mahattan plot with gwas output.

Rscript manhattan_plot.r HPA.gwas.result.txt HPA.gwas
