###########################################################################
###																		
###   This Makefile takes in input the link of a downloadble resulting  
###   alignment from a ChIP-Seq experiment in BAM format and gives as   
###	  output a list of regulatory scores between the TF under study		
###   and each gene annotated in the genome of the species of interest  
###																		
###	  The pipeline:													    
###	  1. ChIP-Seq BAM alignment											
###		 	|__ 2. conversion of file format to WIG 					
###				   		|__ 3. execution of TIP algorithm on WIG		
###						   		|__ 4. assignment of regulatory scores  
###																		
###	  With respect to the following project's file system:				
###																		
###				fucoup_directed/										
###				|__data/												
###				|   |-- annotations/									
###				|   |-- db/												
###				|   |   |-- regnet/										
###				|   |   |__ trrust/										
###				|   |-- evidences/										
###				|   |   |__ encode/										
###				|   |       |-- human/									
###				|   |       |   |__ chip_seq/							
###				|   |       |       |__ scores/							
###				|   |       |           |__ 01/							
###				|__local/												
###					|-- bin/																		
###																		
###   From fucoup_directed/:											
###   - RUN `export PATH=$PATH:$PWD/local/bin`							
###	  - RUN `sudo chmod 755 local/src/*`								
###   From fucoup_directed/local/bin:									
###   - RUN `for i in `ls ../src`; 										
###				do ln -s ../src/$i $(echo $i | cut -d '.' -f1); done`	
###   From funcoup_directed/evidences/encode/							
###	  - Run `make` (specify parameters' content if necessary)			

# Define the shell interpreter
shell := /bin/bash -O globstar

# User-defined variables
organism := human
ref_genome := GRCh38
anno_version := 29
q_value := 05

# Define path to directories
project := ~/funcoup_directed/
encode := $(shell echo $$PWD/)
source := $(encode)$(organism)/
chip_seq := $(source)chip_seq/
exp := $(chip_seq)experiments/
scores := $(chip_seq)scores/$(q_value)/

# Define path to files
files := $(source)files.txt
metadata := $(source)metadata_bam.no_rep.tsv
annotation := $(project)data/annotations/$(organism)/$(ref_genome)_GENCODE_V$(anno_version).TIP20000.gtf

################# GOAL ###################

# Use `addsuffix` to expand a list of accession file names
# and define the objective of the pipeline
content := $(shell tail -n +2 $(metadata) | cut -f1 | sort -R)
all: $(addsuffix .tip, $(content))

################# STEP 5 #################
## Execute the TIP algorithm on compressed WIG formatted chromosomes
## to extract TF-gene regulatory scores of interaction. Sort scores and isolate
## just those interactions whose qvalue is higher than a pre-set threshold.
## Copy results and metadata on `scores` folder.
## Create a compressed tarball with all the WIG formatted chromosomes.
%.tip: %.zig
	@tip_on_wiggle $(exp) $* $(annotation) -og $(organism) -g $(ref_genome) -r
	@tail -n +2 $(exp)$*/score_refseq.tsv | awk '$$2>0 && $$3>0 && $$5<0.$(q_value)' | LC_ALL=C sort -grk 2 > $(exp)$*/score.$(q_value).sort.tsv
	@mkdir $(scores)$*; cp $(exp)$*/score.$(q_value).sort.tsv $(exp)$*/metadata.tsv $(scores)$*
	@tar -czvf $*.tar.gz -C $(exp)$* zig; mv $*.tar.gz $(exp)$*/zig.tar.gz; rm -r $(exp)$*/zig

################# STEP 4 #################
## Once the WIG file is ready, split the file in chromosomes
## and compress the information in fewer bites.
%.zig: %.wig
	@compress_wiggle $(exp) $* -og $(organism) -r

################# STEP 3 #################
## Convert the BAM file to WIG format and remove BAM.  
%.wig: %.metadata
	@bam_to_wiggle $(exp)$*/$*.bam --outfile=$(exp)$*/$*.wig
	@rm $(exp)$*/$*.bam $(exp)$*/$*.bam.bai &> /dev/null

################# STEP 2 #################
## Once the BAM file is downloaded, save specific metadata
## in its dedicated folder.
%.metadata: %.bam
	@head -n +1 $(metadata) > $(exp)$*/metadata.tsv
	@grep $* $(metadata) >> $(exp)$*/metadata.tsv

################# STEP 1 #################
## If the list with links and the metadata file is available,
## create a dedicated folder and download a BAM file.
%.bam: $(files) $(metadata)
	@echo $*
	@mkdir -p $(exp)$*/zig
	@wget -q $$(grep $* $(files)) -P $(exp)$*

.PRECIOUS: %.metadata %.zig
