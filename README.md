# Building directed regulatory FunCoup networks from HT-data 

This work is part requirement for the MSc in Bioinformatics at University of Bologna (IT). 
The project is developped at Science For Life Laboratory, Stockholm (SE)

## Overview

If we highlight natural interactions between genes or proteins by integrating numerous and different types of evidences, we get to understand the complex interplay between biomolecules more accurately [1]. Also, using small-scale experiments to find out these interactions, high quality (HQ) data is delivered but this strategy lacks in coverage and flexibility towards the aforementioned variability in datasets and evidences. We started to represent genes or proteins with links that depict respective interactions. These links are supported by the available high-throughput data (HT-data) that are analysed with the framework in use, which is most of the time a statistical approach and less frequently a tool based on a machine-learning scheme. What finally widens the stringent definition of an interaction and introduces a more general term, that is functional association or coupling (FC), is the integration of multiple dataset and/or different evidence types.

The power of evidences’ integration is taken into account by various methods, such as the FunCoup functional association network [2]. In particular, FunCoup employs a redundancy weighted Bayesian integration framework to combine functional association HT-data in order to build comprehensive biological networks. To assign a log likelihood ratio score to different types of evidence, it uses HQ gold standard data, presenting known functionally associated gene pairs. Thanks to this, those types of evidence can then be integrated. Currently, all such evidences are undirected in their nature and, by exploiting online accessible HT-data, this project aims at building new networks with directional links. The mentioned data include regulatory links, i.e. links from transcription factors (TFs) to their target genes.


## Prerequisites
Python 3.8

## Installation
1. Install direnv:
    * run `brew install direnv (macOS)`
    * run `sudo apt-get install direnv (linux)`
2. Install venv (inspired by https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/):
    * run `sudo apt-get install python3-venv`
3. Setup of .bash_profile and .config/direnv (inspired by https://github.com/direnv/direnv/wiki/Python):
    * run `cat direnv_allow.txt >> ~/.bashrc; echo "export PATH=$PATH:./local/bin" > "./.envrc"`
    * run `cat venv_allow.txt >> ~/.bashrc; echo "layout python-venv python3.5"  >> "./.envrc"`
    * run `cat direnvrc_config.txt > ~/.config/direnv/direnvrc`
    * run `source ~/.bashrc; direnv allow ./`
4. Install python packages:
    * run `pip install -r requirements.txt`
5. Install from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/) the program in `local/bin`:
    * `wigToBigWig`

## Project File System 
```bash
__data/
|   |-- annotations/
|   |   |__ human/
|   |-- db/
|   |   |-- fc4.1/
|   |   |   |__ human/
|   |   |-- regnet/
|   |   |   |__ human/
|   |   |__ trrust/
|   |       |__ human/
|   |-- evidences/
|   |   |__ encode/
|   |       |__ human/
|   |           |__ chip_seq/
|   |               |-- experiments/
|   |               |-- info/
|   |               |__ scores/
|   |                   |-- 01/
|   |                   |__ 05/
|   |__ stat/
|       |__ human/
|           |-- regnet/
|           |__ trrust/
|__local/
    |-- bin/
    |   |__ symb_links_to_src
    |-- config/
    |-- docs/
    |-- makefiles/
    |__ src
        |__ python3/bash scripts
    
```

## Key files
Each program has its own header with detailed documentation and examples of usage. Each function in every script is introduced by some rows of explanation. 

### The network in format gene1 gene2 LLR
- `data/stat/human/trrust/05/weighted_funcoup.ensemblID.fl-False.tsv`: FC_directed trained on TRRUSTv2
- `data/stat/human/regnet/05/weighted_funcoup.ensemblID.fl-False.tsv`: FC_directed trained on RegNetwork

### The tarball with GS and evidences
- `fc_directed.tar.gz`: Ready to be shared.

### Build Gold Standards, preprocess annotation files
- `local/src/Dataset.py`: Builds positive/negative gold standards[3-4].
- `local/docs/human_annotation_files.ipynb`: Pre-processes raw original annotation files.
- `local/src/encode_api.py`: Gets in touch with ENCODE-API[5] to download ChIP-Seq metadata of BAM files. 

### Download, convert and extract scores from ChIP-Seq BAM files
- `local/makefiles/chip_seq.mk`: Takes a list of BAM links in input and exectues in order:
    1. `local/src/bam_to_wiggle.py`: Converts BAM files to Wiggle file format.
    2. `local/src/compress_wiggle.py`: Reworks the content of a Wiggle file to save space.
    3. `local/src/tip_on_wiggle.py`: Executes TIP[6] on a Wiggle file format.

### Compute some final statistics
- `local/src/db_overlap.py`: Computes overlap of nodes and links between databases'network. 
- `local/src/tip_statistics.py`: Computes some naïve statistics on TIP execution's results.

### Build the FC_directed network and compare it with FC4.1
- `local/src/funcoup_directed.py`: Computes and assign LLR to TF-gene links.
- `local/docs/FC41-FC_directed(trrust).ipynb`: Looks at overlap and coorellation with actual FC4.1 assignments.

- **RUN `results human` to execute db_overlap.py. tip_statistics.py and funcoup_directed.py in a shot.**

### Project documentation
- `local/docs/regulatory_links-report.pdf`: Reads a detailed overview of the project.
- `local/docs/regulatory_links-presentation.pdf`: Shows a detailed overview of the project.

## Key directories
In https://bitbucket.org/sonnhammergroup/directed_links/:

- `data/annotations/`: Contains genome annotation files. 
- `data/db/`: Contains raw and processed datasets.
- `data/evidences/encode/`: Contains ChIP-Seq related files.
- `data/evidences/encode/human/chip_seq/scores/`: Contains potential regulatory scores to ChIP-Seq experiments.
- `data/stat/`: Contains all the statistics done. 

Exclusively in Octa2:

- `data/evidences/encode/human/chip_seq/experiments/`: Contains 1468 experiments (see below for an example).
- `data/evidences/encode/human/chip_seq/nan/`: Contains 156 experiments which gave errors during TIP execution.

**Example of experiment stored in Octa2** --> `data/evidences/encode/human/chip_seq/experiments/ENCFF002PPW/`:

1. **metadata.tsv**: the metadata of ENCFF002PPW ChIP-Seq experiment with BHLHE40 as target TF;
2. **score_refseq.tsv**: the list of 19927 potential BHLHE40-genes regulatory interaction;
3. **score.05.sort.tsv**: the sorted list of significant BHLHE40-gene regulatory interactions (qvalue < 0.05)
4. **score.01.sort.tsv**: the sorted list of significant BHLHE40-gene regulatory interactions (qvalue < 0.01)
5. **weights.npy**: the BHLHE40 experimental binding profile, stored as numpy array in binary format. 
6. **zig.tar.gz**: the compressed archive of chromosomes' binding signals. RUN `tar xf zig.tar.gz` to extract the compressed archive into the current directory. Eventually RUN `local/src/tip_on_wiggle.py` to tune the size of the window centered into genes'TSS.  


## References
[1] Guala, D., Ogris, C., Müller, N., & Sonnhammer, E. L. L. (2019). Genome-wide functional association networks: background, data & state-of-the-art resources. Briefings in Bioinformatics, 00(February), 1–14. https://doi.org/10.1093/bib/bbz064

[2] Alexeyenko A, Sonnhammer EL. “Global networks of functional coupling in eukaryotes from comprehensive data integration.” Genome Res. 2009 Jun;19(6):1107-16 

[3] Vaquerizas, J. M., Kummerfeld, S. K., Teichmann, S. A., & Luscombe, N. M. (2009). A census of human transcription factors: Function, expression and evolution. Nature Reviews Genetics, 10(4), 252–263. https://doi.org/10.1038/nrg2538

[4] Davis, C. A., Hitz, B. C., Sloan, C. A., Chan, E. T., Davidson, J. M., Gabdank, I., ... Cherry, J. M. (2018). The Encyclopedia of DNA elements (ENCODE): Data portal update. Nucleic Acids Research, 46(D1), D794–D801. https://doi.org/10.1093/nar/gkx1081

[5] Han, H., Cho, J. W., Lee, S., Yun, A., Kim, H., Bae, D., ... Lee, I. (2018). TRRUST v2: An expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic Acids Research, 46(D1), D380–D386. https://doi.org/10.1093/nar/gkx1013

[6] Cheng, C., Min, R., & Gerstein, M. (2011). TIP: A probabilistic method for identifying transcription factor target genes from ChiP-seq binding profiles. Bioinformatics, 27(23), 3221– 3227. https://doi.org/10.1093/bioinformatics/btr552