
# HyPlas
HyplAs is a tool aimed at assembling plasmids from hybid short-read and long-read sequencing data for bacerial isolates.
HyPlAs main novlty is to incorporate a plasmid classification tools (Such as platon) on short-read assembled contigs to aid plasmidic long-read selection and performs hybrid plasmids assembly.
HyPlAs has been desiged to work with single genome sequencing data, and has not been tested on metagenomics data.  

## Installation #TODO add setup.py

### Bioconda #TODO
### Build hyplass binaries
```
conda create -f environment.yml
conda activate hyplass-env
make
```

Or with virtualenv 
```
source module_load.sh #For cedar. should be installed if not available
python -m venv hyplass_env
source hyplass_env/bin/activate
python3 build.py hyplass_env
```

## Overview

HyPlAs  is a pipeline combinig existing tools and specific Python scripts and C++ programs. HyPlAs is composed of
the following steps (see figure below): 
1. Reads preprocessing;  
	1.a. Short reads are preprocessed with <a href="https://github.com/OpenGene/fastp">fastp</a>,  
    	1.b. Long reads are preprocessed with <a href="https://github.com/wdecoster/chopper">chopper</a>;  
3. Short reads are assembled using <a href="https://github.com/rrwick/Unicycler">Unicycler</a>;
4. The detection of putative plasmidic long reads is done in four stages:  
	3.a. the plasmid contigs classification tool <a href="https://github.com/oschwengers/platon">Platon</a> is used to detect plasmidic short-read contigs,  
	3.b. long reads are mapped to the assembly graph using <a href="https://github.com/lh3/minigraph">minigraph</a>,  
   	3.c. long-read mapping to short-read contigs and platon results are used to select an initial set of putative plasmidic long reads,  
   	3.d. the set of putative plasmidic long reads is augmented by iteratively detecting overlapping long reads;  
5. The full short-read assembly graph generated in step 2 is refined with the plasmidic long reads selected during step 3, using Unicycler.

![HyPlAs](resources/HyPlAs_pipeline.png?raw=true)

## Usage
```
python src/hyplas.py --platon-db db -s sr_*.fastq -l lr.fastq.gz -o hyplass-out/ -t threads -p prop_rounds      
```
### Input
 - --platon-db: Database used by Platon (<a href="https://zenodo.org/record/4066768/files/db.tar.gz">https://zenodo.org/record/4066768/files/db.tar.gz</a>)
 ```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
# Move the database to a suitable location
```
- -s space separated short read files
- -l long reads file (required to be gzipped)
- -o output folder
- -p number of long read recovery rounds to be executed (Recommend 2 rounds)

### Output
HyPlAs creates in the output folder the following files and directories:  
- plasmids.final.fasta:   
	- assembled plasmids, in FASTA format;  
- unicycler_sr (directory):  
	- short-read-only assembly by Unicycler;  
- classify (directory):  
	- classify/result.log: Platon log file,  
	- classify/result.json: Platon classification details of the short-read assembly contigs in json format,  
	- classify/result.tsv: Platon classification details of the short-read assembly contigs in tsv format,  
	- classify/result_p.tsv: List of contigs predicted by Platon as plasmidic or chromosomal;  
- lr2assembly.gaf: Graph alignment of long reads to the short-read-only assembly cotigs;  
- plasmid_long_reads/plasmid.fastq.gz: Long-reads classified as plasmidic by HyPlAs, in FASTQ format;  
- prop_lr/ (directory):
	- prop_lr/lr.round.[0-9]+.paf: Mappings of the known plasmid long-reads to unknown long-reads, the integer suffix indicates the iteration round of plasmidic long-reads augmentation (step 3.d),  
	- prop_lr/lr.round.[0-9]+.fastq.gz: Plasmidic long-read sequences recovered in augmentation round X (X in [0-9]);  
- unicycler_lr (directory):  
	- Output of the final Unicycler hybrid assembly using the short-read--only assembly and the predicted plasmdic long reads.  
