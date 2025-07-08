
# HyPlas
Hyplas is an automated plasmid discovery method works on short- long- hybrid sequencing data. 
It incorporates plasmid classification tools (Such as platon) on short-read assembled contigs to aid plasmidic long-read selection and performs hybrid assembly.

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
## Usage
```
python src/hyplass.py --platon-db db -s sr_*.fastq -l lr.fastq.gz -o hyplass-out/ -t threads -p prop_rounds      
```
### Input explanation
 - --platon-db: Database for the platon
 ```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
#Move to a suitable location
```
- -s space separated short read files
- -l long reads file (required to be gzipped)
- -o output folder
- -p number of long read recovery rounds to be executed (Recommend 2 rounds)

### Output explanation
- assembly.final.fasta: 
	- Fully circularized plasmid sequences
- unicycler_sr (directory):
	- - Output of the short-read-only assembly by Unicycler
- classify (directory):
	- classify/result.log: Output log of the platon
	- classify/result.json: Classification details of the short-read assembly contigs in json format
	- classify/result.tsv: Classification details of the short-read assembly contigs in tsv format
	- classify/result_p.tsv: List of contigs predicted as plasmidic or chromosomal
- lr2assembly.gaf: Graph alignment of long reads to the short-read-only assembly.
- plasmid_long_reads/plasmid.fastq.gz: Plasmid predicted long-reads
- plasmid_long_reads/plasmid.fastq.gz: Unknown origin long-reads
- prop_lr/ (directory):
	- prop_lr/lr.round.[0-9]+.paf: Mappings of the known plasmid long-reads to unknown long-reads
	- prop_lr/lr.round.[0-9]+.fastq.gz: Plasmidic long-read sequences recovered in this round
- unicycler_lr (directory):
	- Output of the long hybrid assembly by Unicycler
