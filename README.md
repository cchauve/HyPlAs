
## Installation

### Build hyplass binaries
```
conda create -f environment.yml
conda activate hyplass-env
make
```

Or with virtualenv
```
python -m venv hyplass_env
source hyplass_env/bin/activate
python3 build.py hyplass_env
source module_load.sh #For cedar. should be installed if not available
```

### Install Platon database (https://github.com/oschwengers/platon)
```
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz
#Move to a suitable location
```

## Usage

```
python src/hyplass.py --platon-db db -s sr_*.fastq -l lr.fastq.gz -o hyplass-out/ -t threads -p prop_rounds      
```
