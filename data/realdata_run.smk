
import os
import sys

DEBUG=True
if len(config) == 0:
    configfile: 'config.yaml'



plasmid_reference_dir  = "ncbi-reference-genomes/eskapee/plasmids"

outpath = config['outpath']
rawdata = config['rawdata']
firstN = -1 if "firstN" not in config else config["firstN"]
PLASMID_SAVE_ITER = 2

plassembler_db_path = f"{outpath}/plassembler-db/"

def generate_fasta(fasta):
    name = "NULL"
    seq = list()
    if fasta.endswith(".gz"):
        f = gzip.open(fasta, "rt")
    else:
        f = open(fasta, "r")
    header = "NULL"
    for _, l in enumerate(f):
        l = l.rstrip("\n")
        if l[0] == ">":
            if len(seq) == 0:
                name = l[1:].split(" ")[0]
                header = " ".join(l.split(" ")[1:])
                continue
            yield (name, header,"".join(seq))
            seq = list()
            name = l[1:].split(" ")[0]
            header =  " ".join(l.split(" ")[1:])

        else:
            seq.append(l)
    yield (name, header,"".join(seq))

def absolute_file_paths(directory):
    path = os.path.abspath(directory)
    return [entry.path for entry in os.scandir(path) if entry.is_file()]

def load_from_tsv(tsv_file):
    import pandas as pd
    blacklisted = set()
    for line in open("blacklist.txt"):
        blacklisted.add(line.strip())
    df = pd.read_csv(tsv_file, sep='\t')
    if 'samples' not in config:
        config['samples'] = {}
    for i, row in df.iterrows():
        nom = row['Assembly BioSample Accession']
        if nom in blacklisted:
            continue
        if pd.isna( row['SRA Accession']):
            continue
        p2a = {x.split(":")[0]:x.split(":")[1] for x in row['SRA Accession'].split("; ")[1:]}

        if 'OXFORD_NANOPORE' not in p2a:
            continue
        if 'ILLUMINA' not in p2a:
            continue
        if not os.path.isfile(f"{rawdata}/ILLUMINA/{nom}/{p2a['ILLUMINA']}_2.fastq.gz"):
            continue
        if not os.path.isfile(f"{plasmid_reference_dir}/{nom}.fasta"):
            continue
        config['samples'][nom] = {
            'lr': f"{rawdata}/OXFORD_NANOPORE/{nom}/{p2a['OXFORD_NANOPORE']}_1.fastq.gz",
            'sr': [f"{rawdata}/ILLUMINA/{nom}/{p2a['ILLUMINA']}_1.fastq.gz",
                    f"{rawdata}/ILLUMINA/{nom}/{p2a['ILLUMINA']}_2.fastq.gz"
                   ],
            'sr-path' : f"{rawdata}/ILLUMINA/{nom}/",
            'lr-path' : f"{rawdata}/OXFORD_NANOPORE/{nom}/"
        }
        if firstN > 0 and i > firstN:
            break


if "tsv" in config:
    load_from_tsv(config['tsv'])

def get_sample_lr(wildcards):
    return config["samples"][wildcards.sample]["lr"]
def get_sample_sr(wildcards):
    return config["samples"][wildcards.sample]["sr"]
def get_sample_sr1(wildcards):
    return sorted(absolute_file_paths(config["samples"][wildcards.sample]["sr-path"]))[0]
def get_sample_sr2(wildcards):
    return sorted(absolute_file_paths(config["samples"][wildcards.sample]["sr-path"]))[1]
def get_sample_sr(wildcards):
    return absolute_file_paths(config["samples"][wildcards.sample]["sr-path"])[:2]

def get_sample_lr(wildcards):
    return absolute_file_paths(config["samples"][wildcards.sample]["lr-path"])[0]

if DEBUG:
    def pipe(X):
        return X
    def temp(X):
        return X

rule all:
    input:
        expand(f"{outpath}/{{sample}}/hyplass_b_cov/plasmid_selected_i{{iter}}.coverage.tsv", sample=config["samples"], iter=[0,1,2]),


rule all_plassembler:
    input:
        expand(f"{outpath}/{{sample}}/plassembler_flye/plassembler_plasmids.fasta",sample=config["samples"]),

rule all_stats:
    input:
        expand(f"{outpath}/{{sample}}/qc/toplasmids.coverage.tsv",sample=config["samples"]),
        expand(f"{outpath}/{{sample}}/plassembler{{flye}}-qc/toplasmids.coverage.tsv",sample=config["samples"], flye=["","_flye"]),




rule hyplas_rule:
    input:
        src = "src/main.py",
        lr = f"{outpath}/{{sample}}/qc/trim.lr.fastq.gz",
        sr = [f"{outpath}/{{sample}}/qc/{{sample}}.sr1.fastq.gz",
f"{outpath}/{{sample}}/qc/{{sample}}.sr2.fastq.gz"],
    output:
        dir=directory(f"{outpath}/{{sample}}/hyplass_b/"),
        fa=[f"{outpath}/{{sample}}/hyplass_b/assembly.final.it{it}.fasta" for it in range(3)],
    threads:
        32
    params:
        iter=2,
        db="../../annotation/db"
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 - 1,
        mem_mb=1024*88, 
    shell:
        "python {input.src}"
        " -l {input.lr}"
        " -o {output.dir}"
        " --platon-db {params.db}"
        " -p {params.iter}"
        " -t {threads}"
        " -s {input.sr}"
        " --force"




rule Plassembler_flye:
    input:
        lr = get_sample_lr,
#        sr1 = get_sample_sr1,
#        sr2 = get_sample_sr2,
        sr = get_sample_sr,
        db = f"{plassembler_db_path}"
    output:
        out_dir = temp(directory(f"{outpath}/{{sample}}/plassembler_flye")),
        fasta = f"{outpath}/{{sample}}/plassembler_flye/plassembler_plasmids.fasta",
        summary = f"{outpath}/{{sample}}/plassembler_flye/plassembler_summary.tsv",
        plasmid_fastq = f"{outpath}/{{sample}}/plassembler_flye/plasmid_fastqs/plasmids_long.fastq",
    params:
        par = "-c 1000000",
        sr_param = lambda wildcards: " ".join([f"-{i} {x}" for i,x in enumerate(get_sample_sr(wildcards), start=1)])
    threads:
        64
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 - 1,
        mem_mb=1024*700, 
    benchmark:
        f"{outpath}/{{sample}}/plassembler_flye.benchmark.tsv"
    shell:
        "plassembler run {params.sr_param}"
        " -l {input.lr}"
        " -d {input.db}"
        " -o {output.out_dir} --force"
        " -t {threads}"
        " {params.par}"
        " --keep_fastqs"



rule Plassembler:
    input:
        lr = get_sample_lr,
        sr = get_sample_sr,
        db = f"{plassembler_db_path}"
    output:
        out_dir = temp(directory(f"{outpath}/{{sample}}/plassembler")),
        fasta = f"{outpath}/{{sample}}/plassembler/plassembler_plasmids.fasta",
        summary = f"{outpath}/{{sample}}/plassembler/plassembler_summary.tsv",
        plasmid_fastq = f"{outpath}/{{sample}}/plassembler/plasmid_fastqs/plasmids_long.fastq",
    params:
        par = "-c 1000000 --use_raven",
        sr_param = lambda wildcards: " ".join([f"-{i} {x}" for i,x in enumerate(get_sample_sr(wildcards), start=1)])
    threads:
        64
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 -1,
        mem_mb=1024*700, 
    benchmark:
        f"{outpath}/{{sample}}/plassembler.benchmark.tsv"
    shell:
        "plassembler run {params.sr_param}"
        " -l {input.lr}"
        " -d {input.db}"
        " -o {output.out_dir} --force"
        " -t {threads}"
        " {params.par}"
        " --keep_fastqs"



rule plassembler_download_db:
    output:
        directory(f"{plassembler_db_path}"),
    shell:
        "plassembler download -d {output}"



rule fastp:
    input:
        sr1 = get_sample_sr1,
        sr2 = get_sample_sr2,
    output:
        sr1 = f"{outpath}/{{sample}}/qc/{{sample}}.sr1.fastq.gz",
        sr2 = f"{outpath}/{{sample}}/qc/{{sample}}.sr2.fastq.gz",
        unpaired = f"{outpath}/{{sample}}/qc/{{sample}}.unpaired.fastq.gz",
    params:
        par = "",
    threads:
        16
    resources:
        runtime=60,
        mem_mb=1024*8,

    shell:
        "fastp --in1 {input.sr1} --in2 {input.sr2}"
        " --out1 {output.sr1} --out2 {output.sr2}"
        " --unpaired1 {output.unpaired}"
        " --unpaired2 {output.unpaired}"
        " --thread {threads}"
        " {params.par}"


rule chopper_lr:
    input:
        get_sample_lr,
    output:
        lr = f"{outpath}/{{sample}}/qc/trim.lr.fastq.gz",
    params:
        minqual=9,
        minlen=500,
        headcrop=75,
        tailcrop=75,
    threads:
        32
    shell:
        "chopper"
        " -q {params.minqual}"
        " --threads {threads}"
        " -l {params.minlen}"
        " --headcrop {params.headcrop}"
        " --tailcrop {params.tailcrop}"
        " --input {input}" 
        " | pigz > {output.lr}"

rule bwa_index:
    input:
        f"{outpath}/{{genome}}.fasta"
    output:
        idx=multiext(f"{outpath}/{{genome}}.{{alg}}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{genome}.{alg}.log",
    params:
        extra="",
    wrapper:
        "v3.13.8/bio/bwa/index"



rule samtoolscoverageforgt:
    input:
        bam=f"{outpath}/{{file}}.bam",
    output:
        tsv=f"{outpath}/{{file}}.coverage.tsv",
    shell:
        "samtools coverage {input.bam} > {output.tsv}"

rule minimapPlassembler2GT:
    input:
        filter_bin = "src/filter_plasmid_sam.py", 
        ref= f"{plasmid_reference_dir}/{{sample}}.fasta",
        ctgs = f"{outpath}/{{sample}}/plassembler{{is_flye}}/plasmid_fastqs/plasmids_long.fastq",
    output:
        sam=f"{outpath}/{{sample}}/plassembler{{is_flye}}-qc/toplasmids.sam",
        samf=f"{outpath}/{{sample}}/plassembler{{is_flye}}-qc/toplasmids_filtered.sam",
        bam=f"{outpath}/{{sample}}/plassembler{{is_flye}}-qc/toplasmids.bam",
    threads:
        64
    wildcard_constraints: is_flye=".*"
    shell:
        "minimap2"
        " {input.ref}"
        " {input.ctgs}"
        " -t {threads}"
        " -x map-ont"
        " -a -o {output.sam}"
        " && python {input.filter_bin} {output.sam} {output.samf}"
        " && samtools sort {output.samf} -@ 8 -o {output.bam}"
rule minimapAll2GT:
    input:
        filter_bin = "src/filter_plasmid_sam.py", 
        ref= f"{plasmid_reference_dir}/{{sample}}.fasta",
        ctgs = f"{outpath}/{{sample}}/qc/trim.lr.fastq.gz",
    output:
        sam=f"{outpath}/{{sample}}/qc/toplasmids.sam",
        samf=f"{outpath}/{{sample}}/qc/toplasmids_filtered.sam",
        bam=f"{outpath}/{{sample}}/qc/toplasmids.bam",
    threads:
        64
    shell:
        "minimap2"
        " {input.ref}"
        " {input.ctgs}"
        " -t {threads}"
        " -x map-ont"
        " -a -o {output.sam}"
        " && python {input.filter_bin} {output.sam} {output.samf}"
        " && samtools sort {output.samf} -@ 8 -o {output.bam}"


from glob import glob
def get_extension_lrs(wc):
    fmt =f"{outpath}/{wc.sample}/hyplass_b/prop_lr/lr.round.{{}}.fastq.gz" 
    return [fmt.format(str(i)) for i in range(int(wc.iter)) if os.path.isfile( fmt.format(str(i))) ]
rule minimap2GT_single_script:
    input:
        filter_bin = "src/filter_plasmid_sam.py", 
        ref=f"{plasmid_reference_dir}/{{sample}}.fasta",
        lr_plasmid = f"{outpath}/{{sample}}/hyplass_b/plasmid_long_reads/plasmid.fastq.gz",
        #lr_plasmid2 = lambda w:[f"{outpath}/{w.sample}/hyplass_b/prop_lr/lr.round.{i}.fastq.gz" for i in range(int(w.iter))],
        dir=directory(f"{outpath}/{{sample}}/hyplass_b/"),
    output:
        sam=temp(f"{outpath}/{{sample}}/hyplass_b_cov/plasmid_selected_i{{iter}}.sam"),
        samf=temp(f"{outpath}/{{sample}}/hyplass_b_cov/plasmid_selected_i{{iter}}_filtered.sam"),
        bam=f"{outpath}/{{sample}}/hyplass_b_cov/plasmid_selected_i{{iter}}.bam",
    threads:
        32
    params:
        extension=get_extension_lrs
    shell:
        "minimap2"
        " {input.ref}"
        " {input.lr_plasmid}"
        " {params.extension}"
        " -t {threads}"
        " -x map-ont"
        " -a -o {output.sam}"
        " && python {input.filter_bin} {output.sam} {output.samf}"
        " && samtools sort {output.samf} -@ 8 -o {output.bam}"






rule bam2fq:
    input:
        f"{{sample}}.bam",
    output:
        f"{{sample}}.fq.gz",
    shell:
        "samtools fastq"
        " -1 /dev/null -2 /dev/null -s /dev/null"
        " -n {input}"
        " | pigz > {output}"


