
import pandas as pd
from glob import glob

outpath = config["outpath"] if "outpath" in config else "plassembler_benchmark_simulation_results"
plas_db_path = config["plassembler_db_path"] if "plassembler_db_path" in config else "plassembler-db"
hyplas_path = config["hyplas_path"] if "hyplas_path" in config else "hyplas"
def get_samples(sample_csv):
    df = pd.read_csv(sample_csv,header=None, usecols=[0])
    print(df)
    return df[0].to_list()
samples = refs = [x[1+ x.rfind("/"):x.rfind(".fasta")] for x in glob ("data/ref/*.fasta")] 

def get_sample_sr(wildcards):
    return [f"data/SR/{wildcards.sample}_R1.fq",f"data/SR/{wildcards.sample}_R2.fq"]
rule all:
    input:
        plassembler_rav = [outpath + f"/results/plassembler_raven/{sample}.fasta" for sample in samples],
        plassembler = [outpath + f"/results/plassembler_flye/{sample}.fasta" for sample in samples],
        hyplas = [outpath + f"/results/hyplas-{i}/{sample}.fasta" for sample in samples for i in range(3)],

rule plassembler_assembly_raven:
    input:
        sr = get_sample_sr,
        lr = "data/LR/{sample}.fq.gz",
    output:
        out_dir = directory(outpath + "/plassembler_raven/{sample}"),
        contigs = outpath + "/results/plassembler_raven/{sample}.fasta",
    params:
        db = plas_db_path,
        # par = "-c 1000000",
        par = "",
        sr_param = lambda wildcards: " ".join([f"-{i} {x}" for i,x in enumerate(get_sample_sr(wildcards), start=1)]),
    threads:
        32
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 - 1,
        mem_mb=1024*120, 
    shell:                
            "plassembler run {params.sr_param}"
            " -l {input.lr}"
            " -d {params.db}"
            " -o {output.out_dir}"
            " -t {threads}"
            " --use_raven"        
            " {params.par}"       
            " --keep_fastqs"
            " && ln -s  {output.out_dir}/plassembler_plasmids.fasta {output.contigs}"      

rule plassembler_assembly:
    input:
        sr = get_sample_sr,
        lr = "data/LR/{sample}.fq.gz",
    output:
        out_dir = directory(outpath + "/plassembler_flye/{sample}"),
        contigs = outpath + "/results/plassembler_flye/{sample}.fasta",
    params:
        db = plas_db_path,
        # par = "-c 1000000",
        par = "",
        sr_param = lambda wildcards: " ".join([f"-{i} {x}" for i,x in enumerate(get_sample_sr(wildcards), start=1)]),
    threads:
        32
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 - 1,
        mem_mb=1024*120, 
    shell:                
            "plassembler run {params.sr_param}"
            " -l {input.lr}"
            " -d {params.db}"
            " -o {output.out_dir}"
            " -t {threads}"        
            " {params.par}"       
            " --keep_fastqs"
            " && ln -s  {output.out_dir}/plassembler_plasmids.fasta {output.contigs}"         

rule fastp:           
    input: 
        sr1 = lambda wc : get_sample_sr(wc)[0],      
        sr2 = lambda wc : get_sample_sr(wc)[1],                            
    output:       
        sr1 = "data/SR/{sample}_fastp_R1.fq",
        sr2 = "data/SR/{sample}_fastp_R2.fq",
        unpaired = "data/SR/{sample}_fastp_UN.fq"
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
        "data/LR/{sample}.fq.gz"
    output:
        lr = "data/LR/{sample}.trim.fastq.gz",
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


rule hyplas_assembly:
    input:
        src = hyplas_path,
        sr = ["data/SR/{sample}_fastp_R1.fq","data/SR/{sample}_fastp_R2.fq"],
        lr = "data/LR/{sample}.trim.fastq.gz",
    output:
        contigs=[outpath + f"/hyplass/{{sample}}/assembly.final.it{i}.fasta" for i in range(3)],
        dir=directory(f"{outpath}/hyplass/{{sample}}/"),
    params:
        iter=2,
        db="../../annotation/db",
    threads:
        32
    resources:
        runtime=lambda wc, attempt: attempt * 6 * 60 - 1,
        mem_mb=1024*120, 
    shell:
        "python {input.src}"
            " -l {input.lr}"
            " -o {output.dir}"
            " --platon-db {params.db}"
            " -p {params.iter}"
            " -t {threads}"
            " -s {input.sr}"
            " --force"

rule link_hyplas:
    input:
        contig=outpath + "/hyplass/{sample}/assembly.final.it{it}.fasta",
    output:
        outpath + "/results/hyplas-{it}/{sample}.fasta"
    shell:
        "ln -s {input.contig} {output}"
    
