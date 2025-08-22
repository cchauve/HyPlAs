
from glob import glob 
import random

refpath=  config["refpath"] if "refpath" in config else "data/ref"
refs = [x[:x.rfind(".fasta")] for x in glob ("*.fasta", root_dir = refpath)]

tksm_path="/home/fatih/workspace/tksm-testing/tksm/build/bin/tksm"
analysispath="data/tksm"
srpath = "data/SR"
lrpath = "data/LR"
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
                header = l[1:]
                continue
            yield (name, header,"".join(seq))
            seq = list()
            name = l[1:].split(" ")[0]
            header = l[1:]

        else:
            seq.append(l)
    yield (name, header,"".join(seq))
def get_chr_depth(wc):
    return random.randint(30, 50)
def get_plasmid_depth(wc):
    return random.randint(30, 50) * random.randint(1, 50)

rule all:
    input:
        [f"{srpath}/{ref}_R1.fq" for ref in refs],
        [f"{lrpath}/{ref}.fq.gz" for ref in refs]


rule index_faidx:
    input:
        "{sample}.fasta"
    output:
        "{sample}.fasta.fai"
    shell:
        "samtools faidx {input}"

rule tksm_random_wgs_chr:
    input:
        tksm_bin = tksm_path,
        ref_index= analysispath + "/{ref}_chr.fasta.fai",
        ref = analysispath + "/{ref}_chr.fasta"
    output:
        analysispath + "/{ref}_chr.mdf"
    params:
        depth = lambda wc: get_chr_depth(wc) 
    shell:
        "{input.tksm_bin}"
        " random-wgs"
        " -r {input.ref}"
        " -o {output}"
        " --frag-len-dist \"normal 350 50\""
        " --circular"
        " --depth {params.depth}"

rule tksm_random_wgs_plas:
    input:
        tksm_bin = tksm_path,
        ref_index= analysispath + "/{ref}_plasmid.fasta.fai",
        ref = analysispath + "/{ref}_plasmid.fasta"
    output:
        analysispath + "/{ref}_plasmid.mdf"
    params:
        depth = lambda wc : get_plasmid_depth(wc) 
    shell:
        "{input.tksm_bin}"
        " random-wgs"
        " -r {input.ref}"
        " -o {output}"
        " --frag-len-dist \"normal 350 50\""
        " --circular"
        " --depth {params.depth}"


rule split_plasmid_and_chr_from_ref:
    input:
        ref = refpath + "/{ref}.fasta"
    output:
        chr = analysispath + "/{ref}_chr.fasta",
        plasmid = analysispath + "/{ref}_plasmid.fasta"
    run:
        with open(output.chr, "w") as chr_file, open(output.plasmid, "w") as plasmid_file:
            for i, (name, header, seq) in enumerate(generate_fasta(input.ref)):
                if i == 0:
                    print(f">{name} {header}\n{seq}", file = chr_file)
                else:
                    print(f">{name} {header}\n{seq}", file = plasmid_file)

rule tksm_sequence_perfect:
    input:
        tksm_bin = tksm_path,
        ref= refpath + "/{ref}.fasta",
        plas_mdf= analysispath + "/{ref}_plasmid.mdf",
        chr_mdf = analysispath + "/{ref}_chr.mdf"
    output:
        perf = analysispath + "/{ref}_fragments.fasta",
    threads:
        32
    shell:
        "cat {input.plas_mdf} {input.chr_mdf} | "
        "{input.tksm_bin}"
        " sequence"
        " -r {input.ref}"
        " -i /dev/stdin"
        " --perfect {output.perf}"
        " -O fasta"
        " -t {threads}"

rule art_short_reads_from_perfect:
    input:
        perf = analysispath + "/{ref}_fragments.fasta",
    output:
        r1 = srpath + "/{ref}_R1.fq",
        r2 = srpath + "/{ref}_R2.fq",
    params:
        param= " -ss HS20 -amp -mp  -na  -l 100 -f 1  -s 10 -qs 10 -qs2 10",
        prefix= "data/SR/{ref}_R",
    shell:
        "art_illumina"
        " {params.param} "
        " -i {input.perf} -o {params.prefix}"

rule badread_for_lr:
    input:
        ref = refpath + "/{ref}.fasta"
    output:
        r = lrpath + "/{ref}.fq.gz"
    threads:
        8
    params:
        depth = lambda wc: f"{random.randint(25,35)}x"
    shell:
        "badread simulate"
        " --reference {input.ref}"
        " --small_plasmid_bias"
        " --quantity {params.depth}"
        " | pigz > {output.r}"
