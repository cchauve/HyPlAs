import argparse
import os
import sys
from contextlib import contextmanager
import subprocess
import tempfile
import shutil
from packaging.specifiers import SpecifierSet
import logging
logger = logging.getLogger(__name__)
import re


import pandas as pd
import numpy as np 
#TODO pull version specs back to minimum working versions
UNICYCLER_VERSION_SPEC = SpecifierSet(">=0.5.1")
PLATON_VERSION_SPEC    = SpecifierSet(">=0.17")
MINIMAP2_VERSION_SPEC  = SpecifierSet(">=2.26")
MINIGRAPH_VERSION_SPEC  = SpecifierSet(">=0.21")




@contextmanager
def temp_fifo():
    """Context Manager for creating named pipes with temporary names."""
    tmpdir = tempfile.mkdtemp()
    filename = os.path.join(tmpdir, 'fifo')  # Temporary filename
    os.mkfifo(filename)  # Create FIFO
    try:
        yield filename
    finally:
        os.unlink(filename)  # Remove file
        os.rmdir(tmpdir)  # Remove directory


def line_count(file):
    with open(file, "rb") as f:
        return  sum(1 for _ in f)

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


def validate_tool(tool_name, vspec, version_cmd="--version", version_split_lambda=lambda x:x.split()[1]):
    tool_path = shutil.which(tool_name)
    if tool_path == None:
        logger.error(f"Cannot find {tool_name}!") 
        return -1
    cmd = [tool_path, version_cmd]
    ret = subprocess.run(cmd, capture_output=True)
    
    tool_version = version_split_lambda(ret.stdout.decode()).strip()

    if tool_version not in vspec:
        logger.warning(f"{tool_name} version:{tool_version} is not supported. {vspec} is required!")
   
    return 0


def run_platon_classifier(args):

    if validate_tool("platon", PLATON_VERSION_SPEC):
        exit(1)
    
    unicycler_fasta_path = f"{args.output_directory}/unicycler_sr/assembly.fasta"
    platon_path = f"{args.output_directory}/classify"

    if not args.force and os.path.isfile(f"{platon_path}/result.tsv"): #TODO implement timestamp checker
        logger.warning(f"Platon file exist at {platon_path} not running it again. Delete the file or use --force")
        return platon_path
    platon_cmd = [
            "platon",
            "-c",
            "--db", args.platon_db,
            "--threads", str(args.threads),
            "--prefix", "result",
            "--output", platon_path,
            unicycler_fasta_path
            ]
    logger.info(f"Running {' '.join(platon_cmd)}")
    ret = subprocess.run(platon_cmd)
    if ret.returncode != 0:
        logger.error(f"Platon failed to finish. Please check its logs at {platon_path}/result.log")
    #    exit(-1)
    return platon_path
def run_unicycler_sr_assembly(args):
    #Check unicycler
    if validate_tool("unicycler", UNICYCLER_VERSION_SPEC):
        exit(1)
    #Check files
    for file in args.short_reads:
        if not os.path.isfile(file):
            logger.error (f"Cannot find {file} !")
            exit(-1)
    unicycler_sr_path = f"{args.output_directory}/unicycler_sr"
    os.makedirs(unicycler_sr_path, exist_ok=True)
    unicycler_cmd = [
            "unicycler",
            "-o", unicycler_sr_path, 
            "-t", str(args.threads),
            "--kmers", "51",
            "-1", args.short_reads[0],
    ]
    if len(args.short_reads) > 1:
        unicycler_cmd.append("-2")
        unicycler_cmd.append(args.short_reads[1])
    
    ret = subprocess.run(unicycler_cmd, capture_output=False) #TODO output capture and tee

    if ret.returncode != 0:
        logger.error(f"Unicycler failed to finish. Please check its logs at {unicycler_sr_path}/unicycler.log")
        exit(-1)
def run_unicycler_lr_assembly(args, plasmid_files_list):
    #Unicycler_runner = "/groups/hachgrp/projects/dev-plasmids/code/Unicycler/unicycler-runner.py"
    Unicycler_runner = "unicycler"
    #Check unicycler
    if validate_tool(Unicycler_runner, UNICYCLER_VERSION_SPEC):
        exit(1)
    #Check files
    cat_reads_cmd = [
        "zcat",
        *plasmid_files_list
    ]


    unicycler_sr_path = f"{args.output_directory}/unicycler_sr"
    unicycler_lr_path = f"{args.output_directory}/unicycler_lr"

    shutil.copytree(unicycler_sr_path, unicycler_lr_path, dirs_exist_ok=True) #TODO only copy required items
    os.unlink(f"{unicycler_lr_path}/assembly.fasta")
    os.unlink(f"{unicycler_lr_path}/assembly.gfa")

    tfw = tempfile.NamedTemporaryFile(delete=False) #TODO switch to tempfile.mkstemp
    ret = subprocess.run(cat_reads_cmd, stdout=tfw, stderr=sys.stderr)
    tfw.close()

    unicycler_cmd = [
            Unicycler_runner,
            "--verbosity", "1",
            "--keep", "3",
            "-o", unicycler_lr_path, 
            "-t", str(args.threads),
            "-l", tfw.name,
    ]

    if len(args.short_reads) == 0:
        #Create empty sr file to fool unicycler

        tfq_fd, tfq_name = tempfile.mkstemp(suffix='.fq')
        tfq_fd2, tfq_name2 = tempfile.mkstemp(suffix='.fq')
        os.fdopen(tfq_fd, 'w').close() 
        os.fdopen(tfq_fd2, 'w').close() 
        unicycler_cmd += ["-1", tfq_name]
        unicycler_cmd += ["-2", tfq_name2]
    else:
        unicycler_cmd += ["-1", args.short_reads[0]]

    print(" ".join(unicycler_cmd), file=sys.stderr)

    if len(args.short_reads) > 1:
        unicycler_cmd.append("-2")
        unicycler_cmd.append(args.short_reads[1])
    
    ret = subprocess.run(unicycler_cmd, capture_output=False) #TODO output capture and tee
    #os.unlink(tfw.name)

    if ret.returncode != 0:
        logger.error(f"Unicycler failed to finish. Please check its logs at {unicycler_sr_path}/unicycler.log")
        exit(-1)
    elif len(args.short_reads) == 0:

        os.unlink(tfq_name)
        os.unlink(tfq_name2)
    return f"{unicycler_lr_path}/assembly.fasta"

def run_minigraph_longreads_to_sr_assembly(args):
    if validate_tool("minigraph", MINIGRAPH_VERSION_SPEC, version_split_lambda=lambda x:x):
        exit(1)
   
    short_read_draft_assembly_graph = f"{args.output_directory}/unicycler_sr/assembly.gfa"
    
    graph_aligment_output_file = f"{args.output_directory}/lr2assembly.gaf"

    if not args.force and os.path.isfile(graph_aligment_output_file):
        logger.warning(f"{graph_aligment_output_file} not running it again. Delete the file or use --force")
        return graph_aligment_output_file
    minigraph_cmd = [
            "minigraph",
            short_read_draft_assembly_graph,
            *args.long_reads,
            "-t", str(args.threads),
            "-x", "lr",
            "-c"
        ]

    with open(graph_aligment_output_file, "w") as hand:
        ret = subprocess.run(minigraph_cmd, stdout=hand, stderr=sys.stderr)
    if ret.returncode != 0:
        logger.error(f"Minigraph failed to finish. Please check its logs at {args.output_directory}/minigraph.log")
        exit(-1)
    return graph_aligment_output_file
def run_long_read_selection(args, prediction_path, graph_alignment_path):
    if validate_tool("split_plasmid_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    
    plasmid_long_reads_path = f"{args.output_directory}/plasmid_long_reads"
    os.makedirs(plasmid_long_reads_path, exist_ok=True)

    plasmid_output = f"{plasmid_long_reads_path}/plasmid.fastq.gz"
    unknown_output = f"{plasmid_long_reads_path}/unknown.fastq.gz"

    if not args.force and os.path.isfile(plasmid_output) and os.path.isfile(unknown_output):
        logger.warning(f"{plasmid_output} and {unknown_output} exists!. not running it again. Delete the files or use --force")
        return plasmid_output, unknown_output

    with tempfile.NamedTemporaryFile(delete=False) as tfw:
        cat_lr_cmd = [
            "zcat",
            *args.long_reads
            ]

        subprocess.run(cat_lr_cmd, stdout=tfw, stderr=sys.stderr)
        tfw.close()

        with open(tfw.name, mode='r') as tfr:
            split_plasmid_read_cmd = [
                    "split_plasmid_reads",
                    graph_alignment_path,
                    tfw.name,
                    prediction_path,
                    plasmid_output,
                    "skip",
                    unknown_output
            ]
            ret = subprocess.run( split_plasmid_read_cmd)

    if ret.returncode != 0:
        logger.error(f"split_plasmid_reads failed to finish. Please check its logs at {args.output_directory}/split_plasmid_reads.log\nWas running {' '.join(split_plasmid_read_cmd)}")
        exit(-1)
    else:
        os.unlink(tfw.name)
    return plasmid_output, unknown_output

def process_platon_output(args, platon_path, max_chr_rds=-7.9, min_plasmid_rds=0.7):
    df = pd.read_csv(f"{platon_path}/result.tsv", sep="\t")

    chr_pos = df.RDS<=max_chr_rds
    pls_pos = df.RDS>=min_plasmid_rds
    rest = ~np.logical_or(chr_pos, pls_pos)
    can_circl = df.Circular=="yes"
    incopat   = df["Inc Type(s)"]>0
    repmob    = (df["# Replication"] + df["# Mobilization"]) > 0
    orit      = df["# OriT"] > 0
    hit       = np.logical_and(np.logical_and(df.RDS > 0.5, df["# Plasmid Hits"] >0), df["# rRNAs"] == 0)
    outpath=f"{platon_path}/result_p.tsv"
    with open(outpath, 'w') as hand:
        print("ID\tPREDICTION", file=hand)
        for i,d in df.loc[pls_pos].iterrows():
            print(f"{d.ID}\tplasmid", file=hand)
        for i,d in df.loc[rest & (can_circl | incopat | repmob | orit | hit)].iterrows():
            print(f"{d.ID}\tplasmid", file=hand)
        for i,d in df.loc[df.RDS<=max_chr_rds].iterrows():
            print(f"{d.ID}\tchromosome", file=hand)

    return outpath

def find_missing_long_reads(args, plasmid_files_list): #TODO Convert processes to unblocking
    
    minimap_output_path = f"{args.output_directory}/prop_lr/lr.round.{len(plasmid_files_list)-1}.paf"
    if not args.force and os.path.isfile(minimap_output_path):
        logger.warning(f"{minimap_output_path} exists!. not running it again. Delete the files or use --force")
        return minimap_output_path
    if validate_tool("innotin", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    if validate_tool("minimap2", MINIMAP2_VERSION_SPEC, version_split_lambda=lambda x:x):
        exit(1)
    tfw = tempfile.NamedTemporaryFile(delete=False)
    cat_reads_cmd = [
        "zcat",
        *plasmid_files_list
    ]


    ret = subprocess.run(cat_reads_cmd, stdout = tfw, stderr=sys.stderr)

    tfw.close()

    innotin_cmd = [
            "innotin",
            *args.long_reads,
            tfw.name
    ]

    tfw2 = tempfile.NamedTemporaryFile(delete=False)
    ret = subprocess.run(innotin_cmd, stdout = tfw2, stderr=sys.stderr)
    os.unlink(tfw.name)

    tfw2.close()
    
    os.makedirs(f"{args.output_directory}/prop_lr", exist_ok=True)



    minimap2_cmd = [
            "minimap2",
            tfw2.name,
            *plasmid_files_list,
            "-o", minimap_output_path,
            "-t", str(args.threads)
    ]
    
    ret = subprocess.run(minimap2_cmd)
    if ret.returncode != 0:
        logger.error(f"Minimap failed to finish. Please check its logs at {args.output_directory}/prop_lr/minimap2.log")
        exit(-1)
    return minimap_output_path
def extract_missing_long_reads(args, plasmid_alignment, graph_alignment_path, prediction_tsv_path):

    extracted_lr_fastq = '.fastq.gz'.join(plasmid_alignment.rsplit('.paf', 1))
    if not args.force and os.path.isfile(extracted_lr_fastq):
        logger.warning(f"{extracted_lr_fastq} exists!. not running it again. Delete the files or use --force")
        return extracted_lr_fastq
    if validate_tool("select_missing_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    tfw =  tempfile.NamedTemporaryFile(delete=False)
    cat_lr_cmd = [
        "zcat",
        *args.long_reads
        ]

    subprocess.run(cat_lr_cmd, stdout=tfw, stderr=sys.stderr)
    tfw.close()


    smr_cmd = [
            "select_missing_reads",
            plasmid_alignment,
            graph_alignment_path,
            tfw.name,
            extracted_lr_fastq,
            prediction_tsv_path
    ]

    ret = subprocess.run(smr_cmd)
    os.unlink(tfw.name)

    if ret.returncode != 0:
        logger.error(f"Long read extraction failed!")
        exit(ret.returncode)
    return extracted_lr_fastq

def save_forgotten_short_read_only_circular_contigs(args):

    #Step 1 Map Circular contigs in the short-read assembly to circular contigs in the long-read assembly.
    unicycler_sr_path = f"{args.output_directory}/unicycler_sr/assembly.fasta"
    

    sr_tfw = tempfile.NamedTemporaryFile(delete=False)
    for name, header, seq in generate_fasta(unicycler_sr_path):
        if "circular" in header:
            sr_tfw.write(f">{name} {header}\n{seq}\n".encode())
    sr_tfw.close()

    unicycler_lr_path = f"{args.output_directory}/unicycler_lr/assembly.fasta"
    
    lr_tfw = tempfile.NamedTemporaryFile(delete=False)
    for name, header, seq in generate_fasta(unicycler_lr_path):
        if "circular" in header:
            lr_tfw.write(f">{name} {header}\n{seq}\n".encode())
    lr_tfw.close()
    minimap_output_path = f"{args.output_directory}/save_missed_contigs/src2lrc.paf" 

    os.makedirs(f"{args.output_directory}/save_missed_contigs", exist_ok=True)
    minimap2_cmd = [
            "minimap2",
            lr_tfw.name,
            sr_tfw.name,
            "-o", minimap_output_path,
            "-t", str(args.threads)
    ]
    ret = subprocess.run(minimap2_cmd)
    #Step 2 Find Skipped Contigs in the mapping file
    
    with open(minimap_output_path, 'r') as hand:
        for line in hand:
            line = line.rstrip().split("\t")
            print(line)
    os.unlink(sr_tfw.name)
    os.unlink(lr_tfw.name)
def save_forgotten_miniasm_only_circular_contigs(args, min_chromosomal_length=1_000_000, min_cov_ratio=.8, min_allowed_sr_coverage=.8):
    #Step 1 Check if a circular miniasm unitig exists

    import mappy as mp
    miniasm_gfa_path = f"{args.output_directory}/unicycler_lr/miniasm_assembly/12_unitig_graph.gfa"
    
    if not os.path.isfile(miniasm_gfa_path):
        return []
    names = []
    maseqs = []
    links = []
    for line in open(miniasm_gfa_path, 'r'):
        if line[0] == "L":
            line=line.rstrip().split("\t")
            if line[1] == line[3]:
                links.append(line[1])

    for line in open(miniasm_gfa_path, 'r'):
        if line[0] == "S":
            line = line.split("\t")
            if True:#line[1] in links:
                maseqs.append(line[2])
    #Step 2 Check if there is a incomplete of the same sequence in the unicycler output 
    unicycler_lr_path = f"{args.output_directory}/unicycler_lr/assembly.fasta"
    to_recover = []
    uni_seqs = []
    for n, h, s in generate_fasta(unicycler_lr_path):
        if "circular" not in h:
            continue
        if len(s) > min_chromosomal_length:
            continue

        uni_seqs.append(s)

    for seq in maseqs:
        mapper = mp.Aligner(seq=seq+seq)
        for s in uni_seqs:
            maximal_mapping = -1
            for alig in mapper.map(s):
                if alig.mlen > maximal_mapping:
                    maximal_mapping = alig.mlen
            if maximal_mapping > min_cov_ratio * len(seq):
                break
        else:
            #Uncovered

            to_recover.append(seq)

    #Step 3 Using the mappings break apart the unitig
    if len(to_recover) == 0:
        return []
    unicycler_sr_gfa_path = f"{args.output_directory}/unicycler_sr/assembly.gfa"

    tfw, tfw_name = tempfile.mkstemp()

    with os.fdopen(tfw, 'w') as hand:
        for i, seq in enumerate(to_recover):
            print(f">{i}\n{seq}",file=hand)

    minigraph_cmd = [
        "minigraph",
        unicycler_sr_gfa_path,
        tfw_name,
        "-t", str(args.threads),
        "-x", "asm"
    ]
    ret = subprocess.run(minigraph_cmd, capture_output=True)
    if ret.returncode != 0:
        logger.error(f"minigraph failed to finish. Please check its logs at {platon_path}/result.log")
        exit(1)
    print(" ".join(minigraph_cmd))

    sr_contigs = {}
    with open(unicycler_sr_gfa_path) as hand: 
        for line in hand:
            if line[0] == "S":
                line = line.rstrip().split("\t")
                sr_contigs[line[1]] = line[2]

#0       120417  5       120401  +       >33>95>484>74<230<119<71<283>63<368<48<484<95<53>283<661>381>523>321>668<216<248>228<664>291<694<383<206<164>685>395<123>632>119>173<179>57>93>117<93<69>217>67<417>62>217>553<417>33   171147  39974       160988  115301  121107  60      tp:A:P  cm:i:17636      s1:i:113914     s2:i:0  dv:f:0.0114
    maximal_mappings = np.zeros(len(to_recover))
    aligs = {}
    output =  ret.stdout.decode().rstrip()
    if output:
        print(output)
        for line in output.split("\n"):
            line = line.rstrip().split("\t")
            cid = int(line[0])
            _s = int(line[2])
            _e = int(line[3])
            _l = int(line[1])
            if (_e - _s) / _l > maximal_mappings[cid]:
                maximal_mappings[cid] = (_e - _s) / _l
                aligs[cid] = line
    else:
        os.unlink(tfw_name)
        return []
    #Step 4 Replace parts of the unitig with unicycler sr contigs

    recovered = []
    for i,m in enumerate(maximal_mappings):
        if m > min_allowed_sr_coverage:
            new_seq = ""
            contigs = re.findall("[<,>][0-9]+", aligs[i][5])
            if len(contigs) < 2 or contigs[0][1:] != contigs[-1][1:]:
                continue
            for octg in  contigs:
                orient = octg[0]
                ctg = octg[1:]
                if orient == ">":
                    new_seq = new_seq + sr_contigs[ctg]
                elif orient == "<":
                    new_seq = new_seq + mp.revcomp(sr_contigs[ctg])
            _s = int(aligs[i][7])
            _e = int(aligs[i][8])
            recovered.append(new_seq[_s:_e])
    
    os.unlink(tfw_name)

    return recovered
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--long-reads",help="long reads fastq file", nargs="+", default=[])
    parser.add_argument("-s","--short-reads",help="short reads fastq files",nargs="+",default=[])
    parser.add_argument("--sr-assembly",help="short reads assembly graph")
    parser.add_argument("-o","--output-directory",required=True)
    parser.add_argument("-p", "--propagate-rounds", help="Number of rounds to propagate plasmid long read information", type=int, default=0)
    parser.add_argument("--platon-db",required=True)
    parser.add_argument("-t","--threads",type=int, default=16)
    parser.add_argument("--verbosity", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")
    parser.add_argument("--force", action='store_true')
    args = parser.parse_args()
    logging.basicConfig( level=args.verbosity)
    min_fasta_seq_len = 200



    if len(args.short_reads) == 0:
        print("Short reads are not provided!", file=sys.stderr)
    if len(args.long_reads) == 0:
        print("Long reads are not provided!", file=sys.stderr)

    logger.info(args)
    
    if args.sr_assembly is not None:
        _gfa_path =f"{args.output_directory}/unicycler_sr/002_depth_filter.gfa"
        _gfa_path2 =f"{args.output_directory}/unicycler_sr/assembly.gfa"
        _fasta_path = f"{args.output_directory}/unicycler_sr/assembly.fasta"
        unicycler_sr_path = f"{args.output_directory}/unicycler_sr"
        os.makedirs(unicycler_sr_path, exist_ok=True)
        shutil.copyfile(args.sr_assembly, _gfa_path)
        shutil.copyfile(args.sr_assembly, _gfa_path2)
        with open (_gfa_path,'r') as hand, open(_fasta_path, 'w') as whand:
            for line in hand:
                if line[0] == "S":
                    line=line.rstrip().split("\t")
                    if len(line[2]) > min_fasta_seq_len:
                        print(f">{line[1]}\n{line[2]}", file=whand)

    else:
        run_unicycler_sr_assembly(args)

    platon_path = run_platon_classifier(args)
    
    prediction_tsv_path = process_platon_output(args, platon_path)
    
    graph_alignment_path = run_minigraph_longreads_to_sr_assembly(args)
    
    #If there are no alignments
    # Check the platon output to grab circularized sr plasmids
    # and terminate

    plasmid_reads_file, unknown_reads_file = run_long_read_selection(args, prediction_tsv_path, graph_alignment_path)
    
    plasmid_files = [unknown_reads_file, plasmid_reads_file]
    for i in range(args.propagate_rounds):
        plasmid_alignment = find_missing_long_reads(args, plasmid_files)
        if line_count(plasmid_alignment) == 0:
            logger.info(f"{plasmid_alignment} has no alignments. Stopping the propagation!")
            break
        plasmid_reads_i = extract_missing_long_reads(args, plasmid_alignment, graph_alignment_path, prediction_tsv_path)
        if line_count(plasmid_reads_i) == 0:
            logger.info(f"{plasmid_reads_i} has no reads. Stopping the propagation!")
            break
        plasmid_files.append(plasmid_reads_i)
    
    lr_assembly_path = run_unicycler_lr_assembly(args, plasmid_files)

    final_assembly_path = f"{args.output_directory}/assembly.final.fasta"

    with open(final_assembly_path, 'w') as hand:
        for n, h, s in generate_fasta(lr_assembly_path):
            if "circular" in h:
                print(f">{n} {h}\n{s}", file=hand)
        

        #for i,seq in enumerate(save_forgotten_miniasm_only_circular_contigs(args)):
        #    print(f">l_{i} circular=True\n{seq}", file=hand)
    
    #save_forgotten_short_read_only_circular_contigs(args)

    return 0

if __name__ == "__main__":
    exit(main())
