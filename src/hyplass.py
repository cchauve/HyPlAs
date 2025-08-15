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
from collections import defaultdict
import pathlib
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
    """
        equivalent to `wc -l`
    """
    with open(file, "rb") as f:
        return  sum(1 for _ in f)

def generate_fasta(fasta):
    """
        Given path to a fasta file, return a generator (name, header, sequence) of each record in the fasta file.
        :param fasta: path to a fasta
        :returns: generator of tuple[str, str, str]
        @example:
            for name, header, seq in generate_fasta(fasta):
                # do work 
    """
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
    """
        Validate the tool and its version
        :param tool_name: the name of tool
        :param vspec: Version specification
        :param version_cmd: the command to get the tool's version
        :param version_split_lambda: the lambda function to split the version string
        :return: 0 for success, -1 otherwise
        @example:
            validate_tool("platon", SpecifierSet(">=0.17"))
    """
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
    """
        Run Platon classifier on unicycler short-read-only assembly.
        :param args: argparse.Namespace object containing:
            platon_db: path to platon database
            threads: number of threads allowed to used
            output_directory: path to output directory
            force: force to run Platon even if output exists
        :returns: path to Platon output directory
    """
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
    """
        Run unicycler assembly on short-reads
        :param args: argparse.Namespace object containing: 
            short_reads: list of paths to fastq files
            output_directory: path to directory where output files will be stored
            threads: number of threads to use
        :return: path to directory containing assembly results
    """
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
            #"--kmers", "51",
            "-1", args.short_reads[0],
            "--min_component_size", "10"
    ]
    if len(args.short_reads) > 1:
        unicycler_cmd.append("-2")
        unicycler_cmd.append(args.short_reads[1])
    
    ret = subprocess.run(unicycler_cmd, capture_output=False) #TODO output capture and tee

    if ret.returncode != 0:
        logger.error(f"Unicycler failed to finish. Please check its logs at {unicycler_sr_path}/unicycler.log")
        exit(-1)
    return unicycler_sr_path

def run_unicycler_lr_assembly(args, plasmid_files_list, it):
    """
        Run Unicycler for LR assembly
        :param args: argparse.Namespace object containing: 
            output_directory: path to directory where output files will be stored
            threads: number of threads to use
            short_reads: list of paths to short reads
        :param plasmid_files_list: list of paths to plasmid files
        :param it: iteration number
        :return: path to assembly fasta file
    """
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
    unicycler_lr_path = f"{args.output_directory}/unicycler_lr_{it}"

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
    """
        Run minigraph to assemble long reads into a single assembly
        :param args: argparse.Namespace object containing:
                output_directory: path to the output directory
                short_reads: path to the short reads
                threads: number of threads to use
                force: whether to force the execution of the command
        :return: str: path to the graph alignment file
    """
    if validate_tool("minigraph", MINIGRAPH_VERSION_SPEC, version_split_lambda=lambda x:x):
        exit(1)
   

    short_read_draft_assembly_graph = f"{args.output_directory}/unicycler_sr/assembly.gfa"
    short_read_draft_assembly_graph_fix = f"{args.output_directory}/unicycler_sr/assembly_segfix.gfa"
    fix_gfa_empty_segments(short_read_draft_assembly_graph, short_read_draft_assembly_graph_fix) 
    graph_aligment_output_file = f"{args.output_directory}/lr2assembly.gaf"

    if not args.force and os.path.isfile(graph_aligment_output_file):
        logger.warning(f"{graph_aligment_output_file} not running it again. Delete the file or use --force")
        return graph_aligment_output_file
    minigraph_cmd = [
            "minigraph",
            short_read_draft_assembly_graph_fix,
            args.long_reads,
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
    """
        Run long read selection
        :param args: argparse.Namespace object containing:
                output_directory: path to the output directory
                force: rerun the tool even if results are already present
                long_reads: path to the long reads
        :param prediction_path: path to the plasmid classification results
        :param graph_alignment_path: path to the graph alignment file
        :return: tuple[str, str, str, str] paths to plasmid, unknown(neither plasmid or chr), unknown(both plasmid and chr) and unmapped reads.
    """
    if validate_tool("split_plasmid_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    
    plasmid_long_reads_path = f"{args.output_directory}/plasmid_long_reads"
    os.makedirs(plasmid_long_reads_path, exist_ok=True)

    plasmid_output = f"{plasmid_long_reads_path}/plasmid.fastq.gz"
    unknown_output_both = f"{plasmid_long_reads_path}/unknown_both.fastq.gz"
    unknown_output_neit = f"{plasmid_long_reads_path}/unknown_neither.fastq.gz"
    unmapped_output = f"{plasmid_long_reads_path}/unmapped.fastq.gz"

    if not args.force and os.path.isfile(plasmid_output) and os.path.isfile(unknown_output_both):
        logger.warning(f"{plasmid_output} and {unknown_output_both} exists!. not running it again. Delete the files or use --force")
        return plasmid_output, unknown_output_both, unknown_output_neit, unmapped_output




    split_plasmid_read_cmd = [
            "split_plasmid_reads",
            graph_alignment_path,
            args.long_reads,
            prediction_path,
            plasmid_output,
            "skip",
            unknown_output_neit,
            unknown_output_both,
            unmapped_output,
    ]
    print(" ".join(split_plasmid_read_cmd), file=sys.stderr)
    ret = subprocess.run( split_plasmid_read_cmd)

    if ret.returncode != 0:
        logger.error(f"split_plasmid_reads failed to finish. Please check its logs at {args.output_directory}/split_plasmid_reads.log\nWas running {' '.join(split_plasmid_read_cmd)}")
        exit(-1)

    return plasmid_output, unknown_output_neit, unknown_output_both, unmapped_output

def process_platon_output(args, platon_path, max_chr_rds=-7.9, min_plasmid_rds=0.7):
    """
        Process platon output and write plasmid and chromosome predicted contigs to a file
        :param args: argparse.Namespace object (Not used)
        :param platon_path: Path to the platon output (classification) directory
        :param max_chr_rds: Maximum RDS value for a chromosome read to be considered as a chromosome. Default is -7.9.
        :param min_plasmid_rds: Minimum RDS value for a plasmid read to be considered as a plasmid. Default is 0.7.
        :return: Path to the predictions
    """
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

def find_missing_long_reads(args, plasmid_files_list, unknown_files): #TODO Convert processes to unblocking
    """
        Maps the unknown reads to the plasmid reads and recovers missed long reads
        :param args: argparse.Namespace object containing:
                output_directory: str path to the directory
                force: run the script even if the output file exists
                threads: number of threads to use with minimap2
        :param plasmid_files_list: list of paths to the plasmid files. Starts with single one and each iteration of
                find_missing_long_reads() and extract_missing_long_reads() adds one more plasmid file
        :param unknown_files: lists of unkown reads including unmapped, both plasmid+chr, neither plasmid and chr
        :return: Path to minimap2 output paf file path
    """
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
        *unknown_files
    ]


    ret = subprocess.run(cat_reads_cmd, stdout = tfw, stderr=sys.stderr)

    tfw.close()

    innotin_cmd = [
            "innotin",
            tfw.name,
            *plasmid_files_list
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
            "-t", str(args.threads),
    ]
    
    ret = subprocess.run(minimap2_cmd)
    if ret.returncode != 0:
        logger.error(f"Minimap failed to finish. Please check its logs at {args.output_directory}/prop_lr/minimap2.log")
        exit(-1)
    return minimap_output_path

def extract_missing_long_reads(args, plasmid_alignment, graph_alignment_path, prediction_tsv_path, unknown_reads): #TODO optimize to only use unknown reads.
    """
        Extracts the missing long reads from the minimap2 alignment.
        :param args: argparse.Namespace object containing:
            output_directory: The path to the output directory.
            long_reads: path to long reads.
            force: If True, will overwrite the output file if it exists.
        :param graph_alignment_path: path to the long read to short-read-assembly graph alignment.
        :param plasmid_alignment: path to the minimap2 alignment of the plasmid reads to unknown reads
        :param prediction_tsv_path: path to the platon prediction file.
        :return Path to the extracted long reads fastq file
    """
    extracted_lr_fastq = '.fastq.gz'.join(plasmid_alignment.rsplit('.paf', 1))
    if not args.force and os.path.isfile(extracted_lr_fastq):
        logger.warning(f"{extracted_lr_fastq} exists!. not running it again. Delete the files or use --force")
        return extracted_lr_fastq
    if validate_tool("select_missing_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    

    cat_reads_cmd = [
        "zcat",
        *unknown_reads
    ]

    tfw = tempfile.NamedTemporaryFile(delete=False)
    ret = subprocess.run(cat_reads_cmd, stdout = tfw, stderr=sys.stderr)
    tfw.close()
    smr_cmd = [
            "select_missing_reads_nofilter",
            plasmid_alignment,
            graph_alignment_path,
            tfw.name,
            extracted_lr_fastq,
            prediction_tsv_path
    ]

    ret = subprocess.run(smr_cmd)


    if ret.returncode != 0:
        logger.error(f"Long read extraction failed!")
        exit(ret.returncode)
    return extracted_lr_fastq


def fix_gfa_empty_segments(gfa_path, out_path):
    """
        Fixes gfa files with empty segments
        :param gfa_path: gfa file path
        :param out_path: output path
    """
    empty_segments = set()
    incoming = defaultdict(list)
    outgoing = defaultdict(list)
    with open(out_path, 'w') as hand:
        for line in open(gfa_path, 'r'):
            fields = line.split("\t")
            if line[0] == "S":
                if fields[2] == "":
                    empty_segments.add(fields[1])
                else:
                    print(line.rstrip(), file=hand)
            elif line[0] == "L":
                if fields[1] in empty_segments:
                    outgoing[fields[1]].append((fields[3],fields[4]))
                elif fields[3] in empty_segments:
                    incoming[line[3]].append((fields[1],fields[2]))
                else:
                    print(line.rstrip(), file=hand)
        for seg in empty_segments:
            for i in incoming[seg]:
                for o in outgoing[seg]:
                    print(f"L\t{i[0]}\t{i[1]}\t{o[0]}\t{o[1]}\t0M", file=hand)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--long-reads",help="long reads fastq file")
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
    if not args.long_reads:
        print("Long reads are not provided!", file=sys.stderr)

    logger.info(args)


    if args.sr_assembly is not None:
        _gfa_path =f"{args.output_directory}/unicycler_sr/002_depth_filter.gfa"
#        _gfa_path2 =f"{args.output_directory}/unicycler_sr/assembly.gfa"
        _fasta_path = f"{args.output_directory}/unicycler_sr/assembly.fasta"
        unicycler_sr_path = f"{args.output_directory}/unicycler_sr"
        os.makedirs(unicycler_sr_path, exist_ok=True)
        shutil.copyfile(args.sr_assembly, _gfa_path)
        #fix_gfa_empty_segments(unicycler_sr_path, _gfa_path)
#        shutil.copyfile(args.sr_assembly, _gfa_path2)
        subprocess.run(["unicycler", "-s", "src/mock.fq", "-o", unicycler_sr_path]) 
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
    final_assembly_path = f"{args.output_directory}/assembly.final.it{{}}.fasta"

    plasmid_reads_file, unknown_reads_file_both, unknown_reads_file_neither, unmapped_reads_file = run_long_read_selection(args, prediction_tsv_path, graph_alignment_path)
    with open(plasmid_reads_file, "rb") as f:
        num_lines = sum(1 for _ in f)  
    if num_lines < 4:
        print("No plasmid long reads are found",file=sys.stderr)
        unicycler_fasta_path = f"{args.output_directory}/unicycler_sr/assembly.fasta"
        with open(final_assembly_path.format(0), 'w') as hand:
            for n, h, s in generate_fasta(unicycler_fasta_path):
                if "circular" in h:
                    print(f">{n} {h}\n{s}", file=hand)
        for i in range(0,args.propagate_rounds):
            pathlib.Path(final_assembly_path.format(i+1)).symlink_to(final_assembly_path.format(0))
        return 0
    plasmid_files = [plasmid_reads_file]
    lr_assembly_path = run_unicycler_lr_assembly(args, plasmid_files, 0)
    with open(final_assembly_path.format(0), 'w') as hand:
        for n, h, s in generate_fasta(lr_assembly_path):
            if "circular" in h:
                print(f">{n} {h}\n{s}", file=hand)
    
    for i in range(args.propagate_rounds):
        plasmid_alignment = find_missing_long_reads(args, plasmid_files, [unmapped_reads_file, unknown_reads_file_both, unknown_reads_file_neither])
        #plasmid_alignment = find_missing_long_reads(args, args.long_reads, plasmid_files)
        if line_count(plasmid_alignment) == 0:
            for it in range(i, args.propagate_rounds):
                pathlib.Path(final_assembly_path.format(it+1)).symlink_to(final_assembly_path.format(i))
            logger.info(f"{plasmid_alignment} has no alignments. Stopping the propagation!")
            break
        plasmid_reads_i = extract_missing_long_reads(args, plasmid_alignment, graph_alignment_path, prediction_tsv_path, [unknown_reads_file_both, unknown_reads_file_neither, unmapped_reads_file])
        if line_count(plasmid_reads_i) == 0:
            for it in range(i, args.propagate_rounds):
                pathlib.Path(final_assembly_path.format(it+1)).symlink_to(final_assembly_path.format(i))
            logger.info(f"{plasmid_reads_i} has no reads. Stopping the propagation!")
            break
        plasmid_files.append(plasmid_reads_i)
    
        lr_assembly_path = run_unicycler_lr_assembly(args, plasmid_files, i+1)
        with open(final_assembly_path.format(i+1), 'w') as hand:
            for n, h, s in generate_fasta(lr_assembly_path):
                if "circular" in h:
                    print(f">{n} {h}\n{s}", file=hand)
    
    return 0

if __name__ == "__main__":
    exit(main())
