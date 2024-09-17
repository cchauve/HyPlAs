import argparse
import os
import sys
import subprocess
import tempfile
from shutil import which
from packaging.specifiers import SpecifierSet
import logging
logger = logging.getLogger(__name__)

import pandas as pd

#TODO pull version specs back to minimum working versions
UNICYCLER_VERSION_SPEC = SpecifierSet(">=0.5.1")
PLATON_VERSION_SPEC    = SpecifierSet(">=0.17")
MINIMAP2_VERSION_SPEC  = SpecifierSet(">=2.26")
MINIGRAPH_VERSION_SPEC  = SpecifierSet(">=0.21")

def line_count(file):
    with open(file, "rb") as f:
        return  sum(1 for _ in f)


def validate_tool(tool_name, vspec, version_cmd="--version", version_split_lambda=lambda x:x.split()[1]):
    tool_path = which(tool_name)
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

    if not args.force and os.path.isfile(f"{platon_path}/result.tsv"):
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
            "-1", args.short_reads[0],
    ]
    if len(args.short_reads) > 1:
        unicycler_cmd.append("-2")
        unicycler_cmd.append(args.short_reads[1])
    
    ret = subprocess.run(unicycler_cmd, capture_output=False) #TODO output capture and tee

    if ret.returncode != 0:
        logger.error(f"Unicycler failed to finish. Please check its logs at {unicycler_sr_path}/unicycler.log")
        exit(-1)

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
    if validate_tool("./src/split_plasmid_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
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
        print(tfw.name)
        with open(tfw.name, mode='r') as tfr:
            split_plasmid_read_cmd = [
                    "./src/split_plasmid_reads",
                    graph_alignment_path,
                    tfw.name,
                    prediction_path,
                    plasmid_output,
                    "skip",
                    unknown_output
            ]
            ret = subprocess.run( split_plasmid_read_cmd)
        os.unlink(tfw.name)
    if ret.returncode != 0:
        logger.error(f"split_plasmid_reads failed to finish. Please check its logs at {args.output_directory}/split_plasmid_reads.log")
        exit(-1)
    return plasmid_output, unknown_output

def process_platon_output(args, platon_path, max_chr_rds=-100000, min_plasmid_rds=10):
    df = pd.read_csv(f"{platon_path}/result.tsv", sep="\t")
    outpath=f"{platon_path}/result_p.tsv"
    with open(outpath, 'w') as hand:
        print("ID\tPREDICTION", file=hand)
        for i,d in df.iterrows():
            if d.RDS <= max_chr_rds:
                print(f"{d.ID}\tchromosome", file=hand)
            elif d.RDS >= min_plasmid_rds:
                print(f"{d.ID}\tplasmid", file=hand)

    return outpath

def find_missing_long_reads(args, plasmid_files_list): #TODO Convert processes to unblocking
    if validate_tool("./src/innotin", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
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
            "./src/innotin",
            *args.long_reads,
            tfw.name
    ]

    tfw2 = tempfile.NamedTemporaryFile(delete=False)
    ret = subprocess.run(innotin_cmd, stdout = tfw2, stderr=sys.stderr)
    os.unlink(tfw.name)

    tfw2.close()
    
    os.makedirs(f"{args.output_directory}/prop_lr", exist_ok=True)

    minimap_output_path = f"{args.output_directory}/prop_lr/lr.round.{len(plasmid_files_list)-1}.paf"

    minimap2_cmd = [
            "minimap2",
            tfw2.name,
            plasmid_files_list[-1],
            "-o", minimap_output_path,
            "-t", str(args.threads)
    ]
    
    ret = subprocess.run(minimap2_cmd)
    if ret.returncode != 0:
        logger.error(f"Minimap failed to finish. Please check its logs at {args.output_directory}/prop_lr/minimap2.log")
        exit(-1)
    return minimap_output_path
def extract_missing_long_reads(args, plasmid_alignment, graph_alignment_path, prediction_tsv_path):
    if validate_tool("./src/select_missing_reads", SpecifierSet(">0"), version_cmd="", version_split_lambda=lambda x:"1"):
        exit(1)
    tfw =  tempfile.NamedTemporaryFile(delete=False)
    cat_lr_cmd = [
        "zcat",
        *args.long_reads
        ]

    subprocess.run(cat_lr_cmd, stdout=tfw, stderr=sys.stderr)
    tfw.close()

    extracted_lr_fastq = '.fastq.gz'.join(plasmid_alignment.rsplit('.paf', 1))
    smr_cmd = [
            "./src/select_missing_reads",
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
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--long-reads",help="long reads fastq file", nargs="+", default=[])
    parser.add_argument("-s","--short-reads",help="short reads fastq files",nargs="+",default=[])
    parser.add_argument("-o","--output-directory",required=True)
    parser.add_argument("-p", "--propagate-rounds", help="Number of rounds to propagate plasmid long read information", type=int, default=0)
    parser.add_argument("--platon-db",required=True)
    parser.add_argument("-t","--threads",type=int, default=16)
    parser.add_argument("--verbosity", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")
    parser.add_argument("--force", action='store_true')
    args = parser.parse_args()
    logging.basicConfig( level=args.verbosity)
    
    if len(args.short_reads) == 0:
        print("Short reads are not provided!", file=sys.stderr)
    if len(args.long_reads) == 0:
        print("Long reads are not provided!", file=sys.stderr)

    logger.info(args)
    
#    run_unicycler_sr_assembly(args)
    
    platon_path = run_platon_classifier(args)
    
    prediction_tsv_path = process_platon_output(args, platon_path)
    
    graph_alignment_path = run_minigraph_longreads_to_sr_assembly(args)

    plasmid_reads_file, unknown_reads_file = run_long_read_selection(args, prediction_tsv_path, graph_alignment_path)
    
    plasmid_files = [plasmid_reads_file]
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

    return 0

if __name__ == "__main__":
    exit(main())
