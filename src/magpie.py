import sys
import argparse
import os
import shutil
from pathlib import Path
import subprocess
import gfapy
import pandas as pd
import sourmash
import gzip
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
    
def hash_ref(ref_path, n, ksize, scaled):
    ref_hashes = {}
    ref_lens   = {}
    for name,_,seq in generate_fasta(ref_path):
        ref_hashes[name] = sourmash.MinHash(n=n,ksize=ksize,scaled=scaled)
        ref_hashes[name].add_sequence(seq,force=True)
        ref_lens[name] = len(seq)
    return ref_hashes, ref_lens
"""
usage: raven [options ...] <sequences> [<sequences> ...]

  # default output is to stdout in FASTA format
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -k, --kmer-len <int>
      default: 15
      length of minimizers used to find overlaps
    -w, --window-len <int>
      default: 5
      length of sliding window from which minimizers are sampled
    -f, --frequency <double>
      default: 0.001
      threshold for ignoring most frequent minimizers
    -i, --identity <double>
      default: 0
      threshold for overlap between two reads in order to construct an edge between them
    -o, --kMaxNumOverlaps <long unsigned int>
      default: 32
      maximum number of overlaps that will be taken during FindOverlapsAndCreatePiles stage
    -M, --use-micromizers
      if this is enabled micromizers will be used instead of mimizers in graph construction
      (performance will increase slightly and memory consumption will decrease but results could be slightly worse)
    -p, --polishing-rounds <int>
      default: 2
      number of times racon is invoked
    -m, --match <int>
      default: 3
      score for matching bases
    -n, --mismatch <int>
      default: -5
      score for mismatching bases
    -g, --gap <int>
      default: -4
      gap penalty (must be negative)
    -u, --min-unitig-size <int>
      minimal uniting size (default 9999)
    -F --graphical-fragment-assembly <string>
      prints the assembly graph in GFA format
    -U --unitig-graphical-fragment-assembly <string>
      prints the unitig assembly graph in GFA format
    --resume
      resume previous run from last checkpoint
    --disable-checkpoints
      disable checkpoint file creation
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage
"""


def raven_main(args):

    with open (args.output, 'w') as hand:
        cmd = [
            "raven",
            "--threads",
            str(args.threads),
            "-F",
            args.graphical_fragment_assembly
        ] + args.reads
        #print(cmd)
        return subprocess.run(cmd,stdout=hand)

def flye_main(args):
    outpath = Path(args.output)
    basename, filename = os.path.split(outpath)
    filebase = filename.split(".")[0]
    if basename:
        tempdir = f"{basename}/{filebase}_dir"
    else:
        tempdir = f"{filebase}_dir"
    os.makedirs(tempdir, exist_ok=True)


    cmd = [
        "flye",
        f"--{args.read_type}",
        "--out-dir",
        tempdir,
        "--threads",
        str(args.threads),
    ] + args.reads
    result = subprocess.run(cmd)
    if result.return_code != 0:
        #print(f"Flye on {args.reads} failed Program Output:\n{result.stdout}\nProgram Error:{result.stderr}", file=sys.stderr)
        print(f"Flye on {args.reads} failed", file=sys.stderr)
        return result.return_code

    shutil.copyfile(f"{tempdir}/assembly.fasta", args.output)
    shutil.copyfile(f"{tempdir}/assembly_graph.gfa", args.graphical_fragment_assembly)

def main():
    parser = argparse.ArgumentParser()
   
    parser.add_argument("-t","--threads", type=int)
    sbp = parser.add_subparsers(required=True)
    
    raven_parser = sbp.add_parser('raven')
    raven_parser.add_argument("reads", nargs="*")

    raven_parser.add_argument("-F", "--graphical-fragment-assembly", help="prints the assembly graph in GFA format")
    raven_parser.add_argument("-O", "--output", help="Output Fasta")
    raven_parser.set_defaults(prog=raven_main)
    
    flye_parser = sbp.add_parser('flye')
    flye_parser.set_defaults(prog=flye_main)
    flye_parser.add_argument("reads")
    flye_parser.add_argument("--read-type", default="nano-raw")
    flye_parser.add_argument("-F", "--graphical-fragment-assembly", help="prints the assembly graph in GFA format")
    flye_parser.add_argument("-O", "--output", help="Output Fasta")
    args = parser.parse_args()
    if not isinstance(args.reads, list):
        args.reads = [args.reads]
    args.prog(args);


    return 0

if __name__ == "__main__":
    exit(main())
