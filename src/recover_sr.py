import sys
import argparse

import subprocess
import gfapy
import pandas as pd
import sourmash
import gzip
from collections import defaultdict
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

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-t","--threads", type=int)
    parser.add_argument("-G", "--draft-gfa", help="Draft GFA assembly")
    parser.add_argument("-L", "--lr-assembly", help="Input Fasta")
    parser.add_argument("-SO", "--sr-output", help="Output Fasta")
    parser.add_argument("-GO", "--output-gfa", help="Output Gfa")
    parser.add_argument("-T", "--plasmid-prediction", help="Platon output")
    parser.add_argument("--only_pick_plasmid_predicted_circular", action='store_true')
    parser.add_argument("--similarity-threshold", type=float, default=0.5)
    parser.add_argument("--contig-length-threshold", type=int, default=1000)
    args = parser.parse_args()
    plasmid_predictions = pd.read_csv(args.plasmid_prediction, sep="\t",dtype=defaultdict(lambda :"str"))
    circular_plasmids = {v.ID for k,v in plasmid_predictions.iterrows()}



    with open(args.sr_output, 'w') as hand:
        lr_ctg_hashes = hash_ref(args.lr_assembly, n=100, ksize=31, scaled=0)

        gfa_sr = gfapy.gfa.Gfa.from_file(args.draft_gfa,vlevel=0)
        gfa_out = gfapy.gfa.Gfa()
        for e in gfa_sr.edges:
            if e.from_name == e.to_name:
                print(f"Testing {e.from_name}",file=sys.stderr)
                node = gfa_sr.segment(e.from_name)
                if  ( args.only_pick_plasmid_predicted_circular and e.from_name not in circular_plasmids):
                    print(f"ID {e.from_name} is not in {circular_plasmids}",file =sys.stderr)
                    continue
                if len(node.sequence) < args.contig_length_threshold:
                    print(f"Too short: {len(node.sequence)}", file=sys.stderr)
                    continue
                sh = sourmash.MinHash(n=100,ksize=31,scaled=0)
                sh.add_sequence(node.sequence)
                for k,v in lr_ctg_hashes[0].items():
                    sim = sh.jaccard(v) >= args.similarity_threshold
                    if sim:
                        print(f"Skipping {e.from_name} because it has similarity of {sim} with lr-contig {k}", file=sys.stderr)
                        break
                else:
                    gfa_out.add_line(str(node))
                    gfa_out.add_line(str(e))
                    print(f"{node.name} is accepted", file=sys.stderr)
                    print(f">{node.name}\n{node.sequence}", file=hand)

    gfa_out.to_file(args.output_gfa)
    return 0

if __name__ == "__main__":
    exit(main())
