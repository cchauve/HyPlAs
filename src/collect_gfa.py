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

    lr_gfa = sys.argv[1]
    lr_pil = sys.argv[2]
    sr_gfa = sys.argv[3]
    out_gfa = sys.argv[4]

    
    gfa_out = gfapy.Gfa.from_file(lr_gfa)
    
    for n,h,s in generate_fasta(lr_pil):
        for s in gfa_out.segments:
            if s.name == n:
                s.sequence = s
                break

    for line in open(sr_gfa,'r'):
        gfa_out.add_line(line.rstrip())

    gfa_out.to_file(out_gfa)
    return 0

if __name__ == "__main__":
    exit(main())
