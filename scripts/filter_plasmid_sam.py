import pysam
import sys
from collections import defaultdict
import numpy as np
import argparse

def process_records(afile, records):
    read_length = records[0].infer_read_length()



    mappings = defaultdict(list)
    referenc = defaultdict(list)

    if records[0].is_unmapped:
        return 0
    for r in records:
        qs = r.query_alignment_start
        qe = r.query_alignment_end
        if r.cigartuples[0][0] in {5}: #Clip
            qs += r.cigartuples[0][1]
            qe += r.cigartuples[0][1]

        
        rs = r.reference_start
        re = r.reference_end


        mappings[r.reference_name].append((qs,qe))
        referenc[r.reference_name].append((rs,re))
    mappings_flat = {}
    for k,v in mappings.items():

        vv = sorted(v, key=lambda x:x[0])
        x = [vv[0]]
        for i in vv[1:]:
            if x[-1][1] < i[0]:
                x.append(i)
            elif i[1] > x[-1][1]:
                x[-1] = (x[-1][0], i[1])
        mappings_flat[k] = x
    #print(mappings_flat, lens)
    mapping_percentages = {k: np.sum([v2-v1 for v1,v2 in v ])/read_length for k,v in mappings_flat.items()}
    
    best_mapping_percentage = -1
    for k,v in mapping_percentages.items():
        if best_mapping_percentage < v:
            best_mapping_percentage = v
    
    if best_mapping_percentage > 1 or best_mapping_percentage < 0:
        print(r.query_name, best_mapping_percentage, read_length, mappings_flat, mappings,referenc, sep="\t")
        #for r in records:
        #    print(r.cigartuples)
        #print("------")
    return best_mapping_percentage


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--stats")
    parser.add_argument("-m", "--minimum-mapping-threshold", type=float, default=.90)

    args = parser.parse_args()

    isam = pysam.AlignmentFile(args.input, "r",check_sq=False)
    osam = pysam.AlignmentFile(args.output, "wb", template=isam)
    threshold = args.minimum_mapping_threshold

    prev_records = []

    for record in isam.fetch():
        if prev_records:
            if prev_records[-1].query_name == record.query_name:
                prev_records.append(record)
            else:
                rate = process_records(isam, prev_records)
                #print(rate)
                if rate > threshold:
                    for rr in prev_records:
                        osam.write(rr)
                prev_records = [record]

        else:
            prev_records = [record]

    if prev_records and process_records(isam, prev_records) > threshold:
        #print( process_records(isam, prev_records) )
        for rr in prev_records:
            osam.write(rr)

    isam.close()
    osam.close()
if __name__ == "__main__":
    exit(main())
