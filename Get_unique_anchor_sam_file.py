"""
usage: call_indel_from_sam.py [-h] --sam_file SAM_FILE --ref_genome REF_GENOME
                              --query_genome QUERY_GENOME --uniqe_anchor
                              UNIQE_ANCHOR --min_size MIN_SIZE --max_size
                              MAX_SIZE --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --sam_file SAM_FILE   Input sam file. Produced by minmap2.
  --ref_genome REF_GENOME
                        Reference genome file, fasta format. Indexed by
                        samtools
  --query_genome QUERY_GENOME
                        Query species genome file, fasta format. Indexed by
                        samtools
  --uniqe_anchor UNIQE_ANCHOR
                        Unique anchors selected by RaGOO
                        (Assemblytics_uniq_anchor.py). We can only call Indel
                        based on unique anchors.
  --min_size MIN_SIZE   Minimal indel size
  --max_size MAX_SIZE   Maximal indel size
  --prefix PREFIX       Prefix for output files.


"""
import argparse
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("--sam_file", type=str, help="Input sam file. Produced by minmap2.", required=True, default="")
parser.add_argument("--unique_anchor", type=str, help="Col_1-6 for unique anchors", required=True, default="")

args = parser.parse_args()

sam_file = args.sam_file
unique_anchor = args.unique_anchor

''' get anchor infor '''

U_anchor_infor = []

with open(unique_anchor) as infile:
    for line in infile:
        line = line.strip()
        U_anchor_infor.append(line)

U_anchor_infor = set(U_anchor_infor)

''' Filter sam file '''

with open(sam_file) as infile:
    for line in infile:
        line = line.strip()
        cells = line.split("\t")
        if line.startswith("@") or "\t".join(cells[0:6]) in U_anchor_infor:
            print (line)
























#srf
