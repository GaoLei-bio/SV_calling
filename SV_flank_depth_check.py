"""
python SV_filter_summary.py [-h] --bed_file BED_FILE --prefix PREFIX

optional arguments:
  -h, --help           show this help message and exit
  --bed_file BED_FILE  Assemblytics output SV file, bed format
  --prefix PREFIX      Prefix for output files.


"""
import argparse
import subprocess
import re
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("--bed_file", type=str, help="Assemblytics output SV file, bed format", required=True, default="")
parser.add_argument("--bam_file", type=str, help="bam file for unique anchor.", required=True, default="")

args = parser.parse_args()

bed_file = args.bed_file
bam_file = args.bam_file


''' Function '''

def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def check_flank_depth(fields):
    check_out = "Fail"
    chrID = fields[0]
    Sta = fields[1]
    End = fields[2]
    left_region = chrID + ":" + str(Sta-5) + "-" + str(Sta-3)
    left_dep = depth_check(left_region)
    right_region = chrID + ":" + str(End + 3) + "-" + str(End + 5)
    right_dep = depth_check(left_region)
    if left_dep == right_dep == 1:
        check_out = "Pass"
    return check_out

def depth_check(coord):
    depth = 1
    depth_output = subprocess.check_output("samtools depth " + bam_file + " -r " +  coord, shell=True).decode().split("\n")[:-1]
    for line in depth_output:
        cells = convert_int(line.strip().split("\t"))
        if cells[2] != 1:
            depth = cells[2]
    return depth

''' Main Program'''

with open(bed_file) as infile, \
     open(bed_file + ".unique_check.pass", "w") as pass_file, \
     open(bed_file + ".unique_check.fail", "w") as fail_file,  \
     open(bed_file + ".clean_SV.bed", "w") as clean_file:
    for line in infile:
        line = line.strip()
        fields = convert_int(line.split("\t"))
        if line.startswith("#"):
            pass_file.write(line + "\n")
            fail_file.write(line + "\n")
            clean_file.write("\t".join(fields[:11]) + "\n")
        else:            
            flank_depth = check_flank_depth(fields)
            if flank_depth == "Fail":
                fields[-1] = "Non_unique_region"
                fail_file.write("\t".join(map(str,fields)) + "\n")
            else:
                pass_file.write("\t".join(map(str,fields)) + "\n")
                clean_file.write("\t".join(map(str,fields[:11])) + "\n")

















































#srf
