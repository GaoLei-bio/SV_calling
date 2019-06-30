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

parser = argparse.ArgumentParser()
parser.add_argument("--bed_file", type=str, help="Assemblytics output SV file, bed format", required=True, default="")
parser.add_argument("--prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

bed_file = args.bed_file
prefix = args.prefix


''' Function '''
def count_sth(Key,Dict):
    if Key not in Dict:
        Dict[Key] = 1
    else:
        Dict[Key] += 1


''' Main Program'''
redundancy_count = 0
GapSV_count = 0
Close_to_AnchorEnd = 0
fail_blast = 0
pass_blast = 0
No_Good_Pair = 0

found_ref_coord = {}
found_qry_coord = {}

with open(bed_file) as infile:
    for line in infile:
        line = line.strip()
        if not line.startswith("#"):
            fields = line.split("\t")
            ref_coord = "\t".join(fields[0:3])
            count_sth(ref_coord,found_ref_coord)
            qry_coord = re.sub(':[+-]', '', fields[9])
            count_sth(qry_coord,found_qry_coord)

with open(bed_file) as infile, \
     open(prefix + ".filter.all", "w") as whole_file, \
     open(prefix + ".filter.pass", "w") as pass_file, \
     open(prefix + ".filter.fail", "w") as fail_file, \
     open(prefix + ".filter.summary", "w") as summary_file, \
     open(prefix + ".clean_SV.bed", "w") as clean_file:
    for line in infile:
        line = line.strip()
        fields = line.split("\t")
        if line.startswith("#"):
            whole_file.write(line + "\n")
            pass_file.write(line + "\n")
            fail_file.write(line + "\n")
            clean_file.write("\t".join(fields[:11]) + "\n")
        else:
            ref_coord = "\t".join(fields[0:3])
            qry_coord = re.sub(':[+-]', '', fields[9])
            if found_ref_coord[ref_coord] > 1 and found_qry_coord[qry_coord] > 1:
                fields[-1] = "Both_end_Redundant"
                redundancy_count += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif found_ref_coord[ref_coord] > 1:
                fields[-1] = "Ref_Redundant"
                redundancy_count += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif found_qry_coord[qry_coord] > 1:
                fields[-1] = "Qry_Redundant"
                redundancy_count += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif "Gap" in fields[-1]:
                GapSV_count += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif "To_AnchorEnd" in fields[-1]:
                Close_to_AnchorEnd += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif "flank_fail" in fields[-1]:
                fail_blast += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")
            elif fields[-1] == "Pass":
                # pass filter
                pass_blast += 1
                whole_file.write("\t".join(fields) + "\n")
                pass_file.write("\t".join(fields) + "\n")
                clean_file.write("\t".join(fields[:11]) + "\n")
            else:
                No_Good_Pair += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")


    Total_SV = redundancy_count + GapSV_count + Close_to_AnchorEnd + fail_blast + pass_blast + No_Good_Pair
    summary_file.write("A total of " + str(Total_SV) + " raw SVs\n")
    summary_file.write("The following SVs were removed in order:\n")
    summary_file.write("\t1. Redundant SVs:\t" + str(redundancy_count) + "\n")
    summary_file.write("\t2. SVs spanning or close to gaps:\t" + str(GapSV_count) + "\n")
    summary_file.write("\t3. within_alignment SVs are very close to anchor end:\t" + str(Close_to_AnchorEnd) + "\n")
    summary_file.write("\t4. between_alignments SVs without unique anchor pair:\t" + str(No_Good_Pair) + "\n")
    summary_file.write("\t5. SVs without blast supports:\t" + str(fail_blast) + "\n")

    summary_file.write("Finally, SVs passed filtering:\t" + str(pass_blast) + "\n")













































#srf
