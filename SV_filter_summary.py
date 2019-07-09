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

def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields


def check_overlap(temp_pass_table):
    overlaped_ids = []
    sort_table = sorted(temp_pass_table, key = itemgetter(0,1,2))
    sv_chr = sort_table[0][0]
    sv_sta = sort_table[0][1]
    sv_end = sort_table[0][2]
    sv_id = sort_table[0][3]
    sv_size = max(sv_sta,sv_end) - min(sv_sta,sv_end)
    for sv in sort_table[1:]:
        this_chr = sv[0]
        this_sta = sv[1]
        this_end = sv[2]
        this_id = sv[3]
        this_size = max(this_sta,this_end) - min(this_sta,this_end)
        if this_chr == sv_chr and sv_sta <= this_sta <= sv_end:
            # overlaping
            overlap_size = len(set(range(min(sv_sta,sv_end),max(sv_sta,sv_end))).intersection(range(min(this_sta,this_end),max(this_sta,this_end))))
            if overlap_size * 100.0 / sv_size > 50:
                overlaped_ids.append(sv_id)
            if overlap_size * 100.0 / this_size > 50:
                overlaped_ids.append(this_id)
        sv_chr = this_chr
        sv_sta = this_sta
        sv_end = this_end
        sv_id = this_id
        sv_size = this_size
    return overlaped_ids

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

temp_pass_ID_list = []
temp_pass_by_ID = {}

temp_pass_by_Ref = []
temp_pass_by_Qry = []

del_group = set(["Deletion","Deletion/Substitution","Repeat_contraction","Tandem_contraction"])
ins_group = set(["Insertion","Insertion/Substitution","Repeat_expansion","Tandem_expansion"])

with open(bed_file) as infile, \
     open(prefix + ".filter.all", "w") as whole_file, \
     open(prefix + ".filter.temp_pass", "w") as pass_file, \
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
                sv_id = fields[3]
                temp_pass_ID_list.append(sv_id)
                temp_pass_by_ID[sv_id] = fields
                fields = convert_int(fields)
                if fields[6] in del_group:
                    temp_pass_by_Ref.append(fields[0:4])
                else:
                    qry_chr = qry_coord.split(":")[0]
                    qry_sta = int(qry_coord.split(":")[1].split("-")[0])
                    qry_end = int(qry_coord.split(":")[1].split("-")[1])
                    temp_pass_by_Qry.append([qry_chr,min(qry_sta,qry_end),max(qry_sta,qry_end),sv_id])
                #whole_file.write("\t".join(fields) + "\n")
                #pass_file.write("\t".join(fields) + "\n")
                #clean_file.write("\t".join(fields[:11]) + "\n")
            else:
                No_Good_Pair += 1
                whole_file.write("\t".join(fields) + "\n")
                fail_file.write("\t".join(fields) + "\n")

    ''' check overlaping '''
    overlaped_by_Ref = check_overlap(temp_pass_by_Ref)
    overlaped_by_Qry = check_overlap(temp_pass_by_Qry)
    overlaped_IDs = set(overlaped_by_Ref + overlaped_by_Qry)
    for sv_id in temp_pass_ID_list:
        fields = temp_pass_by_ID[sv_id]
        if sv_id in overlaped_IDs:
            fields[-1] = "Overlaped"
            whole_file.write("\t".join(map(str,fields)) + "\n")
            fail_file.write("\t".join(map(str,fields)) + "\n")
        else:
            whole_file.write("\t".join(map(str,fields)) + "\n")
            pass_file.write("\t".join(map(str,fields)) + "\n")
            clean_file.write("\t".join(map(str,fields[:11])) + "\n")

    Total_SV = redundancy_count + GapSV_count + Close_to_AnchorEnd + fail_blast + pass_blast + No_Good_Pair
    summary_file.write("A total of " + str(Total_SV) + " raw SVs\n")
    summary_file.write("The following SVs were removed in order:\n")
    summary_file.write("\t1. Redundant SVs:\t" + str(redundancy_count) + "\n")
    summary_file.write("\t2. SVs spanning or close to gaps:\t" + str(GapSV_count) + "\n")
    summary_file.write("\t3. within_alignment SVs are very close to anchor end:\t" + str(Close_to_AnchorEnd) + "\n")
    summary_file.write("\t4. between_alignments SVs without unique anchor pair:\t" + str(No_Good_Pair) + "\n")
    summary_file.write("\t5. SVs without blast supports:\t" + str(fail_blast) + "\n")
    summary_file.write("\t6. At least 50% of SV region is overlapped by others:\t" + str(len(overlaped_IDs)) + "\n")

    summary_file.write("Finally, SVs passed filtering:\t" + str(pass_blast - len(overlaped_IDs)) + "\n")













































#srf
