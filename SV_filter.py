"""
python SV_filter.py [-h] --bed_file BED_FILE --anchor_file ANCHOR_FILE
                    --ref_genome REF_GENOME --query_genome QUERY_GENOME
                    --gap_flank GAP_FLANK --blast_flank BLAST_FLANK --temp_dir
                    TEMP_DIR --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --bed_file BED_FILE   Assemblytics output SV file, bed format
  --anchor_file ANCHOR_FILE
                        Assemblytics sorted unique anchor file, tab format.
                        e.g.
  --ref_genome REF_GENOME
                        Reference genome file, fasta format. Indexed by
                        samtools
  --query_genome QUERY_GENOME
                        Query species genome file, fasta format. Indexed by
                        samtools
  --gap_flank GAP_FLANK
                        Require SVs are not very close to Gap. This is the
                        size of flank regions not allowed with gaps. For
                        within_alignment indels, this is the minimal distance
                        to anchor end.
  --blast_flank BLAST_FLANK
                        The size for flanking regions for blast check.
  --temp_dir TEMP_DIR   Temporary directory for flanking fasta seqs.
  --prefix PREFIX       Prefix for output files.



"""
import argparse
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("--bed_file", type=str, help="Assemblytics output SV file, bed format", required=True, default="")
parser.add_argument("--anchor_file", type=str, help="Assemblytics sorted unique anchor file, tab format. e.g. ", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query species genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--gap_flank", type=str, help="Require SVs are not very close to Gap. This is the size of flank regions not allowed with gaps. For within_alignment indels, this is the minimal distance to anchor end.", required=True, default="")
parser.add_argument("--blast_flank", type=str, help="The size for flanking regions for blast check.", required=True, default="")
parser.add_argument("--temp_dir", type=str, help="Temporary directory for flanking fasta seqs.", required=True, default="")
parser.add_argument("--prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

bed_file = args.bed_file
anchor_file = args.anchor_file
ref_genome = args.ref_genome
query_genome = args.query_genome
gap_flank = int(args.gap_flank)
blast_flank = int(args.blast_flank)
temp_dir = args.temp_dir
prefix = args.prefix

blast_parameter = "  -perc_identity 90 -dust no -num_threads 1 -outfmt '7 qseqid sseqid qlen slen length qstart qend sstart send sstrand pident nident mismatch gapopen gaps qcovhsp qcovs score bitscore evalue'  -evalue 1e-5 -gapopen 5 -gapextend 3 "

# if temp dir does not exist, creat it
file_list = set(subprocess.check_output("ls", shell=True).decode().split("\n"))
if temp_dir not in file_list:
    subprocess.check_output("mkdir " + temp_dir, shell=True)

''' Functions'''
def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit():
            fields[i] = int(fields[i])
    return fields

def add_to_dict(key,value,dict):
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = []
        dict[key].append(value)

def gap_count(coord,genome):
    ext_seq = "".join(subprocess.check_output("samtools faidx " + genome + " " +  coord, shell=True).decode().split("\n")[1:-1]).upper()
    N_num = ext_seq.count("N")
    return N_num

#def find_host_anchor(ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,anchor_dict):
#    host_anchor = []
#    max_distance = 0
#    for cells in anchor_dict[ref_chr + "-" + qry_chr]:
#        if cells[0] <= ref_sta <= ref_end <= cells[1] and min(cells[2],cells[3]) <= min(qry_sta,qry_end) <= max(qry_sta,qry_end) <= max(cells[2],cells[3]):
#            distance = min(ref_sta - cells[0], cells[1] - ref_end, min(qry_sta,qry_end) - min(cells[2],cells[3]), max(cells[2],cells[3]) - max(qry_sta,qry_end))
#            if distance > max_distance:
#                max_distance = distance
#                host_anchor = [max_distance,cells]
#    return host_anchor

def SV_2_AnchorEnd(ref_sta,ref_end,qry_sta,qry_end,anchor_info):
    anchor_ref_sta = anchor_info[0]
    anchor_ref_end = anchor_info[1]
    anchor_qry_sta = min(anchor_info[2],anchor_info[3])
    anchor_qry_end = max(anchor_info[2],anchor_info[3])
    distance2end = min( ref_sta - anchor_ref_sta, anchor_ref_end - ref_end, min(qry_sta,qry_end) - anchor_qry_sta, anchor_qry_end - max(qry_sta,qry_end))
    return distance2end


def find_neigbor_anchor(ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,anchor_dict):
    neigbor_anchor = []
    left_neigbor = []
    right_neigbor = []
    for cells in anchor_dict[ref_chr + "-" + qry_chr]:
        if ref_sta in set(cells[0:2]) and (qry_sta in set(cells[2:4]) or qry_end in set(cells[2:4])):
            left_neigbor.append(cells)

        elif ref_end in set(cells[0:2]) and (qry_sta in set(cells[2:4]) or qry_end in set(cells[2:4])):
            right_neigbor.append(cells)

    if len(left_neigbor) == 1 and len(right_neigbor) == 1:
        neigbor_anchor = [left_neigbor,right_neigbor]
    elif len(left_neigbor) == 1 and len(right_neigbor) > 1:
        good_right = get_best_pair_by_one_side(left_neigbor,right_neigbor)
        neigbor_anchor = [left_neigbor,good_right]
    elif len(left_neigbor) > 1 and len(right_neigbor) == 1:
        good_left = get_best_pair_by_one_side(right_neigbor,left_neigbor)
        neigbor_anchor = [good_left,right_neigbor]
    elif len(left_neigbor) == 2 and len(right_neigbor) == 0:
        good_left = []
        good_right = []
        for cells in left_neigbor:
            if qry_end in set(cells[2:4]):
                good_left.append(cells)
            elif qry_sta in set(cells[2:4]):
                good_right.append(cells)
        if len(good_left) == 1 and len(good_right) == 1:
            neigbor_anchor = [good_left,good_right]
        else:
            good_left = []
            good_right = []
            for cells in left_neigbor:
                if ref_end in set(cells[0:2]):
                    good_right.append(cells)
                else:
                    good_left.append(cells)
            if len(good_left) == 1 and len(good_right) == 1:
                neigbor_anchor = [good_left,good_right]
            else:
                neigbor_anchor = [left_neigbor,right_neigbor]
    else:
        neigbor_anchor = [left_neigbor,right_neigbor]

    return neigbor_anchor

def get_best_pair_by_one_side(left_neigbor,right_neigbor):
    # left_neigbor is only 1
    # Find one right_neigbor with same direction as left
    good_right = []
    if left_neigbor[0][3] - left_neigbor[0][2] > 0:
        for cells in right_neigbor:
            if cells[3] - cells[2] > 0:
                good_right.append(cells)
    else:
        for cells in right_neigbor:
            if cells[3] - cells[2] < 0:
                good_right.append(cells)
    return good_right


def confirm_neigbor_pair(neigbor_anchor,line):
    left_neigbor = neigbor_anchor[0][0]
    right_neigbor = neigbor_anchor[1][0]
    if left_neigbor[3] > left_neigbor[2]:
        left_str = "plus"
    else:
        left_str = "minus"
    if right_neigbor[3] > right_neigbor[2]:
        right_str = "plus"
    else:
        right_str = "minus"

    if left_str == right_str == "plus":
        if left_neigbor[2] > right_neigbor[2]:
            left_neigbor = neigbor_anchor[1][0]
            right_neigbor = neigbor_anchor[0][0]
        elif left_neigbor[2] == right_neigbor[2]:
            print ("#" *20)
            print (line + "\t" + Note)
            print ("#left_neigbor")
            print ("\t".join(map(str,left_neigbor)))
            print ("#right_neigbor")
            print ("\t".join(map(str,right_neigbor)))
    elif left_str == right_str == "minus":
        if left_neigbor[2] < right_neigbor[2]:
            left_neigbor = neigbor_anchor[1][0]
            right_neigbor = neigbor_anchor[0][0]
        elif left_neigbor[2] == right_neigbor[2]:
            print ("#" *20)
            print (line + "\t" + Note)
            print ("#left_neigbor")
            print ("\t".join(map(str,left_neigbor)))
            print ("#right_neigbor")
            print ("\t".join(map(str,right_neigbor)))
    else:
        Note = str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])) + ":" + left_str + "|" + right_str + ":Diff_direction_anchor"

    if left_str == right_str == "plus":
        if left_neigbor[2] < right_neigbor[2]:
            Note = str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])) + ":" + left_str + "|" + right_str + ":Anchor_Order_good"
        else:
            Note = str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])) + ":" + left_str + "|" + right_str + ":Anchor_Order_bad"
    elif left_str == right_str == "minus":
        if left_neigbor[2] < right_neigbor[2]:
            Note = str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])) + ":" + left_str + "|" + right_str + ":Anchor_Order_bad"
        else:
            Note = str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])) + ":" + left_str + "|" + right_str + ":Anchor_Order_good"


    if "Anchor_Order_bad" in Note:
        print ("@" *20)
        print (line + "\t" + Note)
        print ("#left_neigbor")
        print ("\t".join(map(str,left_neigbor)))
        print ("#right_neigbor")
        print ("\t".join(map(str,right_neigbor)))
    return left_neigbor,right_neigbor,Note



def within_alignment_blast_check(sv_id,ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,anchor_info,anchor_str):
    blast_result = "Both_flank_fail"
    left_check = "Fail"
    right_check = "Fail"
    # extract ref left and right, blast to query genome
    # if fail, extract qry left and right, blast to ref genome
    anchor_ref_sta = anchor_info[0]
    anchor_ref_end = anchor_info[1]
    anchor_qry_sta = min(anchor_info[2],anchor_info[3])
    anchor_qry_end = max(anchor_info[2],anchor_info[3])
    # check left seq
    if ref_sta > gap_flank:
        left_region = ref_chr + ":" + str(max(ref_sta-blast_flank+1, anchor_ref_sta)) + "-" + str(ref_sta)
        left_seq = re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + ref_chr + ".fasta " +  left_region, shell=True).decode().split("\n")[1:-1]))
        with open(temp_dir + "/" + sv_id + ".ref.left.fa", 'w') as fa_file:
            fa_file.write(">" + sv_id + "_L_" + left_region + "\n")
            fa_file.write(left_seq + "\n")
        # blast
        blast_outputs = subprocess.check_output("blastn -query " + temp_dir + "/" + sv_id + ".ref.left.fa " + " -db " +  qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
        if anchor_str == "plus":
            exp_hit_sites = set(range(qry_sta - 3, qry_sta + 4))
        else:
            exp_hit_sites = set(range(qry_end - 3, qry_end + 4))

        for line in blast_outputs:
            if not line.startswith("#"):
                cells = convert_int(line.split("\t"))
                if cells[1] == qry_chr and cells[9] == anchor_str and cells[4] >= 50 and cells[2] - cells[6] < 500:
                    # same strand, end 500
                    if anchor_str == "plus":
                        exp_end = cells[8] + (cells[2] - cells[6])
                    else:
                        exp_end = cells[8] - (cells[2] - cells[6])

                    if exp_end in exp_hit_sites:
                        left_check = "Pass"
                        break
    # check right seq
    right_region = ref_chr + ":" + str(ref_end) + "-" + str(min(ref_end + blast_flank - 1,anchor_ref_end))
    right_seq = re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + ref_chr + ".fasta " +  right_region, shell=True).decode().split("\n")[1:-1]))
    with open(temp_dir + "/" + sv_id + ".ref.right.fa", 'w') as fa_file:
        fa_file.write(">" + sv_id + "_L_" + right_region + "\n")
        fa_file.write(right_seq + "\n")
    # blast
    blast_outputs = subprocess.check_output("blastn -query " + temp_dir + "/" + sv_id + ".ref.right.fa " + " -db " +  qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
    if anchor_str == "plus":
        exp_hit_sites = set(range(qry_end - 3, qry_end + 4))
    else:
        exp_hit_sites = set(range(qry_sta - 3, qry_sta + 4))
    for line in blast_outputs:
        if not line.startswith("#"):
            cells = convert_int(line.split("\t"))
            if cells[1] == qry_chr and cells[9] == anchor_str and cells[4] >= 50 and cells[5] < 500:
                # same strand, end 500
                if anchor_str == "plus":
                    exp_sta = cells[7] - (cells[5] - 1)
                else:
                    exp_sta = cells[7] + (cells[5] - 1)
                if exp_sta in exp_hit_sites:
                    right_check = "Pass"
                    break
    # get result
    if left_check == "Pass" and right_check == "Pass":
        blast_result = "Pass"
    elif left_check == "Pass":
        blast_result = "Right_flank_fail"
    elif right_check == "Pass":
        blast_result = "Left_flank_fail"
    return blast_result


def between_alignment_blast_check(sv_id,ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,qry_direct,left_neigbor,right_neigbor,Note):
    blast_result = "Both_flank_fail"
    left_check = "Fail"
    right_check = "Fail"
    if "plus" in Note:
        anchor_str = "plus"
        left_exp_sits = set(range(qry_sta - 3, qry_sta + 4))
        right_exp_sits = set(range(qry_end - 3, qry_end + 4))
    else:
        anchor_str = "minus"
        left_exp_sits = set(range(qry_end - 3, qry_end + 4))
        right_exp_sits = set(range(qry_sta - 3, qry_sta + 4))


    # left seq
    if ref_sta > gap_flank:
        anchor_ref_sta = min(left_neigbor[0],right_neigbor[0])
        left_region = ref_chr + ":" + str(max(ref_sta-blast_flank+1, anchor_ref_sta,gap_flank)) + "-" + str(ref_sta)
        left_seq = re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + ref_chr + ".fasta " +  left_region, shell=True).decode().split("\n")[1:-1]))
        with open(temp_dir + "/" + sv_id + ".ref.left.fa", 'w') as fa_file:
            fa_file.write(">" + sv_id + "_L_" + left_region + "\n")
            fa_file.write(left_seq + "\n")
        left_blast = subprocess.check_output("blastn -query " + temp_dir + "/" + sv_id + ".ref.left.fa " + " -db " +  qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
        # check hits
        for line in left_blast:
            if not line.startswith("#"):
                cells = convert_int(line.split("\t"))
                if cells[1] == qry_chr and cells[9] == anchor_str and cells[4] >= 50 and cells[2] - cells[6] < 500:
                    # same strand, end 500
                    if anchor_str == "plus":
                        exp_end = cells[8] + (cells[2] - cells[6])
                    else:
                        exp_end = cells[8] - (cells[2] - cells[6])
                    if exp_end in left_exp_sits:
                        left_check = "Pass"
                        break


    # right_seq
    anchor_ref_end = max(left_neigbor[1],right_neigbor[1])
    right_region = ref_chr + ":" + str(ref_end) + "-" + str(min(ref_end + blast_flank - 1,anchor_ref_end))
    right_seq = re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + ref_chr + ".fasta " +  right_region, shell=True).decode().split("\n")[1:-1]))
    with open(temp_dir + "/" + sv_id + ".ref.right.fa", 'w') as fa_file:
        fa_file.write(">" + sv_id + "_L_" + right_region + "\n")
        fa_file.write(right_seq + "\n")
    right_blast = subprocess.check_output("blastn -query " + temp_dir + "/" + sv_id + ".ref.right.fa " + " -db " +  qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
    # right
    for line in right_blast:
        if not line.startswith("#"):
            cells = convert_int(line.split("\t"))
            if cells[1] == qry_chr and cells[9] == anchor_str and cells[4] >= 50 and cells[5] < 500:
                # same strand, end 500
                if anchor_str == "plus":
                    exp_sta = cells[7] - (cells[5] - 1)
                else:
                    exp_sta = cells[7] + (cells[5] - 1)
                if exp_sta in right_exp_sits:
                    right_check = "Pass"
                    break

    # get result
    if left_check == "Pass" and right_check == "Pass":
        blast_result = "Pass"
    elif left_check == "Pass":
        blast_result = "Right_flank_fail"
    elif right_check == "Pass":
        blast_result = "Left_flank_fail"
    return blast_result


''' Get anchors '''
anchor_dict = {}

with open(anchor_file) as infile:
    for line in infile:
        line = line.strip()
        cells = convert_int(line.split("\t"))
        ref_chr = cells[6]
        qry_chr = cells[7]
        add_to_dict(ref_chr + "-" + qry_chr,cells,anchor_dict)

''' Main Program'''
noGapSV_count = 0
GapSV_count = 0
Close_to_AnchorEnd = 0
No_Good_Pair = 0

checked_ref_coord = set()
checked_qry_coord = set()
redundancy_count = 0

pass_blast = 0
fail_blast = 0


with open(bed_file) as infile, open(prefix + ".all_checked.bed", "w") as whole_file:
    for line in infile:
        line = line.strip()
        if line.startswith("#"):
            whole_file.write(line + "\tNote\n")
        else:
            fields = convert_int(line.split("\t"))
            sv_type = fields[6]
            sv_id = fields[3]
            ref_chr = fields[0]
            ref_sta = fields[1]
            ref_end = fields[2]
            ref_size = fields[7]
            ref_coord = fields[0] + ":" + str(max(ref_sta-gap_flank, 1)) + "-" + str(ref_end + gap_flank)
            qry_chr = fields[9].split(":")[0]
            qry_sta = int(fields[9].split(":")[1].split("-")[0])
            qry_end = int(fields[9].split(":")[1].split("-")[1])
            qry_size = fields[8]
            qry_coord = qry_chr + ":" + str(max(qry_sta - gap_flank,1)) + "-" + str(qry_end + gap_flank)
            qry_direct = fields[9].split(":")[2]
            if ref_coord in checked_ref_coord and qry_coord in checked_qry_coord:
                # 1st, remove redundant SVs
                Note = "Both_end_Redundant"
                redundancy_count += 1
                whole_file.write(line + "\t" + Note + "\n")
            elif ref_coord in checked_ref_coord:
                Note = "Ref_Redundant"
                redundancy_count += 1
                whole_file.write(line + "\t" + Note + "\n")
            elif qry_coord in checked_qry_coord:
                Note = "Qry_Redundant"
                redundancy_count += 1
                whole_file.write(line + "\t" + Note + "\n")
            else:
                checked_ref_coord.add(ref_coord)
                checked_qry_coord.add(qry_coord)
                ref_gap = gap_count(ref_coord,ref_genome)
                qry_gap = gap_count(qry_coord,query_genome)
                if ref_gap == 0 and qry_gap == 0:
                    # 2nd, no gap in both
                    # blast check flanking regions
                    noGapSV_count += 1
                    if fields[10] == "within_alignment":
                        anchor_info = fields[12:]
                        if fields[11] == "Forward":
                            anchor_str = "plus"
                        else:
                            anchor_str = "minus"

                        distance2end = SV_2_AnchorEnd(ref_sta,ref_end,qry_sta,qry_end,anchor_info)
                        if distance2end < gap_flank:
                            Close_to_AnchorEnd += 1
                            Note = "To_AnchorEnd_" + str(distance2end) + "bp"
                            whole_file.write(line + "\t" + Note + "\n")
                        else:
                            # SV is far enough to anchor end
                            blast_result = within_alignment_blast_check(sv_id,ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,anchor_info,anchor_str)
                            outline = line + "\t" + blast_result +  "\n"
                            whole_file.write(outline)
                            if blast_result == "Pass":
                                pass_blast += 1
                            else:
                                fail_blast += 1

                    else:
                        neigbor_anchor = find_neigbor_anchor(ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,anchor_dict)
                        if not (len(neigbor_anchor[0]) == 1 and len(neigbor_anchor[1]) == 1):
                            # if cann't get 1 unique pair of anchor, output information for both SV and potential anchors for manually check
                            # should no output
                            No_Good_Pair += 1
                            print ("#" *20)
                            print (str(len(neigbor_anchor[0])) + "|" + str(len(neigbor_anchor[1])))
                            print (line)
                            print ("#left_neigbor")
                            print (neigbor_anchor[0])
                            print ("#right_neigbor")
                            print (neigbor_anchor[1])
                            Note = ("No_Good_anchor_pair")
                            whole_file.write(line + "\tNA" * 9 + "\t" + Note + "\n")
                        else:
                            left_neigbor,right_neigbor,Note = confirm_neigbor_pair(neigbor_anchor,line)
                            if "Diff_direction_anchor" in Note:
                                No_Good_Pair += 1
                                print ("#" *20)
                                print (line)
                                print (Note)
                                print ("#left_neigbor")
                                print (left_neigbor)
                                print ("#right_neigbor")
                                print (right_neigbor)
                                whole_file.write(line + "\tNA" * 9 + "\t" + "Diff_direction_anchor" + "\n")
                            else:
                                blast_result = between_alignment_blast_check(sv_id,ref_chr,ref_sta,ref_end,qry_chr,qry_sta,qry_end,qry_direct,left_neigbor,right_neigbor,Note)
                                if blast_result == "Pass":
                                    pass_blast += 1
                                else:
                                    fail_blast += 1
                                # get anchor pair
                                if "plus" in Note:
                                    anchor_pair = ["Forward"]
                                else:
                                    anchor_pair = ["Reverse"]
                                for i in range(0,len(left_neigbor)):
                                    anchor_pair.append(str(left_neigbor[i]) + "|" + str(right_neigbor[i]))

                                outline = line + "\t" + "\t".join(anchor_pair) + "\t" + blast_result +  "\n"
                                whole_file.write(outline)

                else:
                    GapSV_count += 1
                    if ref_gap > 0 and qry_gap > 0:
                        Note = "Both_Gap"
                    elif ref_gap > 0:
                        Note = "Ref_Gap"
                    elif qry_gap > 0:
                        Note = "Qry_Gap"
                    if fields[10] == "within_alignment":
                        whole_file.write(line + "\t" + Note + "\n")
                    else:
                        whole_file.write(line + "\tNA" * 9 + "\t" + Note + "\n")







# log
# log
print ("#" * 30)
print ("Input SV file:\t" + bed_file)
print ("Reference genome file:\t" + ref_genome)
print ("Query genome file:\t" + query_genome)
print ("Input SVs:\t" + str(redundancy_count + noGapSV_count + GapSV_count))
print ("The following SVs were removed in order:")
print ("\t1. Redundant SVs:\t" + str(redundancy_count))
print ("\t2. SVs spanning or close to gaps:\t" + str(GapSV_count))
print ("\t3. within_alignment SVs are very close to anchor end:\t" + str(Close_to_AnchorEnd))
print ("\t4. between_alignments SVs without unique anchor pair:\t" + str(No_Good_Pair))
print ("\t5. SVs without blast supports:\t" + str(fail_blast))

print ("Finally, SVs passed filtering:\t" + str(pass_blast))





























#srf
