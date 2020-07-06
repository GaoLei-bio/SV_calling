"""

ref_start,ref_end,query_start,query_end,ref_length,query_length,ref,query,tag,alignment_length
24189148,25248284,27054200,28111716,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448
29303211,30236628,31714169,32647582,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448
21780547,22455127,24841632,25516203,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448
32730544,33389287,35113415,35788525,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448
30934593,31517123,33318176,33900568,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448
23619973,24197977,26478099,27056087,51253844,53272422,OXv7ch01,KOv7ch01,unique,1448

"""
import sys
from operator import itemgetter

def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

i = -1
my_table = []
with open(sys.argv[1]) as infile:
    for line in infile:
        cells = convert_int(line.strip().replace("ch00","chXX").split(","))
        i += 1
        if i == 0:
            print "\t".join(map(str,cells))
        else:
            my_table.append(cells)


sorted_table = sorted(my_table, key = itemgetter(6,7,0,1,2,3))

for cells in sorted_table:
    print "\t".join(map(str,cells)).replace("chXX","ch00")





#srf
