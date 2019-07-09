####################################################################
# Version 1.0
# 06/27/2019
# This pipeline is used to detect structural variations (SVs) between 2 genomes.
# Step 1. Genome alignment 				Using minimap2
# Step 2. Unique anchor detection		Using Assemblytics implemented in RaGOO
# Step 3. SV calling					Using Assemblytics implemented in RaGOO for SVs between anchors, and my script for SVs within anchors
# Step 4. Filtering						My script

####################################################################
# Command:
#      bash SV_pipeline.sh  <Reference_genome_file>   <Query_genome_file>  <Prefix_for_outputs> <number of threads>

####################################################################
# Dependencies:
#	1. minimap2		For genome alignment
#	2. RaGOO        For unique anchor detection and SV calling between anchors
#	3. samtools
#	4. blast+
#	5. Python 3

 
####################################################################
# Genomes should be provided as fasta format and indexed by samtools

date 
echo "Now, prepare genome files..."

reference=$1
query=$2
prefix=$3
cpuN=$4

samtools faidx $query
samtools faidx $reference

#makeblastdb -dbtype nucl -in $query
#makeblastdb -dbtype nucl -in $reference

# To speed up blast for SV filtering, make blast database for each chr 
for chr in `cut -f 1 "$reference".fai`
do
	samtools faidx $reference $chr > $chr.fasta 
	samtools faidx $chr.fasta
	makeblastdb -dbtype nucl -in $chr.fasta
done

for chr in `cut -f 1 "$query".fai`
do
	samtools faidx $query $chr > $chr.fasta 
	samtools faidx $chr.fasta
	makeblastdb -dbtype nucl -in $chr.fasta
done


####################################################################
# Step 1. Genome alignment
#  	minimap2 (https://github.com/lh3/minimap2)
echo "#####################################"
date 
echo "Step 1, genome alignment"

minimap2 -ax asm5 --cs -t$cpuN  $reference $query > pm_against_ref.sam 2> pm_contigs_against_ref.sam.log

# -x STR
#    asm5  | sequence divergence below 1%.
#    asm10 | for divergence around a couple of percent
#    asm20 | for divergence not more than 10%.

samtools view -bS pm_against_ref.sam | samtools sort --threads $cpuN -o pm_against_ref.bam /dev/stdin
samtools index pm_against_ref.bam
####################################################################
# Step 2. Unique anchor detection
# Step 2.1 Convert sam file to delta file
#          sam2delta.py is a python script from RaGOO
echo "#####################################"
date 
echo "Step 2, get unique alignments"

#source activate test_py3 

echo "Activate python 3"

sam2delta.py  pm_against_ref.sam   # Output pm_against_ref.sam.delta

# Step 2.2 Get unique alignments as anchors for SV calling
#          Assemblytics_uniq_anchor.py is a python script from RaGOO
Assemblytics_uniq_anchor.py --delta pm_against_ref.sam.delta \
                            --unique-length 10000 \
							--out $prefix \
							--keep-small-uniques

####################################################################
# Step 3. SV calling
# Step 3.1 Detect SVs between unique anchors 
#          Assemblytics_between_alignments.pl is a perl script from RaGOO
#          The following script will call SV with a minimal length of 10 bp and a maximal size of 10000 bp 
echo "#####################################"
date 
echo "Step 3. SV calling"
echo "Step 3.1 Detect SVs between unique anchors "

echo -e "#reference\tref_start\tref_stop\tID\tsize\tstrand\ttype\tref_gap_size\tquery_gap_size\tquery_coordinates\tmethod" > $prefix.variants_between_alignments.bed
Assemblytics_between_alignments.pl $prefix.coords.tab 10 10000 all-chromosomes exclude-longrange bed >> $prefix.variants_between_alignments.bed

# Step 3.2 Detect SVs within unique alignment
#          Assemblytics can also call indels within alignments. However, the coordinates on reference genome for those detected on reverse alignment are incorrect.
#    	   The following script will call indels within alignments based on the sam file produced by minmap2 (Step 1), and unique anchors identified by Assemblytics (Step 2.2)
#          The SVs detected on forward alignments are compatible with those identified by Assemblytics. 
#          For those detected on reverse alignments, the SV number and coordinates on query genome are generally consistent with Assemblytics results. The SV coordinates on reference genome are corrected.

echo "Step 3.2 Detect SVs within unique alignment"

# run Assemblytics_within_alignment.py
# Get this result in case for future check
Assemblytics_within_alignment.py --delta $prefix.Assemblytics.unique_length_filtered_l10000.delta --min 10 > assemblytics_out.variants_within_alignments.bed

# run my indel calling script
python Call_Indel_within_alignment.py --sam_file    pm_against_ref.sam \
                              --ref_genome   $reference \
							  --query_genome $query \
							  --uniqe_anchor $prefix.coords.tab \
							  --min_size     10  \
							  --max_size     100000 \
							  --prefix       $prefix  >  $prefix.SV_within_anchor.log
# Main output: $prefix.variants_within_alignments.bed
#    column 1-11,  SV information, according to the format of assemblytics_out.variants_within_alignments.bed
#	 column 12     anchor direction
#    column 13-20, anchor information, according to $prefix.coords.tab format. This can be directly used for SV filtering.

python Get_unique_anchor_sam_file.py --sam_file  pm_against_ref.sam   --unique_anchor  $prefix.sam_1-6_col.unique  > $prefix.unique.sam

samtools view -bS $prefix.unique.sam | samtools sort --threads $cpuN -o $prefix.unique.bam /dev/stdin
samtools index $prefix.unique.bam
										  
####################################################################
# Step 4. SV filtering
# The following SVs will be removed in order:
#    1. Redundant SVs
#    2. SVs spanning or close to gaps (50bp) 
#    3. within_alignment SVs are very close to anchor end (50bp, same as minimal distance to gap)
#    4. between_alignments SVs without unique anchor pair (Should be no this type)
#    5. SVs without blast supports of at least 1 flanking region
echo "#####################################"
date 
echo "Step 4. SV filtering"

cat  $prefix.variants_within_alignments.bed > $prefix.rawSV.bed
grep -v ^"#" $prefix.variants_between_alignments.bed >> $prefix.rawSV.bed 

lineNum=`cat "$prefix".rawSV.bed | wc -l| awk '{print int($1/"'$cpuN'")+1}'`

split -l $lineNum $prefix.rawSV.bed -d -a 4 Split_


# Step 4.1 SV filtering 
for file in `ls Split_????`
do

python SV_filter.py         --bed_file     $file \
                            --anchor_file  $prefix.coords.tab \
                            --ref_genome   $reference \
							--query_genome $query \
							--gap_flank    50 \
							--blast_flank  5000 \
							--temp_dir     ${file}_temp_dir      \
							--prefix       $file  > $file.log  &

done 

wait

cat Split_????.all_checked.bed   > ${prefix}.checked.all

rm Split_????* -rf

# Step 4.2 Get summary of SV filtering results 
# Further remove redundant SVs
python SV_filter_summary.py --bed_file  ${prefix}.checked.all  --prefix   ${prefix}
# Outputs:
#     ${prefix}.clean_SV.bed       Clean SV list, Assemblytics format. This is the final result
#     ${prefix}.filter.all         All SVs with filtering detail
#     ${prefix}.filter.temp_pass   with anchor information
#     ${prefix}.filter.fail        with fail reason
#     ${prefix}.filter.summary     Summary for SV filtering 
# rm ${prefix}.checked.all

rm temp_dir -rf


lineNum=`cat "$prefix".filter.temp_pass | wc -l| awk '{print int($1/"'$cpuN'")+1}'`

split -l $lineNum $prefix.filter.temp_pass -d -a 4 Split_

for file in `ls Split_????`
do
	python SV_flank_depth_check.py --bed_file  ${file}  --bam_file   $prefix.unique.bam  &
done

wait

pass_uniqe=`cat Split_????.unique_check.pass | wc -l| awk '{print $1-1}'`

echo -e "\t"$pass_uniqe" passed unique anchor check" >> ${prefix}.filter.summary

cat  Split_????.unique_check.pass  >  ${prefix}.filter.pass 
cat  Split_????.unique_check.fail  >> ${prefix}.filter.fail
cat  Split_????.clean_SV.bed       > ${prefix}.clean_SV.bed


####################################################################
# Step 5. Clean temporary files

echo "#####################################"
date 
echo "Step 6. Clean temporary files"



for chr in `cut -f 1 "$reference".fai "$query".fai`
do
	rm $chr.fasta*	
done


#source deactivate


