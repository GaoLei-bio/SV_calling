####################################################################
# Version 1.0
# 06/27/2019
# This pipeline is used to detect structural variations (SVs) between 2 genomes.
# Step 1. Genome alignment 			Using minimap2
# Step 2. Unique anchor detection		Using Assemblytics implemented in RaGOO
# Step 3. SV calling				Using Assemblytics implemented in RaGOO for SVs between anchors, and my script for SVs within anchors
# Step 4. Filtering				My script

####################################################################
# Dependencies:
#	1. minimap2		For genome alignment
#	2. RaGOO        	For unique anchor detection and SV calling between anchors
#	3. samtools
#	4. blast+
#	5. Python 3

####################################################################
Command:
    bash SV_pipeline.sh  <Reference_genome_file>   <Query_genome_file>  <Prefix_for_outputs> <number of threads>


For example: 
    In current directory, We have 2 genomes:
            a reference genome:  SL4.0.genome.fasta
            a query genome:      Pimp_v1.2.fasta

    We want a prefix for outputs: SP2SL
    Number of CUPs we can use: 32

    

The command is:
    bash SV_pipeline.sh  SL4.0.genome.fasta  Pimp_v1.2.fasta   SP2SL   32 > log.txt

Main outputs:
    prefix.clean_SV.bed       Filtered SV list, Assemblytics format. 





