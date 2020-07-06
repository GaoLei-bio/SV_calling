####################################################################
SV calling scripts
By Lei Gao
Version 1.1
Date: 07/04/2020
This pipeline is designed to detect structural variations (SVs) between high-quality genomes of 2 closely related species/variaties.
One is called Reference genome, the other one is Query genome.

Step 1. Genome alignment by minimap2	
Step 2. Unique anchor/alignment detection
Step 3. SV calling and Filtering by illumina reads mapping 
Step 4(Optional SV_PacBio.sh ): Combining indels called by pacbio read mapping

This pipeline calls SVs based on the anchors identified by Assemblytics (https://pubmed.ncbi.nlm.nih.gov/27318204/).

####################################################################
Dependencies:
   1. minimap2 (v2.11 or higher, https://github.com/lh3/minimap2)
   2. sam2delta.py (A python script from RaGOO package, https://github.com/malonge/RaGOO/raw/master/sam2delta.py)
   3. Assemblytics (Assemblytics, https://github.com/MariaNattestad/Assemblytics)
   4. samtools (v1.5 or higher, http://www.htslib.org/)
   5. blast+ 
   6. Python 2.7

####################################################################
# Genome file format requirement:
    1. fasta
    2. Chromosome name should be look like xxxch01, xxxch02 et al. 
    3. The contigs not anchored into pseudochromosome can be merged as xxxch00.
       For example, SL4.0ch00, SL4.0ch01, SL4.0ch02...SL4.0ch12 in Heinz 1706 genome;
                    PIMPch00, SPIMPch01, SPIMPch02...SPIMPch12 in LA2093 genome.
    4. SV calling will run on normal chromosomes first and then ch00.

####################################################################
Because the genome assembly may be incomplete, some unassembled regions might be identified as deletion. To avoid this mistake, this pipeline adopt illumina reads to filter detected SVs.
To do this, please align the illumina reads to the two reference genomes, respectively, by any read mapping program (e.g. bwa), and named the bam file as following.

Inputs:
     Reference_genome_file
     Query_genome_file
     QryRead2Ref.bam         Query illumina reads on Reference genome
     RefRead2Qry.bam         Reference illumina reads on Query genome
     Ref_self.bam            Reference illumina reads on Reference genome
     Qry_self.bam            Query illumina reads on Query genome

Command:      
      bash path_to_SV_calling_script/SV_calling.sh path_to_SV_calling_script <Reference_genome_file>   <Query_genome_file>  <Prefix_for_outputs> <number of threads>  <min_SV_size>   <max_SV_size>  <assemblytics_path>

For example:
      bash path_to_SV_calling_script/SV_calling.sh path_to_SV_calling_script SL4.0.genome.fasta Pimp_v1.4.fasta SP2SL 24 10 1000000 path_to_assemblytics_scripts
Final result for this example in 2 formats:
      SP2SL.Genome_comp_SV.tsv
      SP2SL.NR.bed (Same SVs as SP2SL.Genome_comp_SV.tsv, using Assemblytics output format. This one is used as input file in SV_PacBio.sh)


Note:
      It's a little tricky to set a proper max SV size. When we called a large deletion, we won't try to call any small indels in the "deleted" region. If this large deletion is wrong, we may miss real indels in that region. 
      For Solanum pimpinellifolium project, I firstly run this pipeline allowing up to 5Mb indel, and then manually check the detected SVs > 1Mb. Then, I rerun the pipeline with max_SV_size = 1Mb and combined the resulted indels and the checked ones > 1Mb.


####################################################################
If you have PacBio reads for the two genomes, you may call SVs using these reads, and combined the resulted SVs into the above genome comparison SVs.
Pleae put the vcf files based the two genomes under working directory. A sample vcf file is provided (Example.vcf).

Combine PacBio SVs:
      bash path_to_SV_calling_script/SV_PacBio.sh path_to_SV_calling_script <Reference_genome_file>   <Query_genome_file>  <Prefix_for_outputs> <number of threads> <Ref_base_pbsv_vcf> <Qry_base_pbsv_vcf>

For example:
      bash path_to_SV_calling_script/SV_PacBio.sh path_to_SV_calling_script SL4.0.genome.fasta  Pimp_v1.4.fasta   SP2SL  24  PimpReads2SL4.0.var.vcf  HeinzReads2Pimp.var.vcf

Final result for this example:
      SP2SL.Master_list.tsv



