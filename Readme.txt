####################################################################
SV calling scripts
By Lei Gao (lg397@cornell.edu or leigao@wbgcas.cn)
Fei lab (http://bioinfo.bti.cornell.edu/)
Version 1.1
Date: 07/04/2020
This pipeline is designed to detect structural variations (SVs) between high-quality genomes of 2 closely related species/variaties.
One is called Reference genome, the other one is Query genome.

SV_calling.sh:
Step 1. Genome alignment by minimap2	
Step 2. Unique anchor/alignment detection
Step 3. SV calling based on genome comparison, and then Filtering by illumina reads mapping 

SV_PacBio.sh (Optional): 
Step 4: Combining indels called by pacbio read mapping

####################################################################
How to run:
Please (1) copy all scripts to path_to_SV_calling_script (absolute path);
       (2) put all inputs in current directory;
       (3) install all dependencies

Step 1-3:
  SV_calling.sh
  Command:      
      bash path_to_SV_calling_script/SV_calling.sh \
           path_to_SV_calling_script \
           <Reference_genome_file>   \
           <Query_genome_file>  \
           <Prefix_for_outputs> \
           <number of threads>  \
           <min_SV_size>   \
           <max_SV_size>  \
           <assemblytics_path>
  Inputs:
     Reference_genome_file   Fasta format, see "Genome file format requirement" for detail
     Query_genome_file       Fasta format, see "Genome file format requirement" for detail
     QryRead2Ref.bam         Sorted bam file, Query illumina reads on Reference genome
     RefRead2Qry.bam         Sorted bam file, Reference illumina reads on Query genome
     Ref_self.bam            Sorted bam file, Reference illumina reads on Reference genome
     Qry_self.bam            Sorted bam file, Query illumina reads on Query genome

  For example:
      bash path_to_SV_calling_script/SV_calling.sh \
           path_to_SV_calling_script \
           SL4.0.genome.fasta \
           Pimp_v1.4.fasta \
           SP2SL 24 10 1000000 \
           path_to_assemblytics_scripts

  Final result for this example (Same SVs in 2 formats):
      SP2SL.Genome_comp_SV.tsv
      SP2SL.NR.bed (Assemblytics output format. This one is input file for next step)

======================
Step 4:
  SV_PacBio.sh
  Command:
      bash path_to_SV_calling_script/SV_PacBio.sh \
           path_to_SV_calling_script \
           <Reference_genome_file>   \
           <Query_genome_file>  \
           <Prefix_for_outputs> \
           <number of threads> \
           <Ref_base_pbsv_vcf> \
           <Qry_base_pbsv_vcf>

  Inputs:
     Prefix_for_outputs.NR.bed      Result of SV_calling.sh
     Reference_genome_file          Fasta format
     Query_genome_file              Fasta format
     Ref_base_pbsv_vcf              vcf file. SV calling by Query sample PacBio read alignments on Reference genome
     Qry_base_pbsv_vcf              vcf file. SV calling by Reference sample PacBio read alignments on Query genome
 
For example:
      bash path_to_SV_calling_script/SV_PacBio.sh \
           path_to_SV_calling_script \
           SL4.0.genome.fasta  \
           Pimp_v1.4.fasta   \
           SP2SL  24  \
           PimpReads2SL4.0.var.vcf  \
           HeinzReads2Pimp.var.vcf

Final result for this example:
      SP2SL.Master_list.tsv

####################################################################
Dependencies:
   1. minimap2 (v2.11 or higher, https://github.com/lh3/minimap2)
   2. sam2delta.py (A python script from RaGOO package, https://github.com/malonge/RaGOO/raw/master/sam2delta.py)
   3. Assemblytics (Assemblytics, https://github.com/MariaNattestad/Assemblytics)
   4. samtools (v1.5 or higher, http://www.htslib.org/)
   5. blast+ 
   6. Python 2.7

####################################################################
Genome file format requirement:
    1. fasta
    2. Chromosome name should be look like xxxch01, xxxch02 et al. 
    3. The contigs not anchored into pseudochromosome can be merged as xxxch00.
       For example, SL4.0ch00, SL4.0ch01, SL4.0ch02...SL4.0ch12 in Heinz 1706 genome;
                    PIMPch00, SPIMPch01, SPIMPch02...SPIMPch12 in LA2093 genome.
    4. SV calling will run on normal chromosomes first and then ch00.

####################################################################
Note:
  1. SV_calling.sh
     1.1 This pipeline calls SVs based on the anchors identified by Assemblytics (https://pubmed.ncbi.nlm.nih.gov/27318204/).
     1.2 Because the genome assembly may be incomplete, some unassembled regions might be identified as deletion. To avoid this mistake, this pipeline adopt illumina reads to filter detected SVs.
         To do this, please align the illumina reads to the two reference genomes, respectively, by any read mapping program (e.g. bwa), and named the bam file as following.
     1.3 It's a little tricky to set a proper max SV size. When we called a large deletion, we won't try to call any small indels in the "deleted" region. If this large deletion is wrong, we may miss real indels in that region. 
         For Solanum pimpinellifolium project, I firstly run this pipeline allowing up to 5Mb indel, and then manually check the detected SVs > 1Mb. Then, I rerun the pipeline with max_SV_size = 1Mb and combined the resulted indels and the checked ones > 1Mb.

   2. SV_PacBio.sh
     If you have PacBio reads for the two genomes, you may call SVs using these reads, and combined the resulted SVs into the above genome comparison SVs. 
     Pleae put the vcf files based the two genomes under working directory. A sample vcf file is provided (Example.vcf).


