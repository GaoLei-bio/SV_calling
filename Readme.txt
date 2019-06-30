Command:
    bash SV_pipeline.sh  <Reference_genome_file>   <Query_genome_file>  <Prefix_for_outputs> <number of threads>


For example: 
    In current directory, We have:
            a reference genome:  SL4.0.genome.fasta
            a query genome:      Pimp_v1.2.fasta

    The prefix for outputs: SP2SL
    Number of CUPs we can use: 32

The command is:
    bash SV_pipeline.sh  SL4.0.genome.fasta  Pimp_v1.2.fasta   SP2SL   32 > log.txt

