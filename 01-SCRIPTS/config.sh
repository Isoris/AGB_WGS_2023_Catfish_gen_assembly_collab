# Output directories for each step for the genome assembly of illumina reads (denovo)
declare -A output=( ["fastqc"]="output_fastqc" ["adapter_removal"]="output_adapter_removal" ["fastqc_after_adapter_removal"]="output_fastqc_after_adapter_removal" ["spades"]="output_spades" ["quast"]="output_quast" ["mummer"]="output_mummer" ["progressive_cactus"]="output_progressive_cactus" )

# Steps in the workflow
step=("fastqc" "adapter_removal" "fastqc_after_adapter_removal" "spades" "quast" "mummer" "progressive_cactus")

# Number of CPUs to use
CPUS="set the CPU number here"

# Short names for samples or experiments
shortNames=("Sample1" "Sample2" "Sample3")
Sex=("female" "male")
Species=("gariepinus" "macrocephalus" "hybrid")
Version=("mapped_assembly" "denovo" "reference")
Version_num=("Ver_1" "Ver_2" "Ver_3")

# Path to reference genomes
path_prefix_ref_genomes="/path/to/reference/genomes/"  # Reference genomes location
ref_genome_clarias_macrocephalus="${path_prefix_ref_genomes}/clarias_macrocephalus_ref.fasta"
ref_genome_clarias_gariepinus="${path_prefix_ref_genomes}/clarias_gariepinus_ref.fasta"
ref_genome_f1="${path_prefix_ref_genomes}/clarias_f1_ref.fasta"  # If available

# Sample information
declare -A shortNames=( ["Sample1"]="female" ["Sample2"]="male" ["Sample3"]="female" )
declare -A Species=( ["Sample1"]="gariepinus" ["Sample2"]="macrocephalus" ["Sample3"]="hybrid" )
declare -A Version=( ["Sample1"]="mapped_assembly" ["Sample2"]="denovo" ["Sample3"]="reference" )
declare -A Version_num=( ["Sample1"]="Ver_1" ["Sample2"]="Ver_2" ["Sample3"]="Ver_3" )

# Raw and trimmed read filenames
declare -A reads1=( ["Sample1"]="read1_sample1.fq.gz" ["Sample2"]="read1_sample2.fq.gz" ["Sample3"]="read1_sample3.fq.gz" )
declare -A reads2=( ["Sample1"]="read2_sample1.fq.gz" ["Sample2"]="read2_sample2.fq.gz" ["Sample3"]="read2_sample3.fq.gz" )
declare -A reads_trim1=( ["Sample1"]="read1_trim_sample1.fq.gz" ["Sample2"]="read1_trim_sample2.fq.gz" ["Sample3"]="read1_trim_sample3.fq.gz" )
declare -A reads_trim2=( ["Sample1"]="read2_trim_sample1.fq.gz" ["Sample2"]="read2_trim_sample2.fq.gz" ["Sample3"]="read2_trim_sample3.fq.gz" )


# Prefix Paths
#	reads (raw)
path_prefix_reads="/path/to/raw/reads/"  # Raw read files location
#	reads (trimmed - step1)
path_prefix_step1_trimmed="/path/to/step1/"  # Output for step 1
#	assemblies (unpolished)
path_unpolished_assemblies="path/to/assemblies/before/illumina"
#	assemblies (polished - PILON and Freebayes)
path_polished_assemblies_pilon=
path_polished_assemblies_freebayes=
#	programs (general prefix)
path_prefix_programs="/path/to/programs/"  # Executable program files location


#	programs outputs (for each step) -- 							Denovo Illumina
path_prefix_step_0_fastqc_before_trimming="/path/to/step2/"  					# Output for step 0 fastqc
path_prefix_step_0_mash="/path/to/step2/"  							# Output for step 0 mash
path_prefix_step_1_adapterremoval="/path/to/step2/"  						# Output for step 1 adapter removal
path_prefix_step_0_fastqc_after_trimming="/path/to/step2/"  					# Output for step 0 fastqc after adapter removal
path_prefix_step_2_spades="/path/to/step2/"  							# Output for step 2 spades 
path_prefix_step_3_quast= 									# Output for step 3 quast for getting assembly metrics
path_prefix_step_4_mummer=									# Output for step 4 MuMmer+nucmer for one-to-one whole genome alignments 
path_prefix_step_5_progressivecactus=								# Output for step 5 Progressive cactus for whole genome alignments
#	programs outputs (for each step) -- 							Polishing Illumina
path_prefix_step_1_freebayes=									# Output for step 1 Freebayes for variant calling of short reads on cns.
path_prefix_step_1_pilon=									# Output for step 1 pilon for variant calling of short reads on consensus
#	programs outputs (for each step) -- 							Pacbio consensus
path_prefix_step_1_pacbio=									# Output for step 1 Consensus from Pacbio + nanopore assembly (denovo)




# Source the config file
source ./progpaths.sh

echo "Sourced progpaths.sh"
# End of config.sh

