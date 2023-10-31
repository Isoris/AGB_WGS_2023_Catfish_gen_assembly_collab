# Define variables to be called in the following rules
{output_folder_prefix} = ~/scratch/catfish/03-OUTPUTS/



# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

rule run_jellyfish_count:
    input:
        "{sample}.fastq"
    output:
        directory("path/to/output/{sample}_LongQC_output/")
    params:
        jellyfish_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1"  # adjust this to your jellyfish installation path
    conda:
        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        source activate {conda}  # Activating the conda environment
        cd {params.longqc_path}  # Navigating to LongQC directory
        jellyfish count \
        -C {params.preset} \
        -m {params.nproc} \
        -s 1000000000 \
        -o {output} \
        -t {threads} \ 
        {input}.fastq  
        """
        
rule run_jellyfish_histo:
    input:
        "{sample}.jq"
    output:
        directory("path/to/output/{sample}_LongQC_output/")
    params:
        jellyfish_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1"  # adjust this to your jellyfish installation path
    conda:
        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        source activate {conda}  # Activating the conda environment
        cd {params.longqc_path}  # Navigating to jellyfish output directory
        jellyfish histo -t {threads} reads.jf > reads.histo 
        """        

# Role: 
rule run_genomescope:
    input:
        "{sample}.jq"
    output:
        directory("path/to/output/{sample}_LongQC_output/")
    params:
        jellyfish_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1"  # adjust this to your jellyfish installation path
    conda:
        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        source activate {conda}  # Activating the conda environment
        cd {jellyfish_output_path}  # Navigating to jellyfish output directory
         Rscript genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose] 
        """            

# Define rule for running LongQC on raw reads (longCQ is under catfish/01-PROGRAMS/) should run within mamba quality_control environment 
# Role: For the QC of long reads from HiFi

rule run_LongQC:
    input:
        "path/to/input/{sample}.fastq"
    output:
        directory("path/to/output/{sample}_LongQC_output/")
    params:
        preset="pb-sequel",  # specify preset, change accordingly
        nproc="$(nproc)",  # number of processors, can be adjusted
        sample_name="{sample}",  # sample name
        longqc_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1.2.0c"  # adjust this to your LongQC installation path
    conda:
        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        source activate {conda}  # Activating the conda environment
        cd {params.longqc_path}  # Navigating to LongQC directory
        python longQC.py sampleqc \
        -x {params.preset} \
        -p {params.nproc} \
        -o {output} \
        {input}
        """

# Define rule for running nanoQC on nanopore raw data
#conda install -c bioconda nanoqc

rule run_nanoqc:
    input: 
         raw_read = 
    output: 
         path_output_nanoqc = "{output_folder_prefix}/nanoqc"
    shell: 
         """
      mkdir -p {output_folder_prefix}/nanoqc
      nanoQC -o {output} {input}.fastq --minlen 100
         """

# Define rule for running FastQC on raw reads
rule fastqc_raw_reads:
    input:
        raw_read = lambda wildcards: config["samples"][wildcards.sample][wildcards.read_type]
    output:
        report = "fastqc/{sample}_{read_type}_fastqc.zip"
    params:
        output_dir = "fastqc/",
        nthreads = config["CPUS"]
    run:
        output_path = os.path.join(params.output_dir, os.path.basename(output.report))
        if not os.path.exists(output_path):
            shell("fastqc -o {params.output_dir} -t {params.nthreads} {input.raw_read}")

# Define rule for running FastQC on trimmed reads (if they exist)
rule fastqc_trimmed_reads:
    input:
        trimmed_read = lambda wildcards: config["trimmed_samples"][wildcards.sample][wildcards.read_type]
    output:
        report = "fastqc/{sample}_{read_type}_trimmed_fastqc.zip"
    params:
        output_dir = "fastqc/",
        nthreads = config["CPUS"]
    run:
        output_path = os.path.join(params.output_dir, os.path.basename(output.report))
        if os.path.exists(input.trimmed_read):
            if not os.path.exists(output_path):
                shell("fastqc -o {params.output_dir} -t {params.nthreads} {input.trimmed_read}")
                
rule quast_quality_control:
    input:
        spades_assembly = "{path_prefix_step2_spades}/output_spades_{species}_{sex}_{version}/XXX_spades.fasta",
        ref_genome = "{path_prefix_ref_genomes}/{species}_ref.fasta",
        gff = "{path_prefix_ref_genomes}/{species}.gff",
        pe1 = "{path_prefix_reads_trimmed}/{sample}_reads_trim1.fq",
        pe2 = "{path_prefix_reads_trimmed}/{sample}_reads_trim2.fq"
    output:
        report = "Output_QUAST_{sample}_{version_num}"
    threads: config['CPUS']
    shell:
        """
        echo "Start of QUAST quality control of the Spades de-novo assemblies for the illumina data"
        echo "Run QUAST with Busco and silva and gridss and rna+gene finding, report all metrics, create circos"
        echo "QUAST Quality control step starting for {wildcards.sex} {wildcards.species}"

        {config[progQUAST]} \
            {input.spades_assembly} \
            -r {input.ref_genome} \
            -g {input.gff} \
            --pe1 {input.pe1} \
            --pe2 {input.pe2} \
            --circos \
            --threads {threads} \
            --gene-finding \
            --rna-finding \
            --conserved-genes-finding \
            --report-all-metrics \
            -L \
            -o {output.report}

        echo "QUAST Quality control step finished for {wildcards.species} {wildcards.sex} {wildcards.version}"
        """
