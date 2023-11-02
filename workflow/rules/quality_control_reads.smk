# snakemake --use-conda --conda-frontend mamba ...

rule bam_to_fastq:
    input:
        lambda wildcards: config["samples"][f"{wildcards.species}_{wildcards.sex}_{wildcards.method}_{wildcards.orientation}"][".bam"]
    output:
        fastq={path_reads_prefix}/config["samples"][f"{wildcards.species}_{wildcards.sex}_{wildcards.method}_{wildcards.orientation}"]_rawreads.fastq"
    conda:
        "quality_control_env.yaml"  # This YAML file contains the dependencies
    shell:
        "bam2fastq -c 1 -o {output.fastq} {input}"


rule mash_sketch:
    input:
        reads=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz"
    output:
        sketch=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz.msh"
    shell:
        "mash sketch -m 2 -o {output.sketch} {input.reads}"

rule mash_dist:
    input:
        ref_sketch="path/to/combined.msh",
        reads_sketch=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz.msh"
    output:
        distances="results/{{sample}}_distances.tab"
    shell:
        "mash dist -p 8 {input.ref_sketch} {input.reads_sketch} > {output.distances}"

rule mash_screen:
    input:
        ref_sketch="path/to/combined.msh",
        reads=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz"
    output:
        screen="results/{{sample}}_screen.tab"
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"


# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

rule run_jellyfish_count:
    input:
        reads=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz"
    output:
        directory=f"{path_out_prefix}/{{sample}}_rawreads_jellyfishcount.jq"
    params:
        jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish", 
        kmer_size = 21
    conda:
        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        source activate {conda}  # Activating the conda environment
        cd {params.longqc_path}  # Navigating to LongQC directory
        zcat _rawreads.fastq.gz | jellyfish count jellyfish count \
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
         raw_read =config["sample"]+"_mash_sketch.msh", 
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
                
