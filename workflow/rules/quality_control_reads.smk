#

rule bam_to_fastq:
    input:
        "path/to/{sample}.bam"
    output:
        fastq="path/to/{sample}_reads.fastq"
    shell:
        "bam2fastq -o {output.fastq} {input}"

rule mash_sketch:
    input:
        reads="path/to/{sample}_reads.fastq.gz"
    output:
        sketch="path/to/{sample}_reads.fastq.gz.msh"
    shell:
        "mash sketch -m 2 -o {output.sketch} {input.reads}"

rule mash_sketch_hifi:
    input:
        reads_hifi_bam=config["sample"]+".bam"
    output:
        mash_reads_hifi_sketches=config["sample"]+"_mash_sketch.msh",
    shell:
        "mash sketch -p {threads} {input} > {output}"

rule mash_dist:
    input:
        ref_sketch="path/to/combined.msh",
        reads_sketch="path/to/{sample}_reads.fastq.gz.msh"
    output:
        distances="results/{sample}_distances.tab"
    shell:
        "mash dist -p 8 {input.ref_sketch} {input.reads_sketch} > {output.distances}"

rule mash_screen:
    input:
        ref_sketch="path/to/combined.msh",
        reads="path/to/{sample}_reads.fastq.gz"
    output:
        screen="results/{sample}_screen.tab"
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"



# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

rule run_jellyfish_count:
    input:
        "{sample}.fastq"
    output:
        directory("path/to/output/{sample}_LongQC_output/")
    params:
        jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish", 
        kmer_size = 21
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
                
