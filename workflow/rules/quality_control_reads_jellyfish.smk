
### Jellyfish
# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

rule run_jellyfish_count:
    input:
        reads=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz"
    output:
        jellyfish_count=f"{path_out_prefix}/{{sample}}_jellyfishcount.jf"
    params:
        jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish",
        kmer_size = 21
    conda:
        "quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        zcat {input.reads} | {params.jellyfish_path} count -m {params.kmer_size} -s 100M -t {threads} \
        -C -o {output.jellyfish_count} /dev/stdin
        """

rule run_jellyfish_histo:
    input:
        jellyfish_count=f"{path_out_prefix}/{{sample}}_jellyfishcount.jf"
    output:
        histo=f"{path_out_prefix}/{{sample}}_jellyfish_histo.histo"
    params:
         jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish"
    conda:
        "quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        "{params.jellyfish_path} histo -t {threads} {input.jellyfish_count} > {output.histo}"


### GenomeScope
# Role: Genomescope is used for genome profiling based on k-mer analysis
rule run_genomescope:
    input:
        jellyfish_histo=f"{path_out_prefix}/{{sample}}_jellyfish_histo.histo"
    output:
        directory(f"{path_out_prefix}/{{sample}}_genomescope_output/")
    params:
        genomescope_script="scripts/genomescope.R",  # adjust this to your genomescope script path
        kmer_length=21,  # adjust as needed
        read_length=100,  # adjust based on your read length
        kmer_max=None,  # optional, set a value if required
        verbose=False  # set to True if verbose output is needed
    conda:
        "quality_control_reads.yaml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        Rscript {params.genomescope_script} {input.jellyfish_histo} {params.kmer_length} {params.read_length} {output} {params.kmer_max if params.kmer_max else ''} {'--verbose' if params.verbose else ''}
        """

# Role: 
#rule run_genomescope:
#    input:
#        "{sample}.jq"
#    output:
#        directory("path/to/output/{sample}_LongQC_output/")
#    params:
#        jellyfish_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1"  # adjust this to your jellyfish installation path
#    conda:
#        "/quality_control_reads.yaml"  # Replace with the path to your 'quality_control' environment file
#    shell:
#        """
#        cd {jellyfish_output_path}  # Navigating to jellyfish output directory
#         Rscript genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose] 
#        """            

# Define rule for running LongQC on raw reads (longCQ is under catfish/01-PROGRAMS/) should run within mamba quality_control environment 
# Role: For the QC of long reads from HiFi

#rule run_LongQC:
#    input:
#        "path/to/input/{sample}.fastq"
#    output:
#        directory("path/to/output/{sample}_LongQC_output/")
#    params:
#        preset="pb-sequel",  # specify preset, change accordingly
#        nproc="$(nproc)",  # number of processors, can be adjusted
#        sample_name="{sample}",  # sample name
#        longqc_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1.2.0c"  # adjust this to your LongQC installation path
#    conda:
#        "/quality_control.yml"  # Replace with the path to your 'quality_control' environment file
#    shell:
#        """
#        source activate {conda}  # Activating the conda environment
#        cd {params.longqc_path}  # Navigating to LongQC directory
#       python longQC.py sampleqc \
#        -x {params.preset} \
#        -p {params.nproc} \
#        -o {output} \
#        {input}
#        """

# Define rule for running nanoQC on nanopore raw data
#conda install -c bioconda nanoqc

#rule run_nanoqc:
#    input: 
#        raw_read =config["sample"]+"_mash_sketch.msh", 
#   output:
#         path_output_nanoqc = "{output_folder_prefix}/nanoqc"
#    shell: 
#         """
#      mkdir -p {output_folder_prefix}/nanoqc
#      nanoQC -o {output} {input}.fastq --minlen 100
#         """
#
# Define rule for running FastQC on raw reads
#rule fastqc_raw_reads:
#    input:
#        raw_read = lambda wildcards: config["samples"][wildcards.sample][wildcards.read_type]
#    output:
#        report = "fastqc/{sample}_{read_type}_fastqc.zip"
#    params:
#       output_dir = "fastqc/",
#        nthreads = config["CPUS"]
#    run:
#        output_path = os.path.join(params.output_dir, os.path.basename(output.report))
#        if not os.path.exists(output_path):
#            shell("fastqc -o {params.output_dir} -t {params.nthreads} {input.raw_read}")
#
# Define rule for running FastQC on trimmed reads (if they exist)
#rule fastqc_trimmed_reads:
#    input:
#        trimmed_read = lambda wildcards: config["trimmed_samples"][wildcards.sample][wildcards.read_type]
#    output:
#        report = "fastqc/{sample}_{read_type}_trimmed_fastqc.zip"
#    params:
#        output_dir = "fastqc/",
#        nthreads = config["CPUS"]
#    run:
#        output_path = os.path.join(params.output_dir, os.path.basename(output.report))
#        if os.path.exists(input.trimmed_read):
#            if not os.path.exists(output_path):
#                shell("fastqc -o {params.output_dir} -t {params.nthreads} {input.trimmed_read}")

      
                
# snakemake --use-conda --conda-frontend mamba ...
#CG_M_ILLUMINA_PE_FWD.fastq
#CG_M_NANOPORE.fastq
#CG_M_HIFI.bam

#rule example_rule:
#    input:
#        reads=lambda wildcards: f"{wildcards.sample}_{config['method']}_reads.fastq.gz" if config['method'] == 'Illumina' else None
#    output:
#        processed_reads="{sample}_Illumina_processed.fastq.gz"
#    shell:
#        """
#        # Command to process Illumina reads
#       process_illumina_reads -i {input.reads} -o {output.processed_reads}
#        """
#rule example_rule:
#    input:
#        reads="{sample}_{method}_reads.fastq.gz"
#    output:
#        processed_reads="{sample}_{method}_processed.fastq.gz"
#    shell:
#        """
#        if [ "{wildcards.method}" == "Illumina" ]; then
#            # Command for processing Illumina reads
 #           process_illumina_reads -i {input.reads} -o {output.processed_reads}
 #       elif [ "{wildcards.method}" == "Nanopore" ]; then
            # Command for processing Nanopore reads
 #           process_nanopore_reads -i {input.reads} -o {output.processed_reads}
        # ... handle other methods as needed
  #      fi
   #     """