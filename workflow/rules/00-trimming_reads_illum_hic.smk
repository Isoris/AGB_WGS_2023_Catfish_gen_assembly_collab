
# Define variables
### Parse and prepare Illumina data
# Define rule for running FastQC on raw illumina reads
rule fastqc_on_illumina_raw_reads:
    input:
        read_FWD = path_reads_prefix + "/{sample}_ILLUMINA_FWD.fq.gz",
        read_REV = path_reads_prefix + "/{sample}_ILLUMINA_REV.fq.gz"
    output:
        fastqc_out_FWD_html = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_FWD/{sample}_ILLUMINA_FWD_fastqc.html",
        fastqc_out_FWD_zip = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_FWD/{sample}_ILLUMINA_FWD_fastqc.zip",
        fastqc_out_REV_html = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_REV/{sample}_ILLUMINA_REV_fastqc.html",
        fastqc_out_REV_zip = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_REV/{sample}_ILLUMINA_REV_fastqc.zip"
    params:
        path_out_prefix = path_out_prefix
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        fastqc {input.read_FWD} -t {threads} -o {params.path_out_prefix}/00-FASTQC/{wildcards.sample}_ILLUMINA_FWD/ && \
        fastqc {input.read_REV} -t {threads} -o {params.path_out_prefix}/00-FASTQC/{wildcards.sample}_ILLUMINA_REV/
        """

# Define rule for running AdapterRemoval on raw illumina reads to trim them and remove adapters
rule adapter_removal_on_illumina_raw_reads:
    input:
        read1 = path_reads_prefix + "/{species}_{sex}_{method}_FWD.fq.gz",
        read2 = path_reads_prefix + "/{species}_{sex}_{method}_REV.fq.gz"
    output:
        out1 = path_reads_prefix + "/{species}_{sex}_{method}_FWD_trimmed_reads.fastq.gz",
        out2 = path_reads_prefix + "/{species}_{sex}_{method}_REV_trimmed_reads.fastq.gz"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        AdapterRemoval \
            --file1 {input.read1} \
            --file2 {input.read2} \
            --gzip \
            --output1 {output.out1} \
            --output2 {output.out2} \
            --threads {threads} 
        """

rule fastqc_on_illumina_trimmed_reads:
    input:
        read_FWD = path_reads_prefix + "/{sample}_ILLUMINA_FWD_trimmed_reads.fastq.gz",
        read_REV = path_reads_prefix + "/{sample}_ILLUMINA_REV_trimmed_reads.fastq.gz"
    output:
        fastqc_out_FWD_html = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_FWD_trim/{sample}_ILLUMINA_FWD_trimmed_reads_fastqc.html",
        fastqc_out_FWD_zip = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_FWD_trim/{sample}_ILLUMINA_FWD_trimmed_reads_fastqc.zip",
        fastqc_out_REV_html = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_REV_trim/{sample}_ILLUMINA_REV_trimmed_reads_fastqc.html",
        fastqc_out_REV_zip = "{path_out_prefix}/00-FASTQC/{sample}_ILLUMINA_REV_trim/{sample}_ILLUMINA_REV_trimmed_reads_fastqc.zip"
    params:
        path_out_prefix = path_out_prefix
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        fastqc {input.read_FWD} -t {threads} -o {params.path_out_prefix}/00-FASTQC/{wildcards.sample}_ILLUMINA_FWD_trim/ && \
        fastqc {input.read_REV} -t {threads} -o {params.path_out_prefix}/00-FASTQC/{wildcards.sample}_ILLUMINA_REV_trim/
        """

