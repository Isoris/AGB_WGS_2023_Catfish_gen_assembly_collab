
# Define variables
### Parse and prepare Illumina data
# Define rule for running FastQC on raw illumina reads
rule fastqc_on_illumina_raw_reads:
    input:
        read_FWD = path_reads_prefix + "/{species}_{sex}_ILLUMINA_FWD.fq.gz",
        read_REV = path_reads_prefix + "/{species}_{sex}_ILLUMINA_REV.fq.gz"
    output:
        fastqc_out_FWD = path_out_prefix + "/00-FASTQC/{species}_{sex}_ILLUMINA_FWD/",
        fastqc_out_REV = path_out_prefix + "/00-FASTQC/{species}_{sex}_ILLUMINA_REV/"
    conda:
        "envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        mkdir -p {output.fastqc_out_FWD} && \
        mkdir -p {output.fastqc_out_REV} && \
        fastqc {input.read_FWD} -t {threads} -o {output.fastqc_out_FWD} && \
        fastqc {input.read_REV} -t {threads} -o {output.fastqc_out_REV}
        """

# Define rule for running AdapterRemoval on raw illumina reads to trim them and remove adapters
rule adapter_removal_on_illumina:
    input:
        read1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD.fq.gz",
        read2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV.fq.gz"
    output:
        out1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD_trimmed_reads.fastq.gz",
        out2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV_trimmed_reads.fastq.gz"
    conda:
        "envs/quality_control_reads.yaml"
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

# Define rule for running FastQC on trimmed reads (if they exist)
rule fastqc_on_illumina_trimmed_reads:
    input:
        read = path_reads_prefix + "/{species}_{sex}_{method}_{orientation_pe}_trimmed_reads.fastq.gz"
    output:
        fastqc_out = path_out_prefix + "/00-FASTQC/{species}_{sex}_{method}_{orientation_pe}_trimmed/"
    conda:
        "envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        mkdir -p {output.fastqc_out} && \
        fastqc {input.read} -t {threads} -o {output.fastqc_out}
        """        


### Parse and prepare HiFi data
rule bam_to_fastq:
    input:
        bam = path_reads_prefix + "/{species}_{sex}_HIFI_None.bam"
    output:
        fastq = path_reads_prefix + "/{species}_{sex}_HIFI_None_rawreads.fastq"
    conda:
        "envs/quality_control_reads.yaml"
    shell:
        """
        bam2fastq -o {output.fastq} {input.bam} 
        """

rule remove_hifi_adapters:
    input:
        raw_read_hifi = path_reads_prefix + "/{species}_{sex}_HIFI_None_rawreads.fastq"
    output:
        trimmed_hifi_reads = path_reads_prefix + "/{species}_{sex}_HIFI_None_reads.fastq"
    conda:
        "envs/quality_control_reads.yaml"
    shell:
        """ 
        seqtk trimfq -b 20 -e 20 {input.raw_read_hifi} > {output.trimmed_hifi_reads} 
        """

rule gzip_hifi_fastq:
    input:
        trimmed_hifi_reads = path_reads_prefix + "/{species}_{sex}_HIFI_None_reads.fastq"
    output:
        fastq = path_reads_prefix + "/{species}_{sex}_HIFI_None_reads.fastq.gz"
    conda:
        "envs/quality_control_reads.yaml"
    shell:
        """ 
        gzip -5 {input.trimmed_hifi_reads} && \
        mv {input.trimmed_hifi_reads}.gz {output.fastq}
        """

### Parse and prepare ONT nanopore data








### Parse and prepare HiC data


# Define rule for running FastQC on HiC reads (if they exist)
#rule fastqc_on_hic:
#    input:
#        read = path_reads_prefix + "/{species}_{sex}_HIFI_None.fq.gz"
#    output:
#        fastqc_out = path_out_prefix + "/00-FASTQC/{species}_{sex}_HIFI_None/"
#    conda:
#        "envs/quality_control_reads.yaml"
#    shell:
#        """
#        mkdir -p {output.fastqc_out} && \
#        fastqc {input.read} -t {threads} -o {output.fastqc_out}
#        """
