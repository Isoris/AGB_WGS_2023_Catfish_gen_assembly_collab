
# Define variables
### Parse and prepare Illumina data
# Define rule for running FastQC on raw illumina reads
rule fastqc_on_illumina_raw_reads:
    input:
        read_FWD = path_reads_prefix + "/{species}_{sex}_ILLUMINA_FWD.fq.gz",
        read_REV = path_reads_prefix + "/{species}_{sex}_ILLUMINA_REV.fq.gz"
    output:
        fastqc_out_FWD = directory(path_out_prefix + "/00-FASTQC/{species}_{sex}_ILLUMINA_FWD/"),
        fastqc_out_REV = directory(path_out_prefix + "/00-FASTQC/{species}_{sex}_ILLUMINA_REV/")
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        fastqc {input.read_FWD} -t {threads} -o {output.fastqc_out_FWD} && \
        fastqc {input.read_REV} -t {threads} -o {output.fastqc_out_REV}
        """

# Define rule for running AdapterRemoval on raw illumina reads to trim them and remove adapters
rule adapter_removal_on_illumina_raw_reads:
    input:
        read1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD.fq.gz",
        read2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV.fq.gz"
    output:
        out1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD_trimmed_reads.fastq.gz",
        out2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV_trimmed_reads.fastq.gz"
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

# Define rule for running FastQC on trimmed reads (if they exist)
rule fastqc_on_illumina_trimmed_reads:
    input:
        read = path_reads_prefix + "/{species}_{sex}_{method}_{orientation_pe}_trimmed_reads.fastq.gz"
    output:
        fastqc_out = directory(path_out_prefix + "/00-FASTQC/{species}_{sex}_{method}_{orientation_pe}_trimmed/")
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        fastqc {input.read} -t {threads} -o {output.fastqc_out}
        """
        
### Parse and prepare HiFi data
rule bam_to_fastq_on_hifi_raw_reads:
    input:
        bam = path_reads_prefix + "/{species}_{sex}_HIFI_None.bam"
    output:
        fastq = path_reads_prefix + "/{species}_{sex}_HIFI_None_rawreads.fastq"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        bam2fastq -o {output.fastq} {input.bam} 
        """

rule seqtk_trim20bases_hifi_adapters_on_raw_reads:
    input:
        raw_read_hifi = path_reads_prefix + "/{species}_{sex}_HIFI_None_rawreads.fastq"
    output:
        trimmed_hifi_reads = path_reads_prefix + "/{species}_{sex}_HIFI_None_reads.fastq"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """ 
        seqtk trimfq -b 20 -e 20 {input.raw_read_hifi} > {output.trimmed_hifi_reads} 
        """

rule gzip_hifi_trimmed_reads:
    input:
        trimmed_hifi_reads = path_reads_prefix + "/{species}_{sex}_HIFI_None_reads.fastq"
    output:
        fastq = path_reads_prefix + "/{species}_{sex}_HIFI_None_trimmed_reads.fastq.gz"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """ 
        gzip -c -5 {input.trimmed_hifi_reads} > {output.fastq}
        """


# Define rule for running LongQC on HIFI trimmed reads (if they exist)

rule longqc_on_hifi_trimmed_reads:
    input:
        fastq = path_reads_prefix + "/{sample}_HIFI_None_trimmed_reads.fastq.gz"
    output:
        longqc_out = directory(path_out_prefix + "/00-LONGQC/{sample}_HIFI/")
    params:
        longqc_path = "/tarafs/data/home/qandres/catfish/01-PROGRAMS/LongQC-1.2.0c/"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        python {params.longqc_path}longQC.py sampleqc -x pb-hifi -o {output.longqc_out} {input.fastq}.gz
        """

### Parse and prepare ONT nanopore data
# Define rule for running nanoQC on ONT reads (if they exist)
rule nanoqc_on_nanopore_reads:
    input: 
        fastq = path_reads_prefix + "/{sample}_NANOPORE_rawreads.fastq.gz"
    output:
        nanoqc_out = directory(path_out_prefix + "/00-NANOQC/{sample}_NANOPORE/")
    conda:
        "../envs/quality_control_reads.yaml"
    shell: 
        """
        nanoQC -o {output.nanoqc_out} {input.fastq} --minlen 100
        """
