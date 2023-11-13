
### Parse and prepare HiC data
#Define rule for running FastQC on HiC reads (if they exist)
rule fastqc_on_hic_raw_reads:
    input:
        read_FWD = path_reads_prefix + "/{species}_{sex}_HIC_FWD.fq.gz",
        read_REV = path_reads_prefix + "/{species}_{sex}_HIC_REV.fq.gz"
    output:
        fastqc_out_FWD = directory(path_out_prefix + "/00-FASTQC/{species}_{sex}_HIC_FWD/"),
        fastqc_out_REV = directory(path_out_prefix + "/00-FASTQC/{species}_{sex}_HIC_REV/")
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        fastqc {input.read_FWD} -t {threads} -o {output.fastqc_out_FWD} && \
        fastqc {input.read_REV} -t {threads} -o {output.fastqc_out_REV}
        """

# Define rule for running AdapterRemoval on raw HiC  reads to trim them and remove adapters
rule adapter_removal_on_hic_raw_reads:
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

# Define rule for running FastQC on trimmed HiC reads (if they exist)
rule fastqc_on_hic_trimmed_reads:
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
