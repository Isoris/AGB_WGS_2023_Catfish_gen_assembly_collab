
# Define variables
# Define rule for running FastQC on raw reads
rule fastqc_on_illumina_raw_reads:
    input:
        read = path_reads_prefix + "{species}_{sex}_{method}_{orientation}.fq.gz",
    output:
        fastqc_out = path_out_prefix + "/00-FASTQC/{species}_{sex}_{method}_{orientation}/"
    conda:
        "rules/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        mkdir {output.fastqc_out} && \
        fastqc {input.read} -t {threads} -o {output.fastqc_out}
        """

rule adapter_removal_on_illumina:
    input:
        read1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD.fq.gz",
        read2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV.fq.gz"
    output:
        out1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD_trimmed.fq.gz",
        out2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV_trimmed.fq.gz"
    conda:
        "rules/quality_control_reads.yaml"
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
        read = path_reads_prefix + "/{species}_{sex}_{method}_{orientation}_trimmed.fq.gz",
    output:
        fastqc_out= path_out_prefix + "/00-FASTQC/{species}_{sex}_{method}_{orientation}_trimmed/"
    conda:
        "rules/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        mkdir -p {output.fastqc_out} && \
        fastqc {input.read} -t {threads} -o {output.fastqc_out}
        """        
