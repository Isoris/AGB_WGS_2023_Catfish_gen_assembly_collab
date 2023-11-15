### Parse and prepare ONT nanopore data
# Define rule for running nanoQC on ONT reads (if they exist)
rule nanoqc_on_nanopore_reads:
    input: 
        fastq = path_reads_prefix + "/{sample}_NANOPORE_None_rawreads.fastq.gz"
    output:
        nanoqc_out = directory(path_out_prefix + "/00-NANOQC/{sample}_NANOPORE/")
    conda:
        "../envs/quality_control_reads.yaml"
    shell: 
        """
        nanoQC -o {output.nanoqc_out} {input.fastq} --minlen 100
        """

rule filtlong_on_nanopore_reads:
    input:
        fastq = path_reads_prefix + "/{sample}_NANOPORE_None_rawreads.fastq.gz"
    output:
        trimmed_ont_reads = path_reads_prefix + "/{sample}_NANOPORE_None_trimmed_reads.fastq.gz"
    params:
        min_length=1000
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        filtlong --min_length {params.min_length}  {input.fastq} | gzip > {output.trimmed_ont_reads}
        """