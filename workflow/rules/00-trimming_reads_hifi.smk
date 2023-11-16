### Parse and prepare HiFi data for phased assembly 
rule hifi_adapter_filt:
    input:
        bam = path_reads_prefix + "/{sample}_HIFI_None.bam"
    output:
        fastq_clean = path_reads_prefix + "/{sample}_HIFI_None_reads.filt.fastq.gz",
        blast_out = path_reads_prefix + "/{sample}_HIFI_None_reads.contaminant.blastout",
        blocklist = path_reads_prefix + "/{sample}_HIFI_None_reads.blocklist",
        stats = path_reads_prefix + "/{sample}_HIFI_None_reads.stats"
    params:
        path_reads_prefix  = path_reads_prefix 
        prefix = "{sample}_HIFI_None_reads",
        min_len = 44,  # default value
        min_match = 97,  # default value
        threads = 8,  # default value
        outdir = path_reads_prefix
    conda:
        "../envs/hifiadapterfilt.yaml"  # Ensure you have all dependencies in this environment
    shell:
        """
        currpwd = echo"pwd" && \ 
        cd {input.path_reads_prefix}
        bash pbadapterfilt.sh -p {params.prefix} -l {params.min_len} -m {params.min_match} \
            -t {params.threads} -o {params.outdir} {input.bam} && \ 
        cd $currpwd
        """

# Define rule for running LongQC on HIFI trimmed reads before phased assembly (if they exist)
rule longqc_on_hifi_cleaned_reads:
    input:
        fastq = path_reads_prefix + "/{sample}_HIFI_None_reads.filt.fastq.gz"
    output:
        longqc_out = directory(path_out_prefix + "/00-LONGQC/{sample}_HIFI/")
    params:
        longqc_path = "/tarafs/data/home/qandres/catfish/01-PROGRAMS/LongQC-1.2.0c/"
    conda:
        "../envs/quality_control_reads.yaml"
    shell:
        """
        python {params.longqc_path}longQC.py sampleqc -x pb-hifi -o {output.longqc_out} {input.fastq}
        """

