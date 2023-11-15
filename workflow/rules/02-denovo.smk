# Define variables to be called within rules
# Define rules for de-novo assembly of reads

rule run_hifasm_hifi_UL_hic:
    input:
        hifi_reads = path_reads_prefix + "/{sample}_HIFI_None_trimmed_reads.fastq.gz",
        ont_reads = path_reads_prefix + "/{sample}_NANOPORE_None_trimmed_reads.fastq.gz",
        hic1_reads = path_reads_prefix + "/{sample}_HIC_FWD.fastq.gz",
        hic2_reads = path_reads_prefix + "/{sample}_HIC_REV.fastq.gz" # Assuming REV is the reverse read
    output:
        primary_ctg = path_reads_prefix + "/{sample}.p_ctg.gfa",
        alternate_ctg = path_reads_prefix + "/{sample}.a_ctg.gfa",
        hap1_ctg = path_reads_prefix + "/{sample}.hap1.p_ctg.gfa",
        hap2_ctg = path_reads_prefix + "/{sample}.hap2.p_ctg.gfa"
    params:
        hifiasm_path = path_prog_prefix + "/hifiasm-0.19.7/hifiasm",
        out_prefix = "{sample}"
    shell:
        """
        {params.hifiasm_path} \
        --dual-scaf \
        --primary \
        --h1 {input.hic1_reads} \
        --h2 {input.hic2_reads} \
        --ul {input.ont_reads} \
        {input.hifi_reads} \
        -o {params.out_prefix}
        """

rule gfa_to_fasta_after_hifasm:
    input:
        primary_ctg = path_reads_prefix + "/{sample}.p_ctg.gfa",
        alternate_ctg = path_reads_prefix + "/{sample}.a_ctg.gfa",
        hap1_ctg = path_reads_prefix + "/{sample}.hap1.p_ctg.gfa",
        hap2_ctg = path_reads_prefix + "/{sample}.hap2.p_ctg.gfa"
    output:
        primary_ctg_fa = path_reads_prefix + "/{sample}.p_ctg.fa",
        alternate_ctg_fa = path_reads_prefix + "/{sample}.a_ctg.fa",
        hap1_ctg_fa = path_reads_prefix + "/{sample}.hap1.p_ctg.fa",
        hap2_ctg_fa = path_reads_prefix + "/{sample}.hap2.p_ctg.fa"
    shell:
        """
        awk '/^S/{print ">"$2;print $3}' {input.primary_ctg} > {output.primary_ctg_fa} && \ 
        awk '/^S/{print ">"$2;print $3}' {input.alternate_ctg} > {output.alternate_ctg_fa} && \ 
        awk '/^S/{print ">"$2;print $3}' {input.hap1_ctg} > {output.hap1_ctg_fa} && \ 
        awk '/^S/{print ">"$2;print $3}' {input.hap2_ctg} > {output.hap2_ctg_fa}
        """

# Run hifiasm with --dual-scaf and --primary   
# Note: This assumes the output prefix is the same as {wildcards.sample}