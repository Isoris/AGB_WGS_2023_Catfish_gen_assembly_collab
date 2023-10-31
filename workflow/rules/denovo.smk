# Define rules for de-novo assembly of reads

# Define variables to be called within rules


rule run_hifasm_hifi_UL_hic:
    input:
        compressed_hifi_raw_reads = "path/to/HiFi_reads.fastq",
        ont_reads = "path/to/ONT_reads.fastq",
        hic1_reads = "path/to/HiC1.fastq",
        hic2_reads = "path/to/HiC2.fastq"
    output:
        primary_ctg = "path/to/output/{sample}.p_ctg.gfa",
        alternate_ctg = "path/to/output/{sample}.a_ctg.gfa",
        hap1_ctg = "path/to/output/{sample}.hap1.p_ctg.gfa",
        hap2_ctg = "path/to/output/{sample}.hap2.p_ctg.gfa"
    params:
        hifiasm_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/hifiasm-0.19.7/hifiasm",
        z = "{tri}"  # This assumes 'tri' is defined somewhere else in your workflow
    log:
        "logs/hifiasm/{sample}.log"
    shell:
        """
        # Trim the HiFi reads
        seqtk trimfq -b 20 -e 20 {input.compressed_hifi_raw_reads} > {input.compressed_hifi_raw_reads}.trimmed.fastq

        # Run hifiasm with --dual-scaf and --primary
        {params.hifiasm_path} --dual-scaf --primary --h1 {input.hic1_reads} --h2 {input.hic2_reads} --ul {input.ont_reads} {input.compressed_hifi_raw_reads}.trimmed.fastq -z {params.z} -o path/to/output/{wildcards.sample}

        # Note: This assumes the output prefix is the same as {wildcards.sample}
        """

rule run_spades_short_reads_only:    # https://github.com/chhylp123/hifiasm/releases
    input:
       compressed_hifi_raw_reads = 
    output:
       
    params: 
       z = {tri}
    shell:
    "hifiasm -o {output}.asm -z20 --dual-scaf --primary -t {threads} {input}.fq.gz \ "


