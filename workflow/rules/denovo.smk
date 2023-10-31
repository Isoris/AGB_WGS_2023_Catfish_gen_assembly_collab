# Define rules for de-novo assembly of reads

# Define variables to be called within rules


rule run_hifasm_hifi_only:    # https://github.com/chhylp123/hifiasm/releases
    input:
       compressed_hifi_raw_reads = 
    output:
       
    params: 
       z = {tri}
    shell:
    "hifiasm -o {output}.asm -z20 --dual-scaf --primary -t {threads} {input}.fq.gz \ "


rule run_hifasm_hifi_UL:    # https://github.com/chhylp123/hifiasm/releases
    input:
       compressed_hifi_raw_reads = 
    output:
       
    params: 
       z = {tri}
    shell:
    "hifiasm -o {output}.asm -z20 --dual-scaf --primary -t {threads} {input}.fq.gz \ "


rule run_hifasm_hifi_UL_hic:    # https://github.com/chhylp123/hifiasm/releases
    input:
       compressed_hifi_raw_reads = 
    output:
       
    params: 
       z = {tri}
    shell:
    """
    seqtk trimfq -b 20 -e 20 input.fastq > output_trimmed.fastq
    hifiasm --h1 HiC1.fastq --h2 HiC2.fastq --ul ONT.fastq HiFi.fastq \ 
    """

rule run_spades_short_reads_only:    # https://github.com/chhylp123/hifiasm/releases
    input:
       compressed_hifi_raw_reads = 
    output:
       
    params: 
       z = {tri}
    shell:
    "hifiasm -o {output}.asm -z20 --dual-scaf --primary -t {threads} {input}.fq.gz \ "


