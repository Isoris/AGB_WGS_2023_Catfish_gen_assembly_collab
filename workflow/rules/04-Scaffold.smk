rule map_hic_reads_to_assembly:
    input:
        hap1_ctg = path_reads_prefix + "/{sample}.hap1.p_ctg.fa",
        hap2_ctg = path_reads_prefix + "/{sample}.hap2.p_ctg.fa",
        hic1_reads = path_reads_prefix + "/{sample}_HIC_FWD.fastq.gz",
        hic2_reads = path_reads_prefix + "/{sample}_HIC_REV.fastq.gz"
    output:
        mapped_reads = path_out_prefix + "/{sample}_mapped_hic_reads.paf"
    params:
        minimap2_path = path_prog_prefix + "/minimap2/minimap2"
    shell:
        """
        cat {input.hap1_ctg} {input.hap2_ctg} > {path_out_prefix}/{wildcards.sample}_combined_haplotypes.fa
        {params.minimap2_path} -x sr {path_out_prefix}/{wildcards.sample}_combined_haplotypes.fa \
        {input.hic1_reads} {input.hic2_reads} > {output.mapped_reads}
        """

rule create_hic_contact_map:
    input:
        mapped_reads = path_out_prefix + "/{sample}_mapped_hic_reads.paf"
    output:
        hic_contact_map = path_out_prefix + "/{sample}_hic_contact_map.hic"
    params:
        juicer_tools_path = path_prog_prefix + "/juicer_tools/juicer_tools"
    shell:
        """
        {params.juicer_tools_path} pre {input.mapped_reads} {output.hic_contact_map}
        """


rule hic_scaffolding_phasing:
    input:
        combined_haplotypes = path_out_prefix + "/{sample}_combined_haplotypes.fa",
        hifi_reads = path_reads_prefix + "/{sample}_HIFI_None_trimmed_reads.fastq.gz",
        pe1_reads = path_reads_prefix + "/{sample}_HIC_FWD.fastq.gz",
        pe2_reads = path_reads_prefix + "/{sample}_HIC_REV.fastq.gz"
    output:
        phased_assembly = path_out_prefix + "/{sample}_phased_assembly.fa"
    params:
        greenhill_path = path_prog_prefix + "/greenhill/greenhill"
    shell:
        """
        {params.greenhill_path} -cph {input.combined_haplotypes} \
        -p {input.hifi_reads} \
        -IP1 {input.pe1_reads} {input.pe2_reads} 
        """
