# Define rules for running quast lg --large on gfa assemblies from 02-
rule quast_lg_assess_assembly:
    input:
	spades_assembly = "{path_prefix_step2_spades}/output_spades_{species}_{sex}_{version}/XXX_spades.fasta",
        ref_genome = "{path_prefix_ref_genomes}/{species}_ref.fasta",
        gff = "{path_prefix_ref_genomes}/{species}.gff",
        pe1 = "{path_prefix_reads_trimmed}/{sample}_reads_trim1.fq",
        pe2 = "{path_prefix_reads_trimmed}/{sample}_reads_trim2.fq"
    output:
	report = "Output_QUAST_{sample}_{version_num}"
    threads: config['CPUS']
    shell:
	"""
	echo "Start of QUAST quality control of the Spades de-novo assemblies for the data"
        echo "Run QUAST with Busco and silva and gridss and rna+gene finding, report all metrics, create circos"
        echo "QUAST Quality control step starting for {wildcards.sex} {wildcards.species}"

        {config[progQUAST]} \
            {input.spades_assembly} \
            -r {input.ref_genome} \
            -g {input.gff} \
            --pe1 {input.pe1} \
            --pe2 {input.pe2} \
            --circos \
            --threads {threads} \
            --gene-finding \
            --rna-finding \
            --conserved-genes-finding \
            --report-all-metrics \
            -L \
            -o {output.report}

        echo "QUAST Quality control step finished for {wildcards.species} {wildcards.sex} {wildcards.version}"
        """


rule read_genomescope_summary_illum_out:
    input:
        species_summary_genomescope = path_out_prefix + "/01-GENOMESCOPE/{sample}_ILLUMINA/summary.txt"
    output:
        haploid_length = path_out_prefix + "/01-GENOMESCOPE/{sample}_ILLUMINA/{sample}_haploid_length.txt"
    shell: 
    """
    grep 'Genome Haploid Length' {input.summary_file} | cut -d ' ' -f 3 > {output.haploid_length}
    """

rule read_genomescope_summary_hifi_out:
    input:
        species_summary_genomescope = path_out_prefix + "/01-GENOMESCOPE/{sample}_HIFI/summary.txt"
    output:
        haploid_length = path_out_prefix + "/01-GENOMESCOPE/{sample}_HIFI/{sample}_haploid_length.txt"
    shell: 
    """
    grep 'Genome Haploid Length' {input.summary_file} | cut -d ' ' -f 3 > {output.haploid_length}
    """

def get_haploid_length_path(wildcards):
    hifi_path = path_out_prefix + f"/01-GENOMESCOPE/{wildcards.sample}_HIFI/summary.txt"
    illumina_path = path_out_prefix + f"/01-GENOMESCOPE/{wildcards.sample}_ILLUMINA/summary.txt"

    if os.path.exists(hifi_path):
        return hifi_path
    elif os.path.exists(illumina_path):
        return illumina_path
    else:
        raise ValueError(f"No Genomescope data found for {wildcards.sample}")



rule meryl_build_k:
    input:
        haploid_length = path_out_prefix + "/01-GENOMESCOPE/{sample}_ILLUMINA/{sample}_haploid_length.txt"
    output:
        kmer_val = path_out_prefix + "/{sample}_meryl_kmer_val.txt"
    params:
        path_build_k = path_prog_prefix + "merqury-1.3/best_k.sh",
    shell:
        """
        # Use the haploid length in the build_k step
        haploid_length=$(cat {input.haploid_length})
        {params.path_build_k} $haploid_length > {output.kmer_val}
        """


rule meryl_to_find_kmersize:
    output: "{species}_meryl_kmer_val.txt"
    params:
        path_build_k = path_prog_prefix + "merqury-1.3/best_k.sh",
        estimated_genome_size_bp = 980000000
        species = config["samples"]["species"]
    shell: 
        """
        {params.path_build_k} {params.estimated_genome_size_bp}" > {output}
        """


rule meryl_build_db:
    input:
    output:
    params:
    shell:


rule :
    input:
    output:
    params:
    shell:

rule :
    input:
    output:
    params:
    shell:

rule :
    input:
    output:
    params:
    shell:
