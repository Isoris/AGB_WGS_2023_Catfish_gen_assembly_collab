rule quast_quality_control:
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
	echo "Start of QUAST quality control of the Spades de-novo assemblies for the illumina data"
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


