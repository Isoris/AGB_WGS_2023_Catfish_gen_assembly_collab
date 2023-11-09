rule align_with_mummer:
    input:
        assembly = "assembly_unicycler.fasta",
        ref_genome = lambda wildcards: config['ref_genomes'][wildcards.ref]
    output:
        delta = f"{config['output'][config['current_step']]}/{wildcards.ref}_output.delta",
        dnadiff_out = f"{config['output'][config['current_step']]}/{wildcards.ref}_dnadiff.report",
        plot = f"{config['output'][config['current_step']]}/{wildcards.ref}_output.png"
    params:
        outdir = config['output'][config['current_step']]
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run MUMmer to compare assemblies
        nucmer -p {params.outdir}/{wildcards.ref}_output {input.assembly} {input.ref_genome}

        # Run dnadiff for additional statistics (optional)
        dnadiff -p {params.outdir}/{wildcards.ref}_dnadiff {input.ref_genome} {input.assembly}

        # Generate a plot of the comparison
        mummerplot -p {params.outdir}/{wildcards.ref}_output {output.delta}
        """

rule all_align_with_mummer:
    input:
        expand("{ref}_output.delta", ref=config['ref_genomes'].keys()),
        expand("{ref}_dnadiff.report", ref=config['ref_genomes'].keys()),
        expand("{ref}_output.png", ref=config['ref_genomes'].keys())

