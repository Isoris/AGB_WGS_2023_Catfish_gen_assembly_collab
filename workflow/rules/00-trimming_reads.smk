



# Define variables
rule adapter_removal:
    input:
        "paths.sh",
        read1 = "{path_prefix_reads}/{sample}_FWD.fq.gz",
        read2 = "{path_prefix_reads}/{sample}_REV.fq.gz"
    output:
        out1 = "{path_prefix_reads}/{sample}_FWD.trimmed.fq.gz",
        out2 = "{path_prefix_reads}/{sample}_REV.trimmed.fq.gz"
    conda:
        "quality_control_reads.yaml"  # Make sure this points to the correct Conda environment file
    shell:
        """
        source paths.sh
        AdapterRemoval \
            --file1 {input.read1} \
            --file2 {input.read2} \
            --gzip \
            --output1 {output.out1} \
            --output2 {output.out2} \
            --threads {threads}
        """

#rule all_adapter_removal:
#    input:
#       expand(f"{path_prefix_reads}/{{sample}}_trimmed.fq.gz", sample=get_combined_samples())




# Note: Make sure that 'config' is loaded either in your main Snakefile or in this file.

