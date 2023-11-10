# Define variables

rule adapter_removal:
    input:
        read1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD.fq.gz",
        read2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV.fq.gz"
    output:
        out1 = "{path_reads_prefix}/{species}_{sex}_{method}_FWD.trimmed.fq.gz",
        out2 = "{path_reads_prefix}/{species}_{sex}_{method}_REV.trimmed.fq.gz"
    log:
        out1_log = "logs/{species}_{sex}_{method}_FWD.trimming.log",
        out2_log = "logs/{species}_{sex}_{method}_REV.trimming.log"
    conda:
        "rules/quality_control_reads.yaml"  # Make sure this points to the correct Conda environment file
    shell:
        """
        AdapterRemoval \
            --file1 {input.read1} \
            --file2 {input.read2} \
            --gzip \
            --output1 {output.out1} \
            --output2 {output.out2} \
            --threads {threads} > {log.out1_log} 2> {log.out2_log}
        """

#rule all_adapter_removal:
#    input:
#       expand(f"{path_prefix_reads}/{{sample}}_trimmed.fq.gz", sample=get_combined_samples())




# Note: Make sure that 'config' is loaded either in your main Snakefile or in this file.

