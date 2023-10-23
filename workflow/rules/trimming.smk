rule adapter_removal:
    input:
        read1 = "{path_prefix_reads}{sample}_R1.fastq.gz",
        read2 = "{path_prefix_reads}{sample}_R2.fastq.gz"
    output:
        out1 = "{path_prefix_reads_trimmed}{sample}_trimmed_R1.fastq.gz",
        out2 = "{path_prefix_reads_trimmed}{sample}_trimmed_R2.fastq.gz"
    threads: config['CPUS']
    shell:
        """
        mkdir -p {config[path_prefix_step1_trimming]}AdaptRemoval
        {config[path_prefix_programs]}{config[progAdapterRemoval]} \
            --file1 {input.read1} \
            --file2 {input.read2} \
            --gzip \
            --output1 {output.out1} \
            --output2 {output.out2} \
            --threads {threads}
        """

rule all_adapter_removal:
    input:
        expand(f"{config['path_prefix_step1_trimming']}AdaptRemoval/{{sample}}_R1_trimmed.fq.gz", sample=config['shortNames'].keys()),
        expand(f"{config['path_prefix_step1_trimming']}AdaptRemoval/{{sample}}_R2_trimmed.fq.gz", sample=config['shortNames'].keys())



# Note: Make sure that 'config' is loaded either in your main Snakefile or in this file.

