
# Snakefile
#https://stackoverflow.com/questions/46506784/how-do-the-terms-job-task-and-step-relate-to-each-other/46532581#46532581

#SPECIES_IN = [""]
#SEX_IN = [""]
#METHODS_IN = [""]
#ORIENTATIONS_ALL = [""]


# ... [Other parts of your Snakefile] ...

# Rule to write config contents to a file
rule debug_config:
    output:
        debug_config_output="debug_config_output.txt"
    run:
        with open(output.debug_config_output, 'w') as f:
            # Iterate over samples and extract the information
            for species, samples_dict in config['samples'].items():
                for sample_id, sample_info in samples_dict.items():
                    # Write out the details you're interested in
                    f.write(f"Species: {species}, Sample ID: {sample_id}\n")
                    f.write(f"Coverage: {sample_info['coverage']}\n")
                    f.write(f"Description: {sample_info['description']}\n")
                    f.write(f"File Prefix: {sample_info['file_prefix']}\n")
                    f.write(f"Method: {sample_info['method']}\n")
                    f.write(f"Orientation: {sample_info['orientation']}\n")
                    f.write(f"Path Prefix: {sample_info['path_prefix']}\n")
                    f.write(f"Sex: {sample_info['sex']}\n")
                    f.write(f"Species Long: {sample_info['species_long']}\n")
                    f.write("\n")

##### Main rule to orchestrate the workflow #####
rule create_log_folder:
    output:
        log="logs/logs.txt"
    shell: "touch {output.log}"

rule log_reads_directory:
    input:
        "paths.sh"
    output:
        log="logs/reads_log.txt"
    shell:
        """
        source paths.sh
        echo "Listing .fq, .fastq, .gz, .sam, .bam, .bed, .msh, .fasta, .gfa files in READS directory:" > {output.log}
        find $path_reads_prefix -type f | grep -E '(\.fq|\.fastq|\.sam|\.bam|\.bed|\.msh|\.fasta|\.gfa)$' >> {output.log}
        echo "" >> {output.log}
        """

rule log_data_directory:
    input:
        "paths.sh"
    output:
        log="logs/data_log.txt"
    shell:
        """
        source paths.sh
        echo "Listing .fq, .fastq, .gz, .sam, .bam, .bed, .msh, .fasta, .gfa files in DATA directory:" > {output.log}
        find $path_data_prefix -type f | grep -E '(\.fq|\.fastq|\.sam|\.bam|\.bed|\.msh|\.fasta|\.gfa)$' >> {output.log}
        echo "" >> {output.log}
        """

rule log_experiments_directory:
    output:
        log="logs/experiments_log.txt"
    shell:
        """
        source paths.sh
        echo "Listing .fq, .fastq, .gz, .sam, .bam, .bed, .msh, .fasta, .gfa files in EXPERIMENTS directory:" > {output.log}
        find $path_exp_prefix -type f | grep -E '(\.fq|\.fastq|\.sam|\.bam|\.bed|\.msh|\.fasta|\.gfa)$' >> {output.log}
        echo "" >> {output.log}
        """

rule log_programs_directory:
    output:
        log="logs/programs_log.txt"
    shell:
        """
        source paths.sh
        echo "Listing .fq, .fastq, .gz, .sam, .bam, .bed, .msh, .fasta, .gfa files in PROGRAMS directory:" > {output.log}
        find path_prog_prefix -type f | grep -E '(\.fq|\.fastq|\.sam|\.bam|\.bed|\.msh|\.fasta|\.gfa)$' >> {output.log}
        echo "" >> {output.log}
        """

rule concatenate_logs:
    input:
        reads="logs/reads_log.txt",
        data="logs/data_log.txt"
        #outputs="logs/outputs_log.txt"
    output:
        log="logs/file_log.txt"
    shell:
        """
        cat {input.reads} {input.data} > {output.log}
        """
