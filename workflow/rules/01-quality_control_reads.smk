### Mash
rule mash_sketch:
    input:
        fastq_gz = path_reads_prefix + "/{sample}_trimmed_reads.fastq.gz", 
    output:
        sketch = path_data_prefix + "/01-MASH_DB/{sample}_sketch_reads.msh"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file    
    shell:
        "mash sketch -k 21 -s 10000 -r -m 2 -o {output.sketch} {input.fastq_gz}"


# Rule mash_dist: calculates the pairwise Mash distances between a reference sketch and read sketches, 
#  saving the output in a table, 
#  generates a key file that strips away extraneous filename components to produce a clean label for each sample.
rule mash_dist: 
    input:
        ref_sketch = path_data_prefix + "/01-MASH_DB/combined.msh",
        reads_sketch = path_data_prefix + "/01-MASH_DB/{species}_{sex}_{method}_{orientation}_sketch_reads.msh"
    output:
        distances = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_combined.tbl",
        key = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_mash_dist_keyfile.txt"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        mash dist -p {threads} {input.ref_sketch} {input.reads_sketch} > {output.distances} && \
        head -n 1 {output.distances} | \
        awk '{{for (i=2; i <=NF; i++) print $i}}' | \
        awk -F "/" '{{print $NF}}' | \
        sed 's/\\.subreads\\.fast[aq]\\.gz//g' | \
        sed 's/_reads\\.fast[aq]\\.gz//g' | \
        sed 's/\\.fast[aq]\\.gz//g' | \
        sed 's/\\.fast[aq]//g'  > {output.key}
        """


rule mash_dist_plot:
    input:
        distance_file = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_combined.tbl",
        key_file = path_out_prefix +"/01-MASH/{species}_{sex}_{method}_{orientation}_mash_dist_keyfile.txt"
    output:
        plot = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_mash_plot.png"
    params:
        plot_mash_script = "workflow/scripts/plot_mash.R"  # Adjust this to your script path
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        Rscript {params.plot_mash_script} {input.distance_file} {input.key_file} {output.plot}
        """

rule mash_screen: # Run on the reads against RefSeq minimal database 
    input:
        ref_sketch = path_data_prefix + "/01-MASH_DB/combined.msh", # The RefSeq database of mash indexes.
        reads = path_reads_prefix + "/{species}_{sex}_{method}_{orientation}_trimmed_reads.fastq.gz"
    output:
        screen = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_screen.tab"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"

### Jellyfish
# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf


rule run_jellyfish_count_illumina_pe:
    input:
        read_fwd = path_reads_prefix + "/{sample}_FWD_trimmed_reads.fastq.gz",
        read_rev = path_reads_prefix + "/{sample}_REV_trimmed_reads.fastq.gz"
    output:
        jellyfish_count = path_out_prefix + "/01-JELLYFISH/{sample}_ILLUMINA/{sample}_ILLUMINA_trimmed_reads_jellyfish_count.jf"
    params:
        jellyfish_path = path_prog_prefix + "/jellyfish-2.3.0/bin/jellyfish",
        kmer_size = 21
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        {params.jellyfish_path} count -m {params.kmer_size} -s 100M -t {threads} -C \
        -o {output.jellyfish_count} <(zcat {input.read_fwd}) <(zcat {input.read_rev})
        """

rule run_jellyfish_histo_illum_pe:
    input:
        jellyfish_count = path_out_prefix + "/01-JELLYFISH/{sample}_ILLUMINA/{sample}_ILLUMINA_trimmed_reads_jellyfish_count.jf"
    output:
        jellyfish_histo = path_out_prefix + "/01-JELLYFISH/{sample}_ILLUMINA/{sample}_ILLUMINA_trimmed_reads_jellyfish_histo.histo"
    params:
        jellyfish_path = path_prog_prefix + "/jellyfish-2.3.0/bin/jellyfish"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        {params.jellyfish_path} histo -t {threads} {input.jellyfish_count} > {output.jellyfish_histo}
        """


rule run_jellyfish_count_hifi:
    input:
        reads = path_reads_prefix + "/{sample}_HIFI_None_trimmed_reads.fastq.gz"
    output:
        jellyfish_count= path_out_prefix + "/01-JELLYFISH/{sample}_HIFI/{sample}_HIFI_None_trimmed_reads_jellyfish_count.jf"
    params:
        jellyfish_path = path_prog_prefix + "/jellyfish-2.3.0/bin/jellyfish",
        kmer_size = 21
    conda:
        "../envs/quality_control_reads.yaml"    # Replace with the path to your conda environment file
    shell:
        """
        zcat {input.reads} | {params.jellyfish_path} count -m {params.kmer_size} -s 100M -t {threads} \
        -C -o {output.jellyfish_count} /dev/stdin
        """

rule run_jellyfish_histo_hifi:
    input:
        jellyfish_count = path_out_prefix + "/01-JELLYFISH/{sample}_HIFI/{sample}_HIFI_None_trimmed_reads_jellyfish_count.jf"
    output:
        jellyfish_histo = path_out_prefix + "/01-JELLYFISH/{sample}_HIFI/{sample}_HIFI_None_trimmed_reads_jellyfish_histo.histo"
    params:
        jellyfish_path = path_prog_prefix + "/jellyfish-2.3.0/bin/jellyfish"
    conda:
        "../envs/quality_control_reads.yaml"    # Replace with the path to your conda environment file
    shell:
        """
        {params.jellyfish_path} histo -t {threads} {input.jellyfish_count} > {output.jellyfish_histo}
        """

### GenomeScope
# Role: Genomescope is used for genome profiling based on k-mer analysis
rule run_genomescope_hifi:
    input:
        jellyfish_histo = path_out_prefix + "/01-JELLYFISH/{sample}_HIFI/{sample}_HIFI_None_trimmed_reads_jellyfish_histo.histo"
    output:
        directory(path_out_prefix + "/01-GENOMESCOPE/{sample}_HIFI/")
    params:
        genomescope_script="scripts/genomescope.R",  # adjust this to your genomescope script path
        kmer_length=25  # adjust as needed
    conda:
        "../envs/quality_control_reads.yaml"    # Replace with the path to your conda environment file
    shell:
        """
        Rscript {params.genomescope_script} {input.jellyfish_histo} {params.kmer_length} {params.read_length} {output}
        """

rule run_genomescope_illum_pe:
    input:
        jellyfish_histo = path_out_prefix + "/01-GENOMESCOPE/{sample}_ILLUMINA_trimmed_reads_jellyfish_histo.histo"
    output:
        directory(path_out_prefix + "/01-GENOMESCOPE/{sample}_ILLUMINA/")
    params:
        genomescope_script="scripts/genomescope.R",  # adjust this to your genomescope script path
        kmer_length=21  # adjust as needed
    conda:
        "../envs/quality_control_reads.yaml"    # Replace with the path to your conda environment file
    shell:
        """
        Rscript {params.genomescope_script} {input.jellyfish_histo} {params.kmer_length} {output}
        """




#{params.kmer_max if params.kmer_max else ''} {'--verbose' if params.verbose else ''}
# Role: 
#rule run_genomescope:
#    input:
#        "{sample}.jq"
#    output:
#        directory("path/to/output/{sample}_LongQC_output/")
#    params:
#        jellyfish_path="/tarafs/data/home/qandres/scratch/catfish/01-PROGRAMS/LongQC-1"  # adjust this to your jellyfish installation path
#    conda:
#        "/quality_control_reads.yaml"  # Replace with the path to your 'quality_control' environment file
#    shell:
#        """
#        cd {jellyfish_output_path}  # Navigating to jellyfish output directory
#         Rscript genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose] 
#        """            

# Define rule for running LongQC on raw reads (longCQ is under catfish/01-PROGRAMS/) should run within mamba quality_control environment 
# Role: For the QC of long reads from HiFi

# https://github.com/VGP/vgp-assembly/tree/master/pipeline/meryl

# https://github.com/marbl/merqury
# https://github.com/cnag-aat/assembly_pipeline/tree/v2.1.0#1--ont-reads