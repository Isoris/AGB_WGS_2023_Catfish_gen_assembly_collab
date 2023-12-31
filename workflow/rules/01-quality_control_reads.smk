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
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    output:
        distances = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_combined.tbl",
        key = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_mash_dist_keyfile.txt"
    shell:
        """
        mash dist -p {threads} {input.ref_sketch} {input.reads_sketch} > {output.distances} && \
        head -n 1 {output.distances} | \
        awk '{for (i=2; i <=NF; i++) print $i}' | \
        awk -F "/" '{print $NF}' | \
        sed 's/\.subreads\.fast[aq]\.gz//g' | \
        sed 's/_reads\.fast[aq]\.gz//g' | \
        sed 's/\.fast[aq]\.gz//g' | \
        sed 's/\.fast[aq]//g'  > {output.key} 
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
        reads = path_reads_prefix + "/{species}_{sex}_{method}_{orientation}_reads.fastq.gz"
    output:
        screen = path_out_prefix + "/01-MASH/{species}_{sex}_{method}_{orientation}_screen.tab"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"

### Jellyfish
# Role: Used for the quality control of reads by fast, memory-efficient counting of k-mers in DNA https://genome.umd.edu/jellyfish.html 
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

rule run_jellyfish_count:
    input:
        reads=f"{path_reads_prefix}/{{sample}}_rawreads.fastq.gz"
    output:
        jellyfish_count=f"{path_out_prefix}/{{sample}}_jellyfishcount.jf"
    params:
        jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish",
        kmer_size = 21
    conda:
        "quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        zcat {input.reads} | {params.jellyfish_path} count -m {params.kmer_size} -s 100M -t {threads} \
        -C -o {output.jellyfish_count} /dev/stdin
        """

rule run_jellyfish_histo:
    input:
        jellyfish_count=f"{path_out_prefix}/{{sample}}_jellyfishcount.jf"
    output:
        histo=f"{path_out_prefix}/{{sample}}_jellyfish_histo.histo"
    params:
         jellyfish_path = "/tarafs/data/home/qandres/catfish/catfish/01-PROGRAMS/jellyfish-2.3.0/bin/jellyfish"
    conda:
        "quality_control.yml"  # Replace with the path to your 'quality_control' environment file
    shell:
        "{params.jellyfish_path} histo -t {threads} {input.jellyfish_count} > {output.histo}"


### GenomeScope
# Role: Genomescope is used for genome profiling based on k-mer analysis
rule run_genomescope:
    input:
        jellyfish_histo=f"{path_out_prefix}/{{sample}}_jellyfish_histo.histo"
    output:
        directory(f"{path_out_prefix}/{{sample}}_genomescope_output/")
    params:
        genomescope_script="scripts/genomescope.R",  # adjust this to your genomescope script path
        kmer_length=21,  # adjust as needed
        read_length=100,  # adjust based on your read length
        kmer_max=None,  # optional, set a value if required
        verbose=False  # set to True if verbose output is needed
    conda:
        "quality_control_reads.yaml"  # Replace with the path to your 'quality_control' environment file
    shell:
        """
        Rscript {params.genomescope_script} {input.jellyfish_histo} {params.kmer_length} {params.read_length} {output} {params.kmer_max if params.kmer_max else ''} {'--verbose' if params.verbose else ''}
        """

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

