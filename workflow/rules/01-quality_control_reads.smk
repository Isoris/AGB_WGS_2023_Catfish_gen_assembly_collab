### Mash

rule mash_sketch:
    input:
        fastq_gz = path_reads_prefix + "/{species}_{sex}_{method}_{orientation}_reads.fastq.gz", 
    output:
        sketch = path_data_prefix + "/00-MASH_DB/{species}_{sex}_{method}_{orientation}_sketch_reads.msh"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file    
    shell:
        "mash sketch -m 2 -o {output.sketch} {input.fastq_gz}"

rule mash_dist: 
    input:
        ref_sketch = path_data_prefix + "/00-MASH_DB/combined.msh",
        reads_sketch = path_data_prefix + "/00-MASH_DB/{species}_{sex}_{method}_{orientation}_sketch_reads.msh"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    output:
        distances = path_out_prefix + "/00-MASH/{species}_{sex}_{method}_{orientation}_combined.tbl",
        key = path_out_prefix + "/00-MASH/{species}_{sex}_{method}_{orientation}_mash_dist_keyfile.txt"
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
        distance_file = path_out_prefix + "/00-MASH/{species}_{sex}_{method}_{orientation}_combined.tbl",
        key_file = path_out_prefix +"/00-MASH/{species}_{sex}_{method}_{orientation}_mash_dist_keyfile.txt"
    output:
        plot = path_out_prefix + "/00-MASH/{species}_{sex}_{method}_{orientation}_mash_plot.png"
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
        ref_sketch = path_data_prefix + "/01-MASHDB/combined.msh", # The RefSeq database of mash indexes.
        reads = path_reads_prefix + "/{species}_{sex}_{method}_{orientation}_reads.fastq.gz"
    output:
        screen = path_out_prefix + "/00-MASH/{species}_{sex}_{method}_{orientation}_screen.tab"
    conda:
        "../envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"

# Rule mash_dist: calculates the pairwise Mash distances between a reference sketch and read sketches, 
#  saving the output in a table, 
#  generates a key file that strips away extraneous filename components to produce a clean label for each sample.
