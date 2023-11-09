### Parse and prepare ONT nanopore data

### Parse and prepare Illumina data



### Parse and prepare HiC data

rule fastqc_on_hic:
    input:
        hic_zg=lambda wildcards: f"{path_reads_prefix}/{get_sample(wildcards)}.gz"
    output:
        fastqc_out="{path_out_prefix}/00-FASTQC_TMP/{species}_{sex}_{method}_{orientation}/"
    params:
        tmp_dir="{path_tmp_prefix}/00-FASTQC_TMP/{species}_{sex}_{method}_{orientation}/",
        log="logs/fastqc_on_hic_{wildcards.species}_{wildcards.sex}_{wildcards.method}_{wildcards.orientation}.log"
    conda:
        "quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        echo "Starting FastQC on HiC data." > {params.log}
        mkdir -p {params.tmp_dir} {output.fastqc_out} && \
        fastqc {input.hic_zg} -d {params.tmp_dir} -o {output.fastqc_out} &>> {params.log}
        echo "FastQC on HiC data completed." >> {params.log}
        """

### Parse and prepare PacBio data

### Parse and prepare PacBio data
rule bam_to_fastq:
    input:
        bam=lambda wildcards: f"{path_reads_prefix}/{get_sample(wildcards)}.bam"
    output:
        fastq="{path_reads_prefix}/{species}_{sex}_{method}_{orientation}_rawreads.fastq"
    params:
        log="logs/bam_to_fastq_{wildcards.species}_{wildcards.sex}_{wildcards.method}_{wildcards.orientation}.log"
    conda:
        "quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        echo "Converting BAM to FASTQ." > {params.log}
        if [ "{wildcards.method}" == "HIFI" ]; then
            bam2fastq -o {output.fastq} {input.bam} &>> {params.log}
        else
            echo "Unsupported method: {wildcards.method}" >> {params.log}
            exit 1
        fi
        echo "Conversion completed." >> {params.log}
        """

rule remove_hifi_adapters:
    input:
        raw_read_hifi=lambda wildcards: f"{path_reads_prefix}/{get_sample(wildcards)}_rawreads.fastq"
    output:
        trimmed_hifi_reads="{path_reads_prefix}/{sample}_reads.fastq"
    params:
        log="logs/remove_hifi_adapters_{sample}.log"
    conda:
        "envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """ 
        echo "Starting adapter removal for HiFi reads." > {params.log}
        if [ "{wildcards.method}" == "HIFI" ]; then
            seqtk trimfq -b 20 -e 20 {input.raw_read_hifi} > {output.trimmed_hifi_reads} &>> {params.log}
        else
            echo "Unsupported method: {wildcards.method}" >> {params.log}
            exit 1
        fi
        echo "Adapter removal completed." >> {params.log}
        """

rule gzip_hifi_fastq:
    input:
        trimmed_hifi_reads=lambda wildcards: f"{path_reads_prefix}/{get_sample(wildcards)}_reads.fastq"
    output:
        fastq="{path_reads_prefix}/{species}_{sex}_{method}_{orientation}_reads.fastq.gz"
    params:
        log="logs/gzip_hifi_fastq_{wildcards.species}_{wildcards.sex}_{wildcards.method}_{wildcards.orientation}.log"
    conda:
        "quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """ 
        echo "Starting gzipping of trimmed HiFi reads." > {params.log}
        if [ "{wildcards.method}" == "HIFI" ]; then
            gzip -5 {input.trimmed_hifi_reads} &>> {params.log}
        else
            echo "Unsupported method: {wildcards.method}" >> {params.log}
            exit 1
        fi
        mv {input.trimmed_hifi_reads}.gz {output.fastq}
        echo "Gzipping completed." >> {params.log}
        """


### Mash

rule mash_sketch:
    input:
        fastq_gz="{path_reads_prefix}/{sample}_reads.fastq.gz", 
    output:
        sketch="{sample}.msh"    
    shell:
        "mash sketch -m 2 -o {output.sketch} {input.fastq_gz}"


rule mash_screen: # Run on the reads against RefSeq minimal database 
    input:
        ref_sketch="{path_data_prefix}/01-MASHDB/combined.msh", # The RefSeq database of mash indexes.
        reads="{path_reads_prefix}/{sample}_reads.fastq.gz"
    output:
        screen="{path_out_prefix}/00-MASH/{sample}_screen.tab"
    shell:
        "mash screen {input.ref_sketch} {input.reads} > {output.screen}"

# Rule mash_dist: calculates the pairwise Mash distances between a reference sketch and read sketches, 
#  saving the output in a table, 
#  generates a key file that strips away extraneous filename components to produce a clean label for each sample.

rule mash_dist: 
    input:
        ref_sketch="{path_data_prefix}/00-MASH_DB/combined.msh",
        reads_sketch="{path_reads_prefix}/{sample}_reads.fastq.gz.msh"
    output:
        distances="{path_out_prefix}/00-MASH/{sample}_combined.tbl",
        key="{path_out_prefix}/00-MASH/{sample}_mash_dist_keyfile.txt"
    shell:
        """
        mash dist -p {threads} {input.ref_sketch} {input.reads_sketch} > {output.distances}
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
        distance_file="{path_out_prefix}/00-MASH/{sample}_combined.tbl",
        key_file="{path_out_prefix}/00-MASH/{sample}_mash_dist_keyfile.txt"
    output:
        plot="{path_out_prefix}/00-MASH/{sample}_mash_plot.png"
    params:
        plot_mash_script="workflow/scripts/plot_mash.R"  # Adjust this to your script path
    conda:
        "envs/quality_control_reads.yaml"  # Replace with the path to your conda environment file
    shell:
        """
        Rscript {params.plot_mash_script} {input.distance_file} {input.key_file} {output.plot}
        """
