#!/bin/bash
#conda create --name quast
#conda install -c bioconda quast
#quast-download-gridss                                                       
#quast-download-silva                                                        
#quast-download-busco

source /path/to/config.sh  # Adjust the path to your config.sh

echo "Start of QUAST quality control of the Spades de-novo assemblies for the illumina data"

for i in ${!shortNames[@]}; do
    sname=${shortNames[$i]}
    sex=${Sex[$i]}
    species=${Species[$i]}
    version=${Version[$i]}
    version_num=${Version_num[$i]}

    echo "Run QUAST with Busco and silva and gridss and rna+gene finding, report all metrics, create circos"

    echo "QUAST Quality control step starting for $sex $species"

    ${progQUAST} \
        ${path_prefix_step2_spades}/output_spades_${species}_${sex}_${version}/XXX_spades.fasta \
        -r ${path_prefix_ref_genomes}/${species}_ref.fasta \
        -g ${path_prefix_ref_genomes}/${species}.gff \
        --pe1 ${path_prefix_reads_trimmed}/${reads_trim1[$i]} \
        --pe2 ${path_prefix_reads_trimmed}/${reads_trim2[$i]} \
        --circos \
        --threads $CPUS \
        --gene-finding \
        --rna-finding \
        --conserved-genes-finding \
        --report-all-metrics \
        -L \
        -o Output_QUAST_${sname}_${version_num}

    echo "QUAST Quality control step finished for $species $sex $version"
done

echo "QUAST Quality control step finished"

