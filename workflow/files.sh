#! /bin/bash
##### load paths #####

path_folder_prefix=tarafs/scratch/proj5034-AGBKU/catfish
###
path_reads_prefix={path_folder_prefix}/00-DATA/00-READS
path_data_prefix={path_folder_prefix}/00-DATA
path_prog_prefix={path_folder_prefix}/01-PROGRAMS
path_exp_prefix={path_folder_prefix}/02-EXPERIMENTS
path_out_prefix={path_folder_prefix}/03-OUTPUTS
path_tmp_prefix={path_folder_prefix}/04-TEMPORARY  




echo "Listing files in READS directory:"; find {path_reads_prefix} -type f -printf '%P\\n';  \
echo "\\nListing files in DATA directory:"; find {path_data_prefix} -type f -printf '%P\\n';  \
echo "\\nListing files in PROGRAMS directory:"; find {path_prog_prefix} -type f -printf '%P\\n';  \
echo "\\nListing files in EXPERIMENTS directory:"; find {path_exp_prefix} -type f -printf '%P\\n';  \
echo "\\nListing files in OUTPUTS directory:"; find {path_out_prefix} -type f -printf '%P\\n';  \
echo "\\nListing files in TEMPORARY directory:"; find {path_tmp_prefix} -type f -printf '%P\\n';   > files.logs
