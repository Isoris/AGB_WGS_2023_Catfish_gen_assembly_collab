#! /bin/bash/


#################################################### CHECK FOR SOFTWARE PATHS BEFORE SETTING UP CONDA
# Create or empty the progPaths.log file for fresh logging
> progPaths_before.log

# Add a header to the log
echo "Verifying if the programs are available in the HPC and capturing their paths:" >> progPaths_before.log

# Check for all Python installations
echo "All Python installations:" >> progPaths_before.log
which -a python >> progPaths_before.log

# Check for the Python executable in the PATH
echo "Default Python installation:" >> progPaths_before.log
which python >> progPaths_before.log






##################################################### PREPARE THE CONDA ENVIRONMENT 
# install conda forge mamba on the base env of conda in the HPC
conda install -n base -c conda-forge mamba
# create the new environment "snakemake" in conda 
mamba create -c conda-forge -c bioconda -n snakemake snakemake 	
mamba init
source ~/.bashrc
# activate snakemake conda / mamba env.
mamba activate snakemake
# check if everything works as intended
snakemake --help
# Check for packages and versions in the mamba snakemake environment 
> mamba_snakemake_packages.log
echo " packages and versions in the mamba snakemake environment" >> mamba_snakemake_packages.log
mamba list >> mamba_snakemake_packages.log


#################################################### CHECK FOR SOFTWARE PATHS AFTER SETTING UP CONDA (NEED TO BE WITHIN SNAKEMAKE ENV.)
# Create or empty the progPaths.log file for fresh logging
> progPaths_after_conda.log

# Add a header to the log
echo "Verifying if the programs are available in the HPC and capturing their paths:" >> progPaths_after_conda.log

# Check for all Python installations
echo "All Python installations:" >> progPaths_after_conda.log
which -a python >> progPaths_after_conda.log

# Check for the Python executable in the PATH
echo "Default Python installation:" >> progPaths_after_conda.log
which python >> progPaths_after_conda.log











































#################################################### PATHS 
#Path to the folder containing programs
progFolder="./02-PROGRAMS/"

#Path to systems, workflow management tools and dependencies
##Snakemake for workflow management 
progSnakemake=										#	https://snakemake.readthedocs.io/en/stable/tutorial/setup.html
##Dependencies for snakemake
progPysam=										#	https://pysam.readthedocs.io/en/stable/
progPython_v3_5=									#	https://www.python.org/
progMatplotlib= 									#	https://matplotlib.org/
progNetworkx=										#	https://networkx.org/
progGraphviz=										#	https://www.graphviz.org/
progJinja2=										#	https://jinja.palletsprojects.com/en/3.1.x/



# Program Paths for Quality Control of Reads
## Pipeline Mash								#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/mash
progMash=										#	https://github.com/marbl/Mash #calculate mash distances
##------------------------- 
progMercury=										#	https://github.com/marbl/merqury

# Program Paths for Adapter Trimming
progAdapterRemoval="/path/to/AdapterRemoval"						#	https://github.com/MikkelSchubert/adapterremoval #remove adapters and trimm.

# Program Paths for De-novo Illumina Assembly
progSpades="/path/to/Spades"								#	https://github.com/ablab/spades #denovo assembly of short-reads

# Program Paths for Assembly Quality Control
progfastqc="/path/to/FastQC"								#	https://github.com/s-andrews/FastQC #tools to evaluate reads
progQUAST="/path/to/QUAST"								#	https://github.com/ablab/quast #tools to evaluate assemblies 

# Program Paths for Short Reads Illumina Mapping
progBWAMEM="/path/to/BWAMEM"								#	https://github.com/bwa-mem2/bwa-mem2 #efficient short-reads mapping
progBowtie2="/path/to/Bowtie2"								#	https://github.com/BenLangmead/bowtie2 #efficient short-reads mapping

# Program Paths for Working with Alignment Files and Variant Calling
progSamtools="/path/to/Samtools"							#	https://github.com/samtools/samtools #manipulation of SAM files
progBEDTools="/path/to/BEDTools"							#	https://github.com/arq5x/bedtools2 #extract genomic ranges and intervals
progBamtools=										#	https://github.com/pezmaster31/bamtools #manipulation of BAM files
progBcftools="/path/to/Bcftools"							#	https://github.com/samtools/bcftools #variant calling 
progVCFtools="/path/to/VCFtools"							#	https://github.com/vcftools/vcftools #manipulation of vcf files
progVCFutils="/path/to/VCFutils"							#	https://github.com/druths/vcfutils # manipulation of vcf files
progVcflib=										#	https://github.com/vcflib/vcflib/tree/master  #its also vcffirstheader	
progPILON=										#	https://github.com/broadinstitute/pilon #variant calling alt. freebayes
progPicard=										#	https://github.com/broadinstitute/picard #extract reads from SAM BAM
progSeqtk=										#	https://github.com/lh3/seqtk #Processing of .fasta and .fastq
progGatk= #${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar 			#	
progFastaStats= #FastaStats.jar								#	
 
#### Program Paths for Long Reads Pacbio Assembly
##Pipeline Triocanu								#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/triocanu
progCanu=										#	https://github.com/marbl/canu # assembler for pacbio and for nanopore
progMhap=										#	https://github.com/marbl/MHAP
progArrow=										#	https://github.com/skoren/ArrowGrid
progArrowgrid=										#	https://github.com/skoren/ArrowGrid
progPbbam=  										#	https://github.com/PacificBiosciences/pbbam
progPbgcpp=										#	https://github.com/PacificBiosciences/pbbioconda 
progVerkko=										#	https://github.com/marbl/verkko#verkko #denovo for diploids
progHifiasm=										#	https://github.com/chhylp123/hifiasm #denovo 
progPbsv=										#	https://github.com/PacificBiosciences/pbsv #mapping of pacbio to reference

## Pipeline minimap2 + purgedups						#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups
progPurge_dups_v1_6=									#	https://github.com/dfguan/purge_dups #remove PCR duplicates

#### Program Paths for Haplotype Phasing
progPurgehaplotypes_v1_5=								#	
progFalcon=										#	https://github.com/PacificBiosciences/FALCON-integrate/wiki/Installation
progFalconunzip=									#	https://github.com/PacificBiosciences/FALCON_unzip/wiki/Binaries
#### Program Paths for HiC Assembly step 1
progMeryl=										#	https://github.com/marbl/meryl
### Pipeline ArimaGenomics							#	https://github.com/ArimaGenomics/mapping_pipeline
progSalsa_v2_2=										#	https://github.com/marbl/SALSA 
progNtrimming=										#
### Pipeline Pretext								#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/pretext
progPretexmap=										#	https://github.com/wtsi-hpag/PretextMap  #visualize the HiC mapping
progPretextview=									#	https://github.com/wtsi-hpag/PretextView #visualize the HiC mapping

#### Program Paths for Mitochondrial Genome Assembly from PacBio
progBlasR=										#	https://github.com/BioinformaticsArchive/blasr
progBlastN=										#	https://www.ncbi.nlm.nih.gov/books/NBK52640/
progFreebayes=${progFolder}freebayes-1.3.7						#	https://github.com/freebayes/freebayes  # variant calling 
progMitoS2annot=									#	
progArraytrimming=									#	
progEndtrimming=									#

#### Program Paths for HiC Assembly step 2 / Freebayes polish			#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish
progLongranger=										# 	https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation
progmat=										#
progmitoVGP_v1_6=									#	https://github.com/gf777/mitoVGP
progMedaka=										#	https://github.com/nanoporetech/medaka
progRacon=										#	https://github.com/isovic/racon #polishing of assemblies from nanopore.
progHtslib=										#	https://github.com/samtools/htslib             includes tabix

## Post evaluation
#### Program Paths for Whole Genome Alignments
progMuMmer="/path/to/MuMmer"								#	https://github.com/mummer4/mummer one-to-one genome alignments
progNucmer=									#  	https://github.com/mummer4/mummer 
progMUMmerplot="/path/to/MUMmerplot"							#	https://github.com/mummer4/mummer/blob/master/scripts/mummerplot.pl # -- 
progProgressivecactus="/path/to/ProgressiveCactus"					#	https://github.com/glennhickey/progressiveCactus # faster than MuMmer
progCactus2hal=										#	https://github.com/glennhickey/cactus2hal #convert cactus to hal
progHal=										#	https://github.com/ComparativeGenomicsToolkit/hal
progProgressivemauve=""									#	https://darlinglab.org/mauve/user-guide/progressivemauve.html
### Pipeline for Asset								#	https://github.com/VGP/vgp-assembly/tree/master/pipeline/asset
progAsset=										#	https://github.com/dfguan/asset

#### Program Paths for Visualizations of Assemblies
progIGV="/path/to/IGV"									#	https://software.broadinstitute.org/software/igv/download #viz SAM, BAM
progBandage="/path/to/Bandage"								#	https://github.com/rrwick/Bandage # visualization of GFA files
progBandageng=										#	https://github.com/asl/BandageNG # Same. 

#### Program Paths for Scaffolding and Dealing With Repetitive Sequences
progMinimap2=										#	https://github.com/lh3/minimap2 #long reads aligner
progRepeatmasker=									#	https://github.com/rmhubley/RepeatMasker #tool to flag the repeats

#### Program Paths for Genome Annotations
progAugustus=										#	https://github.com/Gaius-Augustus/Augustus # genome annotation of eukaryotes

#### Program Paths for Visualization of Genome Annotations
ProgCircos=

#### Program Files for Pan-genomics
progvg=											#	https://github.com/vgteam/vg #construction of pangenome graphs
progpggb=										#	https://github.com/pangenome/pggb/tree/master #visualization of pan-gen
progodgi=										#	https://github.com/pangenome/odgi #manipulation of pan-genome graphs


#Dependencies

#	python
#	biopython
#	Python 2.7, BOOST libraries and Networkx(version lower than 1.2). FOR SALSA 2.2


#Dependencies for R (= packages)

#install.packages("argparse")
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("GRIDSS")
#install.packages("tidyverse")
#install.packages("readr")
#install.packages("ggtree")
#install.packages("")

#Other dependencies 
            environment
        Basics: An example workflow
        Advanced: Decorating the example workflow
        Additional features
    Short tutorial
    Snakemake Executor Tutorials
    Best practices

Executing workflows

    Command line interface
    Cluster Execution
    Cloud execution
    Job Grouping
    Between workflow caching
    Interoperability
    Monitoring

Defining workflows

    Writing Workflows
    Snakefiles and Rules
    Configuration
    Modularization
    Remote files
    Utils
    Distribution and Reproducibility
    Reports
    Automatically generating unit tests
    Integrating foreign workflow management systems

Project Info

    Citing and Citations
    More Resources
    Frequently Asked Questions
    Contributing
    Credits
    Changelog
    License

    Docs » Snakemake Tutorial » Setup
    Edit on GitHub

Setup
Requirements

To go through this tutorial, you need the following software installed:

    Python ≥3.5

    Snakemake ≥5.24.1

    BWA 0.7

    SAMtools 1.9

    Pysam 0.15

    BCFtools 1.9

    Graphviz 2.42

    Jinja2 2.11

    NetworkX 2.5
# https://maven.apache.org/  for canu
#jdk
# zlib for purgedups 
# runner for purgedups
# python3 for purgedups
# perl 5.6 for QUAST
# GNU make and ar for QUAST
# GCC 4.7 or higher for QUAST



# path to the different tools 
#progVcfutils=vcfutils.pl 	#OK installed with sudo apt -y install *** 
#progBcftools=bcftools 		#OK installed with sudo apt -y install *** 
#progBamtools=bamtools 		#OK installed with sudo apt -y install *** 
#progBedtools=bedtools 		#OK installed with sudo apt -y install *** 
#progPicard=${progPath}/picard-tools-1.119 				#OK 
#progBowtie2=/home/quentin/miniconda3/envs/bowtie2/bin/bowtie2			#OK

#progAbyss=/home/quentin/miniconda3/envs/amos/bin/abyss #KO but 1.1.1-2
#progAmos=/home/quentin/miniconda3/envs/amos/bin #OK
#progNucmer=nucmer #OK
#progGatk=${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar #
#progSoapdenovo2=/home/quentin/miniconda3/envs/amos/bin #OK
#progRemovShortSeq=${progPath}/RemoveShortSeq.jar #OK
#progGetBlocks=${progPath}/GetBlocks.jar #OK
#progFastaToAmos=${progPath}/FastaToAmos.jar #OK
#progWriteSoapConfig=${progPath}/WriteSoapConfig.jar #OK
#progFastaStats=${progPath}/FastaStats.jar #OK
#progSplitSeqLowCov=${progPath}/SplitSeqLowCov.jar #OK




