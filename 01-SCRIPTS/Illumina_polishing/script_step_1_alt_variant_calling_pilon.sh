#!/bin/bash

#https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish

# Source the config.sh to get variables
config_path="/path/to/config.sh"
source ${config_path}

# Log File
log=${path_prefix_step_1_pilon}/Illumina_polishing_pilon.log

# Create directories for this step
mkdir -p ${path_prefix_step1_trimmed}Your_step_here

# Initialize the array to hold reference names
declare -a shortNamesRef=()

# Go to the directory containing the consensus files and check for indexes
cd ${path_prefix_step_1_pacbio}

# Loop over each fasta file to build index
for fasta in *.fasta; do
  ref_name=$(basename "$fasta" .fasta)
  
  # Assign to variable Consensusfasta
  Consensusfasta=${ref_name}

  # Log and Build index
  if [ ! -f "${Consensusfasta}.1.bt2" ]; then
    echo "Building index for ${Consensusfasta}..." >> $log
    ${progBowtie2}-build ${fasta} ${Consensusfasta}
  fi

  # Append reference name to the array
  shortNamesRef+=("$Consensusfasta")
done

cd -  # Move back to the original directory


# Loop starts here
for y in ${!shortNamesRef[*]}
do  
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mappedAll[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_01_all.sorted.bam
    mapped[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_02_sorted.bam
    unmapped[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_03_unmapped.sorted.bam
    mappedFiltered[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_04_sorted.filtered.bam
    bowtieFailPair[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_05_failPair.fastq
    bowtieFailUnPair1[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_06_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_07_failUnPairR2.fastq  
    mappedAllpair[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_08_allreads.bam
    mappedAllpairfixmate[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}
    mappedAllpairfixmatesorted[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}
    FinalBAM[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_11_Final
    Consensusvcf[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_12_vcf
    Consensusfq[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_13_fastq
    Consensusfasta[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_14_cns
    Consensusfinalfasta[i]=${path_prefix_step_1_pilon}/${shortNames[i]}_${shortNamesRef[y]}_15_cns
    BT_INDEX[i]=${path_prefix_step_1_pilon}/bowtie2_indexes_tmp/${shortNames[i]}_${shortNamesRef[y]}
    newconsensus=${Consensusfasta[i]}
    count=0  
    (

	
# Assign variables for the location of the Bowtie2 index, the input FASTQ files, and the output BAM file
	# index reference file  
	echo "run reference mapping to new consensus before pilon..."  
	echo "run reference mapping to new consensus before pilon..." >> $log
	mkdir -p bowtie2_indexes_tmp
	
	${progBowtie2}-build ${Consensusfasta[i]}.fasta  ${BT_INDEX[i]} 
	
	OUTPUT_BAM=output.improved.bam
	
	${progBowtie2} --sensitive -q --phred33 --no-unal --no-mixed --no-discordant -p ${CPUS} -X ${insHigh[i]} -x ${BT_INDEX[i]} -1 ${reads1trimmed[i]} -2 ${reads2trimmed[i]} | ${progSamtools} view -bhS -@ ${CPUthreads} - | ${progSamtools} sort - -T ${shortNames[i]} -n -o ${OUTPUT_BAM%.bam}_${shortNames[i]}_${shortNamesRef[y]}_improved.sorted.bam -@ ${CPUthreads} >> $log
	
 #filter unmapped reads    
      ${progSamtools} view -@ ${CPUthreads} -bh -F 4 ${mappedAll[i]} > ${mapped[i]}
      ${progSamtools} index ${mapped[i]}
      
      #get unmapped reads
      ${progSamtools} view -@ ${CPUthreads} -bh -f 4 ${mappedAll[i]} > ${unmapped[i]}
      ${progSamtools} view -@ ${CPUthreads} -bh -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard}/SamToFastq.jar INPUT=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${bowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -@ ${CPUthreads} -bh -F 8 -f 64 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair1[i]}  
      ${progSamtools} view -@ ${CPUthreads} -bh -F 8 -f 128 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair2[i]}
    
      ${progBamtools} stats -in ${mappedAll[i]} >> $log
      echo "--> ${mappedAll[i]}" >> $log
    
      #filter for mapping quality >=10    
      ${progSamtools} view -@ ${CPUthreads} -bh -q 10 ${mapped[i]} > ${mappedFiltered[i]}
      ${progBamtools} stats -in ${mappedFiltered[i]} >> $log
      echo "--> ${mappedFiltered[i]}" >> $log
      
      #Sort by read names, add the flags containing informations about the mates for each read in pair and remove PCR Duplicates 
      ${progSamtools} sort -@ ${CPUthreads} -n -o ${mappedAllpair[i]}  -O BAM ${mappedAll[i]}  
      ${progSamtools} fixmate -@ ${CPUthreads} -m ${mappedAllpair[i]} ${mappedAllpairfixmate[i]}_09_fixmate.bam
      ${progSamtools} sort -o ${mappedAllpairfixmatesorted[i]}_10_fixmate_sorted.bam ${mappedAllpairfixmate[i]}_09_fixmate.bam 
      ${progSamtools} markdup -@ ${CPUthreads} -r -s ${mappedAllpairfixmatesorted[i]}_10_fixmate_sorted.bam ${FinalBAM[i]}.bam
      ${progSamtools} index ${FinalBAM[i]}.bam

	${progPilon} --genome ${Consensusfasta[i]}.fasta --bam ${FinalBAM[i]}.bam --output ${Consensusfinalfasta[i]} --verbose --fix all --vcf --tracks >> $log
	# convert fasta to fastq
	java -jar ${progPicard}/NormalizeFasta.jar INPUT=${Consensusfinalfasta[i]}.fasta LINE_LENGTH=80 OUTPUT=${Consensusfinalfasta[i]}_fixedlines.fasta
	#${progSeqtk} seq -aQ64 ${Consensusfinalfasta[i]}.fasta > ${Consensusfq[i]}.fastq
	# convert fastq to fna
	#${progSeqtk} seq -A ${Consensusfq[i]}.fastq > ${Consensusfinalfasta[i]}.fna

	cp ${Consensusfinalfasta[i]}.fasta ${shortNames[i]}_pilon.fasta
	
	mkdir -p OLD/${shortNames[i]}_${shortNamesRef[y]}
	mv ${shortNames[i]}_${shortNamesRef[y]} OLD/${shortNames[i]}_${shortNamesRef[y]}
	) 
	
	 )
  done
done

  done
done
