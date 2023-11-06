# Define the variables to perform the scaffolding of contigs after de-novo assembly.smk
 
rule run_greenhill:
    input:
    output:
    params:
    shell:
       """
      greenhill -cph out.hap1.p_ctg.fa out.hap2.p_ctg.fa -p HiFi.fastq -IP1 PE1.fastq PE2.fastq -HiC HiC1.fastq HiC2.fastq
       """
