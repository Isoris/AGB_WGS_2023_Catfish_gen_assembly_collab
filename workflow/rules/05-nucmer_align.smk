# Define rules to perform one-to-one whole genome alignments

rule mummer_dna_diff:
   input:
      query_genome = path_out_prefix + "02-ASSEMBLIES/{sample}"
      reference_genome = path_data_prefix + "01-REFERENCES/{reference}}
   output:
      $out.dnadiff.1delta"
      $out.dnadiff.plot.1.png
   conda:
   "../envs/alignments.yaml"
   script:


   