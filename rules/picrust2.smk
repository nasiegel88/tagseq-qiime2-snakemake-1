# picrust2 analyis
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

rule asv_reftree:
  input:
    seqs = OUTPUTDIR + "/picrust2/dna-sequences.fasta"
  output:
    picrust2tree = OUTPUTDIR + "/picrust2/out.tre"
  log:
    SCRATCH + "/logs/" + PROJ + "_exportTREE_picrust.log"
  conda:
    "../envs/picrust2-env.yaml"
  params:
    directory(OUTPUTDIR + "/picrust2/intermediate/place_seqs")
  shell:
    "place_seqs.py -s {input.seqs} -o {output.picrust2tree} -p 1 --intermediate {params}"

rule hsp:
  input:
    picrust2tree = OUTPUTDIR + "/picrust2/out.tre"
  output:
    marker = OUTPUTDIR + "/picrust2/marker_predicted_and_nsti.tsv.gz",
    EC = OUTPUTDIR + "/picrust2/EC_predicted.tsv.gz"
  log:
    SCRATCH + "/logs/" + PROJ + "_marker_EC_predictions_picrust.log"
  conda:
    "../envs/picrust2-env.yaml"
  shell:
      """
       hsp.py -i 16S -t {input.picrust2tree} -o {output.marker} -p 1 -n
       hsp.py -i EC -t {input.picrust2tree} -o {output.EC} -p 1
      """

rule metagenome:
  input:
    table_tsv = OUTPUTDIR + "/table/feature-table.biom",
    marker = OUTPUTDIR + "/picrust2/marker_predicted_and_nsti.tsv.gz",
    EC = OUTPUTDIR + "/picrust2/EC_predicted.tsv.gz"
  output:
    metagenome_contrib = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_contrib.tsv.gz",
    metagenome_unstrat = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    seqtab_norm = OUTPUTDIR + "/picrust2/EC_metagenome_out/seqtab_norm.tsv.gz",
    weighted_nsti =OUTPUTDIR + "/picrust2/EC_metagenome_out/weighted_nsti.tsv.gz"
  log:
    SCRATCH + "/logs/" + PROJ + "_marker_EC_predictions_picrust.log"
  conda:
    "../envs/picrust2-env.yaml"
  params: 
    directory(OUTPUTDIR + "/picrust2/EC_metagenome_out")
  shell:
    """
    metagenome_pipeline.py -i {input.table_tsv} -m {input.marker} -f {input.EC} \
                             -o {params} --strat_out      
    """                

rule pl_infer:
  input:
    metagenome_unstrat = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
  output:
    path_abun_unstrat = OUTPUTDIR + "/picrust2/pathways_out/path_abun_unstrat.tsv.gz"
  log:
    SCRATCH + "/logs/" + PROJ + "_pre_metagenome_picrust.log"
  conda:
    "../envs/picrust2-env.yaml"
  params:
    directory(OUTPUTDIR + "/picrust2/pathways_out")
  shell:
    """
    pathway_pipeline.py -i {input.metagenome_unstrat} \
                          -o {params} -p 1
    """

rule add_describe:
  input:
    metagenome_unstrat = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    path_abun_unstrat = OUTPUTDIR + "/picrust2/pathways_out/path_abun_unstrat.tsv.gz"
  output:
    marker_describe = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz",
    path_abun_unstrat_describe = OUTPUTDIR + "/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz"
  log:
    SCRATCH + "/logs/" + PROJ + "_metagenome_description_picrust.log"
  conda:
    "../envs/picrust2-env.yaml"
  shell:
      """
        add_descriptions.py -i {input.metagenome_unstrat} -m EC \
                              -o {output.marker_describe}
        add_descriptions.py -i {input.path_abun_unstrat} -m METACYC \
                              -o {output.path_abun_unstrat_describe}
      """