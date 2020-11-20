# dada2
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

rule dada2:
  input:
    q2_primerRM = SCRATCH + "/qiime2/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza",
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    stats = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_dada2_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input.q2_primerRM} \
        --p-trunc-q {config[truncation_err]} \
        --p-trunc-len-f {config[truncation_len-f]} \
        --p-trunc-len-r {config[truncation_len-r]} \
        --p-max-ee-f {config[quality_err]} \
        --p-max-ee-r {config[quality_err]} \
        --p-n-reads-learn {config[training]} \
        --p-chimera-method {config[chimera]} \
        --o-table {output.table} \
        --o-representative-sequences {output.rep} \
        --o-denoising-stats {output.stats}
    """

rule drop_blanks:
  input:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza"
  output:
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-remove-blanks.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
    metadata  = config['metadata'],
    table_blanks = config['table_blanks'],
    meta_blanks = config['metadata_blanks']
  shell:
    """
    declare -a arr=("{config[remove_blanks]}")
    if [ "${{arr[@]}}" == yes ]; then
    fgrep -v {params.meta_blanks} {params.metadata} > {output.cleaned_metadata}
      echo "All blanks dropped"
      qiime feature-table filter-samples \
      --i-table {input.table} \
      --m-metadata-file {params.metadata} \
      --p-exclude-ids TRUE  \
      --p-where "{params.table_blanks}"  \
      --o-filtered-table {output.cleaned_table}
    elif [ "${{arr[@]}}" == 'no' ]; then
      cp {params.metadata} {output.cleaned_metadata}
      echo "no blanks dropped"
      qiime feature-table filter-samples \
      --i-table {input.table} \
      --m-metadata-file {params.metadata} \
      --p-exclude-ids FALSE \
      --o-filtered-table {output.cleaned_table}
    fi
    """

rule dada2_stats:
  input:
    stats = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qza"
  output:
    stats_viz = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_dada2-stats_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
      qiime metadata tabulate \
          --m-input-file {input.stats} \
          --o-visualization {output.stats_viz}
    """

rule assign_tax:
  input:
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    db_classified = DB_classifier
  output:
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_sklearn_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
    qiime feature-classifier classify-sklearn \
       --i-classifier {input.db_classified} \
       --i-reads {input.rep} \
       --o-classification {output.sklearn}
    """

rule gen_table:
  input:
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    table_tsv = OUTPUTDIR + "/qiime2/asv/table/feature-table.biom"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportBIOM_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/table")
  shell:
    "qiime tools export --input-path {input.cleaned_table} --output-path {params}"

rule convert:
  input:
    OUTPUTDIR + "/qiime2/asv/table/feature-table.biom"
  output:
    OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-exportTSV_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    "biom convert -i {input} -o {output} --to-tsv"

rule gen_tax:
  input:
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza"
  output:
    table_tax = OUTPUTDIR + "/qiime2/asv/tax_assigned/taxonomy.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-exportTAXTSV_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/tax_assigned")
  shell:
    "qiime tools export --input-path {input.sklearn} --output-path {params}"

rule gen_seqs:
  input:
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
  output:
    seqs = OUTPUTDIR + "/qiime2/asv/picrust2/dna-sequences.fasta"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportTAXTSV_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/picrust2")
  shell:
    "qiime tools export --input-path {input.rep} --output-path {params}"