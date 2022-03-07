# dada2
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

rule dada2:
  input:
    q2_primerRM = SCRATCH + "/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    table = OUTPUTDIR + "/" + PROJ + "-asv-table.qza",
    rep = OUTPUTDIR + "/" + PROJ + "-rep-seqs.qza",
    stats = OUTPUTDIR + "/" + PROJ + "-stats-dada2.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_dada2_q2.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
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

rule filter_seqs:
  input:
    rep = OUTPUTDIR + "/" + PROJ + "-rep-seqs.qza",
    sklearn = OUTPUTDIR + "/" + PROJ + "-tax_sklearn.qza"
  output:
    filtered_rep = OUTPUTDIR + "/" + PROJ + "-filtered-rep-seqs.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_filtered-seqs.log"
  conda:
    "../envs/qiime2-2020.8.yaml" 
  shell:
    """
    qiime taxa filter-seqs \
      --i-sequences {input.rep} \
      --i-taxonomy {input.sklearn} \
      --p-include p__ \
      --p-exclude {EXCLUDESEQS} \
      --o-filtered-sequences {output.filtered_rep}
    """

rule filter_table:
  input:
    sklearn = OUTPUTDIR + "/" + PROJ + "-tax_sklearn.qza",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    filtered_table = OUTPUTDIR + "/" + PROJ + "-filtered-asv-table.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_filtered-table.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime taxa filter-table \
      --i-table {input.cleaned_table} \
      --i-taxonomy {input.sklearn} \
      --p-mode exact \
      --p-exclude {EXCLUDESEQS} \
      --o-filtered-table {output.filtered_table}
    """

rule filter_features:
  input:
    filtered_table = OUTPUTDIR + "/" + PROJ + "-filtered-asv-table.qza",
    filtered_sklearn = OUTPUTDIR + "/" + PROJ + "-filtered_tax_sklearn.qza"
  output:
    id_filtered_table = OUTPUTDIR + "/" + PROJ + "-id_filtered-asv-table.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_filtered-features.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime feature-table filter-features \
      --i-table {input.filtered_table} \
      --m-metadata-file {input.filtered_sklearn} \
      --o-filtered-table {output.id_filtered_table}
    """

  
rule metadata:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv"
  output:
    tabulated_cleaned_metadata = HOME + "/tabulated-sample-metadata.qzv"
  log:
    SCRATCH + "/logs/tabulated-metadata.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime metadata tabulate \
      --m-input-file {input.cleaned_metadata} \
      --o-visualization {output.tabulated_cleaned_metadata}
    """

rule drop_blanks:
  input:
    table = OUTPUTDIR + "/" + PROJ + "-asv-table.qza"
  output:
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza",
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv"
  log:
    SCRATCH + "/logs/" + PROJ + "-remove-blanks.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  params:
    metadata  = config['metadata'],
    table_blanks = config['table_blanks']
  shell:
    """
    declare -a arr=("{config[remove_blanks]}")
    if [ "${{arr[@]}}" == yes ]; then
    sed -e '/BLANK*/d' -e '/NS.B6.47455F/d' {params.metadata} > {output.cleaned_metadata} # NS.B6.47455F was used as an additional control
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
    stats = OUTPUTDIR + "/" + PROJ + "-stats-dada2.qza"
  output:
    stats_viz = OUTPUTDIR + "/" + PROJ + "-stats-dada2.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "_dada2-stats_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  shell:
    """
      qiime metadata tabulate \
          --m-input-file {input.stats} \
          --o-visualization {output.stats_viz}
    """

rule assign_tax:
  input:
    rep = OUTPUTDIR + "/" + PROJ + "-rep-seqs.qza",
    db_classified = CLASSIFIER
  output:
    sklearn = OUTPUTDIR + "/" + PROJ + "-tax_sklearn.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_sklearn_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime feature-classifier classify-sklearn \
       --i-classifier {input.db_classified} \
       --i-reads {input.rep} \
       --o-classification {output.sklearn}
    """

rule filtered_assign_tax:
  input:
    filtered_rep = OUTPUTDIR + "/" + PROJ + "-filtered-rep-seqs.qza",
    db_classified = CLASSIFIER
  output:
    filtered_sklearn = OUTPUTDIR + "/" + PROJ + "-filtered_tax_sklearn.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "_sklearn_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime feature-classifier classify-sklearn \
       --i-classifier {input.db_classified} \
       --i-reads {input.filtered_rep} \
       --o-classification {output.filtered_sklearn}
    """

rule gen_table:
  input:
    id_filtered_table = OUTPUTDIR + "/" + PROJ + "-id_filtered-asv-table.qza"
  output:
    table_tsv = OUTPUTDIR + "/table/feature-table.biom"
  log:
    SCRATCH + "/logs/" + PROJ + "_exportBIOM_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  params:
    directory(OUTPUTDIR + "/table")
  shell:
    "qiime tools export --input-path {input.id_filtered_table} --output-path {params}"

rule convert:
  input:
    OUTPUTDIR + "/table/feature-table.biom"
  output:
    OUTPUTDIR + "/" + PROJ + "-asv-table.tsv"
  log:
    SCRATCH + "/logs/" + PROJ + "-exportTSV_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  shell:
    "biom convert -i {input} -o {output} --to-tsv"

rule gen_tax:
  input:
    filtered_sklearn = OUTPUTDIR + "/" + PROJ + "-filtered_tax_sklearn.qza"
  output:
    table_tax = OUTPUTDIR + "/tax_assigned/taxonomy.tsv"
  log:
    SCRATCH + "/logs/" + PROJ + "-exportTAXTSV_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  params:
    directory(OUTPUTDIR + "/tax_assigned")
  shell:
    "qiime tools export --input-path {input.filtered_sklearn} --output-path {params}"

rule gen_seqs:
  input:
    filtered_rep = OUTPUTDIR + "/" + PROJ + "-filtered-rep-seqs.qza",
  output:
    seqs = OUTPUTDIR + "/picrust2/dna-sequences.fasta"
  log:
    SCRATCH + "/logs/" + PROJ + "_exportTAXTSV_q2.log"
  conda:
     "../envs/qiime2-2020.8.yaml"
  params:
    directory(OUTPUTDIR + "/picrust2")
  shell:
    "qiime tools export --input-path {input.filtered_rep} --output-path {params}"