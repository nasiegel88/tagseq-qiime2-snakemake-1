# Phylogeny analyis
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu 

rule alignment:
  input:
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza"
  output:
    masked_align = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-masked_alignment.qza",
    rooted_tree = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-rooted_tree.qza",
    align = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-alignment.qza",
    tree = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-tree.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-aligned.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input.rep} \
        --o-alignment {output.align} \
        --o-masked-alignment {output.masked_align} \
        --o-tree {output.tree} \
        --o-rooted-tree {output.rooted_tree}
    """

rule core_metrics:
  input:
    rooted_tree = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-rooted_tree.qza",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    evenness_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-evenness_vector.qza",
    faith_pd_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-faith_pd_vector.qza",
    observed_otus = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
    rarefied_table = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-rarefied_table.qza",
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon_vector.qza",
    bray_curtis = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_distance_matrix.qza",   
    bray_curtis_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_pcoa_results.qza",
    jaccard_distance = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_distance_matrix.qza",
    jaccard_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_pcoa_results.qza",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza",
    unweighted_unifrac_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ +"-unweighted_unifrac_pcoa_results.qza",
    weighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_distance_matrix.qza",
    weighted_unifrac_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_pcoa_results.qza",
    bray_curtis_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_emperor.qzv",
    jaccard_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_emperor.qzv",
    unweighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_emperor.qzv",
    weighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_emperor.qzv",
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-core-metrics.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny {input.rooted_tree} \
        --i-table {input.cleaned_table} \
        --p-sampling-depth {config[sampling-depth]} \
        --m-metadata-file {input.cleaned_metadata} \
        --o-rarefied-table {output.rarefied_table} \
        --o-faith-pd-vector {output.faith_pd_vector} \
        --o-observed-features-vector {output.observed_otus} \
        --o-shannon-vector {output.shannon_vector} \
        --o-evenness-vector {output.evenness_vector} \
        --o-unweighted-unifrac-distance-matrix {output.unweighted_unifrac_mat} \
        --o-weighted-unifrac-distance-matrix {output.weighted_unifrac_mat} \
        --o-jaccard-distance-matrix {output.jaccard_distance} \
        --o-bray-curtis-distance-matrix {output.bray_curtis} \
        --o-unweighted-unifrac-pcoa-results {output.unweighted_unifrac_pcoa} \
        --o-weighted-unifrac-pcoa-results {output.weighted_unifrac_pcoa} \
        --o-jaccard-pcoa-results {output.jaccard_pcoa} \
        --o-bray-curtis-pcoa-results {output.bray_curtis_pcoa} \
        --o-unweighted-unifrac-emperor {output.unweighted_unifrac} \
        --o-weighted-unifrac-emperor {output.weighted_unifrac} \
        --o-jaccard-emperor {output.jaccard_emperor} \
        --o-bray-curtis-emperor {output.bray_curtis_emperor}
    """

rule richness:
  input:
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ +  "-shannon_vector.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    shannon_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ +  "-shannon-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity alpha-group-significance \
        --i-alpha-diversity {input.shannon_vector} \
        --m-metadata-file {input.cleaned_metadata} \
        --o-visualization {output.shannon_signif}
    """

rule richcorr:
  input:
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon_vector.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    shannon_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance-association.qzv",
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-shannon-significance-association.qzv"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity alpha-correlation \
        --i-alpha-diversity {input.shannon_vector} \
        --m-metadata-file {input.cleaned_metadata} \
        --p-method {config[alpha-div-p-method]} \
        --o-visualization {output.shannon_correl}
    """

rule asv_signif:
  input:
    observed_asv = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    observed_asv_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-observed_otus-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity alpha-group-significance \
        --i-alpha-diversity {input.observed_asv} \
        --m-metadata-file {input.cleaned_metadata} \
        --o-visualization {output.observed_asv_signif}
    """

rule asv_corr:
  input:
    observed_asv = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    observed_asv_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed-otus-significance-association.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-observed_otus-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity alpha-correlation \
        --i-alpha-diversity {input.observed_asv} \
        --m-metadata-file {input.cleaned_metadata} \
        --p-method {config[alpha-div-p-method]} \
        --o-visualization {output.observed_asv_correl}
    """

rule evenness:
  input:
    bray_curtis = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_distance_matrix.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    bray_curtis_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-bray-curtis-group-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-bray-curtis-group-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity beta-group-significance \
        --i-distance-matrix {input.bray_curtis} \
        --m-metadata-file {input.cleaned_metadata} \
        --m-metadata-column {config[metadata_category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.bray_curtis_signif} \
        --p-no-pairwise
    """

rule unifrac:
  input:
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    unweighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-unweighted-unifrac-group-site-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-unweighted-unifrac-group-site-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity beta-group-significance \
        --i-distance-matrix {input.unweighted_unifrac_mat} \
        --m-metadata-file {input.cleaned_metadata} \
        --m-metadata-column {config[metadata_category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.unweighted_unifrac_viz} \
        --p-no-pairwise
    """

rule weighted_unifrac:
  input:
    weighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ +  "-weighted_unifrac_distance_matrix.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    weighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-weighted-unifrac-group-site-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-weighted-unifrac-group-site-significance.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity beta-group-significance \
        --i-distance-matrix {input.weighted_unifrac_mat} \
        --m-metadata-file {input.cleaned_metadata} \
        --m-metadata-column {config[metadata_category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.weighted_unifrac_viz} \
        --p-no-pairwise
    """

rule barplot:
  input:
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza",
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    barplots = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-taxa-bar-plots.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-taxa-bar-plots.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime taxa barplot \
        --i-table {input.cleaned_table} \
        --i-taxonomy {input.sklearn} \
        --m-metadata-file {input.cleaned_metadata} \
        --o-visualization {output.barplots}
    """