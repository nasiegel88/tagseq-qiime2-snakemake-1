# Longitudinal analysis
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu


rule rel_freq:
  input:
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    relative_frequency = OUTPUTDIR  + "/" + PROJ + "-relative_frequency_table.qza"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime feature-table relative-frequency \
      --i-table {input.cleaned_table} \
      --o-relative-frequency-table {output.relative_frequency}
    """

rule longitudinal_pw_diff:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
     pw_diff = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "pairwise-differences.qzv"
  conda:
    "../envs/qiime2-2020.8.yaml"
  log:
    SCRATCH + "/logs/" + PROJ + "pwd1.log"
  shell:
    """
    qiime longitudinal pairwise-differences \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table}\
      --p-metric {METACATEGORY} \
      --p-state-column {STATE} \
      --p-state-1 0 \
      --p-state-2 6 \
      --p-group-column {METACATEGORY}
      --p-individual-id-column {ID} \
      --p-replicate-handling {REPLICATE} \
      --o-visualization {output.pw_diff}
    """

rule longitudinal_pw_dist:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    pw_dist = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_pairwise-distances.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "pwd.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal pairwise-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-group-column {METACATEGORY} \
      --p-state-column {STATE} \
      --p-state-1 0 \
      --p-state-2 6 \
      --p-individual-id-column {ID} \
      --p-replicate-handling {REPLICATE} \
      --o-visualization {output.pw_dist}
    """

rule longitudinal_me:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    lme = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_linear-mixed-effects.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "lme1.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal linear-mixed-effects \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-metric {METACATEGORY} \
      --p-group-columns {GROUPMETACAT} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --o-visualization {output.lme}
    """

rule longitudinal_volitilty:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    vol = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_volatility.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "vol.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal volatility \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-default-metric {METACATEGORY} \
      --p-default-group-column {config[metadata_category]} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --o-visualization {output.vol}
    """

rule longitudinal_first_diff:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    shannon_fd = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_shannon-first-differences.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "longitudinal-dif.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal first-differences \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-metric {METACATEGORY} \
      --p-individual-id-column {ID} \
      --p-replicate-handling {REPLICATE} \
      --o-first-differences {output.shannon_fd}
    """

rule longitudinal_first_diff_beta:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/core-metrics-results/" + PROJ +  "-unweighted_unifrac_distance_matrix.qza"
  output:
    first_diff = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_first_distances.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "first-diff-beta.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal first-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --p-replicate-handling {REPLICATE} \
      --o-first-distances {output.first_diff}
  """

rule longitudinal_track:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    first_dist = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_first_distances_LME.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "lme.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal linear-mixed-effects \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-metric {METACATEGORY} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --p-group-columns {GROUPMETACAT} \
      --o-visualization {output.first_dist}
    """

rule longitudinal_first_diff_static:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    baseline_fd = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_first_distances_baseline_0.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "first-dif.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal first-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --p-replicate-handling {REPLICATE} \
      --p-baseline 0 \
      --o-first-distances {output.baseline_fd}
    """

rule nmit:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    relative_frequency = OUTPUTDIR  + "/" + PROJ + "-relative_frequency_table.qza"
  output:
    nmit = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "nmit.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime longitudinal nmit \
      --i-table {input.relative_frequency} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-individual-id-column {ID} \
      --p-corr-method {ALPHASTATISTIC} \
      --o-distance-matrix {output.nmit}
    """

rule longitudinal_beta:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    nmit = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  output:
    nmit_viz = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_nmit.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "longitudinal-beta.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity beta-group-significance \
      --i-distance-matrix {input.nmit} \
      --m-metadata-file {input.cleaned_metadata} \
      --m-metadata-column {config[metadata_category]} \
      --o-visualization {output.nmit_viz}
    """

rule longitudinal_pcoa:
  input:
    nmit = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  output:
    nmit_pcoa = OUTPUTDIR + "/longitudinal/output/" + PROJ + ANALYSIS + "_nmit-pc.qza"
  log:
    SCRATCH + "/logs/" + PROJ + "longitudinal-pcoa.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  shell:
    """
    qiime diversity pcoa \
      --i-distance-matrix {input.nmit} \
      --o-pcoa {output.nmit_pcoa}
    """

rule longitudinal_emperor:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    nmit_pcoa = OUTPUTDIR + "/longitudinal/output/" + PROJ + ANALYSIS + "_nmit-pc.qza"
  output:
    nmit_emp = OUTPUTDIR + "/longitudinal/output/" + PROJ + ANALYSIS + "_nmit_emperor.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "longitudinal-emp.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  params:

  shell:
    """
    qiime emperor plot \
      --i-pcoa {input} \
      --m-metadata-file {input.cleaned_metadata} \
      --o-visualization {output.nmit_emp}
    """

rule longitudinal_feature_volitilty:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    vol_table = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "filtered_table.qza",
    vol_est = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "sample_estimator.qza",
    vol_imp = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "feature_importance.qza",
    vol_acc = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "accuracy_results.qzv",
    vol_plot = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "volatility_plot.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "longitudinal-vol.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  params:
    directory(OUTPUTDIR + "/longitudinal/mld_feat_volatility") 
  shell:
    """
    qiime longitudinal feature-volatility \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --p-missing-samples {MISSING} \
      --o-filtered-table  {output.vol_table} \
      --o-feature-importance {output.vol_imp} \
      --o-volatility-plot {output.vol_plot} \
      --o-accuracy-results {output.vol_acc} \
      --o-sample-estimator {output.vol_est}
    """

rule maturity_index:
  input:
    cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    maturity_maz =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "maz_scores.qza",
    maturity_emp =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "sample_estimator.qza",
    maturity_imp =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "feature_importance.qza",
    maturity_pred = OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "predictions.qza",
    maturity_acc =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "accuracy_results.qzv",
    maturity_vol =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "volatility_plots.qzv",
    maturity_clus = OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "clustermap.qzv",
    maturity_mod =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "model_summary.qzv"
  log:
    SCRATCH + "/logs/" + PROJ + "maturity.log"
  conda:
    "../envs/qiime2-2020.8.yaml"
  params:
      directory(OUTPUTDIR + "/logitudinal/maturity")
  shell:
    """
    qiime longitudinal maturity-index \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-group-by {METACATEGORY} \
      --p-individual-id-column {ID} \
      --p-control {CONTROL} \
      --p-test-size 0.4 \
      --p-stratify \
      --p-random-state 1010101 \
      --o-feature-importance {output.maturity_imp} \
      --o-predictions {output.maturity_pred} \
      --o-model-summary {output.maturity_mod} \
      --o-accuracy-results {output.} \
      --o-maz-scores {output.maturity_maz} \
      --o-clustermap {output.maturity_clus} \
      --o-volatility-plots {output.maturity_vol}
    """
