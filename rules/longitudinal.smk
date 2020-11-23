# Longitudinal analysis
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

## still in development ###

rule longitudinal_pw_diff:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
     pw_diff = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "pairwise-differences.qzv"
  conda:
    "../envs/qiime2-2019.10.yaml"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "pwd1.log"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    pw_dist = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_pairwise-distances.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "pwd.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    lme = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_linear-mixed-effects.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "lme1.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    vol = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_volatility.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "vol.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    shannon_fd = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_shannon-first-differences.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "longitudinal-dif.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ +  "-unweighted_unifrac_distance_matrix.qza"
  output:
    first_diff = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "first-diff-beta.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    first_dist = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances_LME.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "lme.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    baseline_fd = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances_baseline_0.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "first-dif.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    nmit = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "nmit.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
    qiime longitudinal nmit \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-individual-id-column {ID} \
      --p-corr-method {ALPHASTATISTIC} \
      --o-distance-matrix {output.nmit}
    """

rule longitudinal_beta:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    nmit = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  output:
    nmit_viz = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "longitudinal-beta.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    nmit = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"
  output:
    nmit_pcoa = OUTPUTDIR + "/qiime2/asv/longitudinal/output/" + PROJ + ANALYSIS + "_nmit-pc.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "longitudinal-pcoa.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
    qiime diversity pcoa \
      --i-distance-matrix {input.nmit} \
      --o-pcoa {output.nmit_pcoa}
    """

rule longitudinal_emperor:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    nmit_pcoa = OUTPUTDIR + "/qiime2/asv/longitudinal/output/" + PROJ + ANALYSIS + "_nmit-pc.qza"
  output:
    nmit_emp = OUTPUTDIR + "/qiime2/asv/longitudinal/output/" + PROJ + ANALYSIS + "_nmit_emperor.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "longitudinal-emp.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    vol_table = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "filtered_table.qza",
    vol_est = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "sample_estimator.qza",
    vol_imp = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "feature_importance.qza",
    vol_acc = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "accuracy_results.qzv",
    vol_plot = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "volatility_plot.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "longitudinal-vol.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility") 
  shell:
    """
    qiime longitudinal feature-volatility \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {STATE} \
      --p-individual-id-column {ID} \
      --p-missing-samples {MISSING} \
      --output-dir {params}
    """

rule maturity_index:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    maturity_maz =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "maz_scores.qza",
    maturity_emp =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "sample_estimator.qza",
    maturity_imp =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "feature_importance.qza",
    maturity_pred = OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "predictions.qza",
    maturity_acc =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "accuracy_results.qzv",
    maturity_vol =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "volatility_plots.qzv",
    maturity_clus = OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "clustermap.qzv",
    maturity_mod =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "model_summary.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "maturity.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
  params:
      directory(OUTPUTDIR + "/qiime2/asv/logitudinal/maturity")
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
      --output-dir {params}
    """
