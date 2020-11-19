## add to configuration file
longitudinal_smk = "/mnt/c/Users/noahs/projects/tagseq-qiime2-snakemake-1/Snakefile-long.smk"
perform_longitudinal = 'yes' # yes or no depending on if you have longitudinal data

## add to snakefile
LONGITUDINAL = config['perform_longitudinal']

run:
       if {LONGITUDINAL} == 'yes':
                include: longitudinal_smk
                print("Performing longitudinal analysis")
       else:
                print("no longitudinal analysis")

## snakefile-long

rule longitudinal_pw_diff:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", T1 = config['start_time'], 
    T2 = config['end_time'], 
    STATE = config['state'], 
    REPLICATE = config['replicate'], 
    COLUMN = config['column'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    T1 = config['start_time']
    T2 = config['end_time']
    STATE = config['state']
    REPLICATE = config['replicate']
    COLUMN = config['column']
  shell:
    """
    qiime longitudinal pairwise-differences \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table}\
      --p-metric {shannon} \
      --p-group-column {delivery} \
      --p-state-column {month} \
      --p-state-1 {0} \
      --p-state-2 {12} \
      --p-individual-id-column {studyid} \
      --p-replicate-handling {random} \
      --o-visualization {pairwise-differences.qzv}
    """

rule longitudinal_pw_dist:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", T1 = config['start_time'], 
    T2 = config['end_time'], 
    STATE = config['state'], 
    REPLICATE = config['replicate'], 
    COLUMN = config['column'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    T1 = config['start_time']
    T2 = config['end_time']
    STATE = config['state']
    REPLICATE = config['replicate']
    COLUMN = config['column']
  shell:
    """
    qiime longitudinal pairwise-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-group-column {delivery} \
      --p-state-column {month} \
      --p-state-1 {0} \
      --p-state-2 {12} \
      --p-individual-id-column {studyid} \
      --p-replicate-handling {random} \
      --o-visualization {pairwise-distances.qzv}
    """

rule longitudinal_me:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", STATE = config['state'], GROUP_COL = config['group_col'], id_col = config['animal_number'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    STATE = config['state']
    GROUP_COL = config['group_col']
    id_col = config['animal_number']
  shell:
    """
    qiime longitudinal linear-mixed-effects \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-metric {shannon} \
      --p-group-columns {delivery,diet,sex} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --o-visualization {linear-mixed-effects.qzv}
    """

rule longitudinal_volitilty:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", STATE = config['state'], COLUMN = config['column'], id_col = config['animal_number'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    STATE = config['state']
    COLUMN = config['column']
    id_col = config['animal_number']
  shell:
    """
    qiime longitudinal volatility \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-default-metric {shannon} \
      --p-default-group-column {delivery} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --o-visualization {volatility.qzv}
    """

rule longitudinal_first_diff:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", STATE = config['state'], id_col = config['animal_number'], REPLICATE config['replicate'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    STATE = config['state']
    id_col = config['animal_number']
    REPLICATE config['replicate']
  shell:
    """
    qiime longitudinal first-differences \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {month} \
      --p-metric {shannon} \
      --p-individual-id-column {studyid} \
      --p-replicate-handling {random} \
      --o-first-differences {shannon-first-differences.qza}
    """

rule longitudinal_first_diff_beta:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", STATE = config['state'], id_col = config['animal_number'], REPLICATE = config['replicate'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    STATE = config['state']
    id_col = config['animal_number']
    REPLICATE = config['replicate']
  shell:
    """
    qiime longitudinal first-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --p-replicate-handling {random} \
      --o-first-distances {first_distances.qza}
  """

rule longitudinal_track:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", METRIC =  config['metric'], STATE = config['state'], id_col = , GROUP_COL = config['group_col'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    METRIC =  config['metric']
    STATE = config['state']
    id_col = config['animal_number']
    GROUP_COL = config['group_col']
  shell:
    """
    qiime longitudinal linear-mixed-effects \
      --m-metadata-file {input.cleaned_metadata} \
      --i-table {input.cleaned_table} \
      --p-metric {Distance} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --p-group-columns {delivery,diet} \
      --o-visualization {first_istances_LME.qzv}
    """

rule longitudinal_first_diff_static:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", COLUMN = config['column'], id_col = config['animal_number'], REPLICATE = config['replicate'], STATE = config['state'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    COLUMN = config['column']
    id_col =
    REPLICATE = config['replicate']
    STATE = config['state']
  shell:
    """
    qiime longitudinal first-distances \
      --i-distance-matrix {input.unweighted_unifrac_mat} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --p-replicate-handling {random} \
      --p-baseline 0 \
      --o-first-distances {first_distances_baseline_0.qza}
    """

rule nmit:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", METHOD = config['method'], id_col = config['animal_number'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    METHOD = config['method']
    id_col = config['animal_number']
  shell:
    """
    qiime longitudinal nmit \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-individual-id-column {studyid} \
      --p-corr-method {pearson} \
      --o-distance-matrix {nmit-dm.qza}
    """

rule longitudinal_beta:
  input:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", METHOD = config['method'], id_col = config['animal_number']),
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", COLUMN = config['column'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    COLUMN = config['column']
  shell:
    """
    qiime diversity beta-group-significance \
      --i-distance-matrix {nmit-dm.qza} \
      --m-metadata-file {input.cleaned_metadata} \
      --m-metadata-column {delivery} \
      --o-visualization {nmit.qzv}
    """

rule longitudinal_pcoa:
  input:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", METHOD = config['method'], id_col = config['animal_number'])
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/output/", METHOD = config['method'], id_col = config['animal_number'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    """
    qiime diversity pcoa \
      --i-distance-matrix {nmit-dm.qza} \
      --o-pcoa {nmit-pc.qza}
    """

rule longitudinal_emperor:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/output/", METHOD = config['method'], id_col = config['animal_number'])
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/output/", METHOD = config['method'], id_col = config['animal_number'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:

  shell:
    """
    qiime emperor plot \
      --i-pcoa {nmit_pc.qza} \
      --m-metadata-file {input.cleaned_metadata} \
      --o-visualization {nmit_emperor.qzv}
    """

rule longitudinal_feature_volitilty:
  input:
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/", 
    STATE = config['state'],
    id_col = ,
    ESTIMATORS = config['estimators'],
    RANDOM_STATE = config['random_state'])
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    STATE = config['state']
    id_col = config['animal_number']
    ESTIMATORS = config['estimators']
    RANDOM_STATE = config['random_state']
    MISSING = config['missing_samples']
    directory(OUTPUTDIR + "/qiime2/asv/longitudinal/volitility") 
  shell:
    """
    qiime longitudinal feature-volatility \
      --i-table {input.cleaned_table} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {month} \
      --p-individual-id-column {studyid} \
      --p-n-estimators {10} \
      --p-random-state {17} \
      --p-missing-samples {} \
      --output-dir {ecam_feat_volatility}
    """

rule maturity_index:
  input:
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  output:
    expand(OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/", STATE = config['state'],
    GROUP_BY = config['group_by'],
    CONTROL = config['control'],
    COLUMN = config['column'] )
  log:
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
      directory(OUTPUTDIR + "/qiime2/asv/logitudinal/maturity")
      STATE = config['state']
      GROUP_BY = config['group_by']
      CONTROL = config['control']
      COLUMN = config['column']
  shell:
    """
    qiime longitudinal maturity-index \
      --i-table {ecam-table.qza} \
      --m-metadata-file {input.cleaned_metadata} \
      --p-state-column {month} \
      --p-group-by {delivery} \
      --p-individual-id-column {studyid} \
      --p-control {Vaginal} \
      --p-test-size 0.4 \
      --p-stratify \
      --p-random-state 1010101 \
      --output-dir maturity
    """