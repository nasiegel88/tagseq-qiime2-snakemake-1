# Snakemake file - input raw fastq reads to generate ASVs
## Last updated 8-18-2020 NS
configfile: "config.yaml"

import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#
PROJ = config["proj_name"]
ANALYSIS = config["analysis_type"]
STATISTIC = config["stat_type"]
INPUTDIR = config["raw-data"]
SCRATCH = config["scratch"]
OUTPUTDIR = config["outputDIR"]
HOME = config["home"]
LONGITUDINAL = config['perform_longitudinal']
DROP_BLANKS = config['remove_blanks']
EXCLUDESEQS =config['seqs-to-exclude']

SUF = config["suffix"]
R1_SUF = str(config["r1_suf"])
R2_SUF = str(config["r2_suf"])

# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(INPUTDIR + "/{sample}_{num}" + SUF)
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

# Uncomment to diagnose if snakemake is reading in wildcard variables properly
#print(SAMPLE_SET)
#print(SET_NUMS)
#print(SUF)
#print(R1_SUF)

# Create final manifest file for qiime2
MANIFEST = pd.read_csv(config["manifest"]) #Import manifest
MANIFEST['filename'] = MANIFEST['absolute-filepath'].str.split('/').str[-1] #add new column with only file name
MANIFEST.rename(columns = {'absolute-filepath':'rawreads-filepath'}, inplace = True)
PATH_TRIMMED = "trimmed" # name of directory with trimmed reads
NEWPATH = os.path.join(SCRATCH, PATH_TRIMMED)
MANIFEST['filename'] = MANIFEST['filename'].str.replace(SUF, "_trim.fastq.gz")
MANIFEST['absolute-filepath'] = NEWPATH+ "/" + MANIFEST['filename']    
MANIFEST[['sample-id','absolute-filepath','direction']].set_index('sample-id').to_csv('manifest-trimmed.txt')
MANIFEST = config["manifest"]

# Database information to assign taxonomy
CLASSIFIER = config['classifier']
DB_classifier = config["database_classified"]

# Phylogeny
META = config["metadata"]
METACATEGORY = config["metadata_category"]
ALPHASTATISTIC = config["alpha-div-p-method"]
BETASTATISTIC = config["beta-div-p-method"]
PERMNUMBER = config["permutations"]

# Longitudinal
ID = config['studyid']          
STATE = config['state']      
GROUPMETACAT = config['group_metadata_category']   
REPLICATE = config['replicate']          
CONTROL = config['control'] 
MISSING = config['missing_samples']

#----DEFINE RULES----#

# fastqc output before trimming
raw_html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
raw_zip = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html",
# Trimmed data output
trimmedData = expand("{scratch}/trimmed/{sample}_{num}_trim.fastq.gz", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS), 
trim_html = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name,

# QIIME2 outputs
q2_import = SCRATCH + "/" + PROJ + "-PE-demux.qza",
q2_primerRM = SCRATCH + "/" + PROJ + "-PE-demux-noprimer.qza",

#vizualization stats
raw = OUTPUTDIR + "/viz/" + PROJ + "-PE-demux.qzv",
primer = OUTPUTDIR + "/viz/" + PROJ + "-PE-demux-noprimer.qzv",

#ASV outputs:
filtered_rep = OUTPUTDIR + "/" + PROJ + "-filtered-rep-seqs.qza",
id_filtered_table = OUTPUTDIR + "/" + PROJ + "-id_filtered-asv-table.qza",
filtered_table = OUTPUTDIR + "/" + PROJ + "-filtered-asv-table.qza",
tabulated_cleaned_metadata = HOME + "/tabulated-sample-metadata.qzv",
table = OUTPUTDIR + "/" + PROJ + "-asv-table.qza",
cleaned_table = OUTPUTDIR + "/" + PROJ + "-no_blanks-asv-table.qza",
cleaned_metadata = HOME + "/noblank-sample-metadata.tsv",
rep = OUTPUTDIR + "/" + PROJ + "-rep-seqs.qza",
stats = OUTPUTDIR + "/" + PROJ + "-stats-dada2.qza",
stats_viz = OUTPUTDIR + "/" + PROJ + "-stats-dada2.qzv",
sklearn = OUTPUTDIR + "/" + PROJ + "-tax_sklearn.qza",
filtered_sklearn = OUTPUTDIR + "/" + PROJ + "-filtered_tax_sklearn.qza",
biom = OUTPUTDIR + "/table/feature-table.biom",
table_tsv = OUTPUTDIR + "/" + PROJ + "-asv-table.tsv",
table_tax = OUTPUTDIR + "/tax_assigned/taxonomy.tsv",
seqs = OUTPUTDIR + "/picrust2/dna-sequences.fasta",
masked_align = OUTPUTDIR + "/" + "mafft-fasttree-output/" + PROJ + "-masked_alignment.qza",
rooted_tree = OUTPUTDIR + "/" + "mafft-fasttree-output/" + PROJ + "-rooted_tree.qza",
align = OUTPUTDIR + "/" + "mafft-fasttree-output/" + PROJ + "-alignment.qza",
tree = OUTPUTDIR + "/" + "mafft-fasttree-output/" + PROJ + "-tree.qza",

# core metrics
evenness_vector = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-evenness_vector.qza",
faith_pd_vector = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-faith_pd_vector.qza",
observed_asv = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
rarefied_table = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-rarefied_table.qza",
shannon_vector = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-shannon_vector.qza",
bray_curtis = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-bray_curtis_distance_matrix.qza",   
bray_curtis_pcoa = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-bray_curtis_pcoa_results.qza",
jaccard_distance = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-jaccard_distance_matrix.qza",
jaccard_pcoa = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-jaccard_pcoa_results.qza",
unweighted_unifrac_mat = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza",
unweighted_unifrac_pcoa = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-unweighted_unifrac_pcoa_results.qza",
weighted_unifrac_mat = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-weighted_unifrac_distance_matrix.qza",
weighted_unifrac_pcoa = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-weighted_unifrac_pcoa_results.qza",
weighted_unifrac_viz = OUTPUTDIR + "/core-metrics-results/" + PROJ + ANALYSIS + "-weighted-unifrac-group-site-significance.qzv",
unweighted_unifrac_viz = OUTPUTDIR + "/core-metrics-results/" + PROJ + ANALYSIS + "-unweighted-unifrac-group-site-significance.qzv",
bray_curtis_emperor = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-bray_curtis_emperor.qzv",
jaccard_emperor = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-jaccard_emperor.qzv",
unweighted_unifrac = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-unweighted_unifrac_emperor.qzv",
weighted_unifrac = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-weighted_unifrac_emperor.qzv",   
shannon_signif = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-shannon-significance.qzv",
shannon_correl = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-shannon-significance-association.qzv",
observed_asv_signif = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-observed_otus-significance.qzv",
observed_asv_correl = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-observed-otus-significance-association.qzv",
bray_curtis_signif = OUTPUTDIR + "/core-metrics-results/" + PROJ + ANALYSIS + "-bray-curtis-group-significance.qzv",
barplots = OUTPUTDIR + "/core-metrics-results/" + PROJ + "-taxa-bar-plots.qzv",

#picrust2
picrust2tree = OUTPUTDIR + "/picrust2/out.tre",
marker = OUTPUTDIR + "/picrust2/marker_predicted_and_nsti.tsv.gz",
EC = OUTPUTDIR + "/picrust2/EC_predicted.tsv.gz",
metagenome_contrib = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_contrib.tsv.gz",
metagenome_unstrat = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
seqtab_norm = OUTPUTDIR + "/picrust2/EC_metagenome_out/seqtab_norm.tsv.gz",
weighted_nsti =OUTPUTDIR + "/picrust2/EC_metagenome_out/weighted_nsti.tsv.gz",
marker_describe = OUTPUTDIR + "/picrust2/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz",
path_abun_unstrat = OUTPUTDIR + "/picrust2/pathways_out/path_abun_unstrat.tsv.gz",
path_abun_unstrat_describe = OUTPUTDIR + "/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz",

# differential analysis
differentials = OUTPUTDIR + "/differential-expression/"  + PROJ + ANALYSIS + "-differential.qza",
viz_differentials = OUTPUTDIR + "/differential-expression/"  + PROJ + ANALYSIS + "-differential.qzv",
sig_diff = OUTPUTDIR + "/differential-expression/"  + PROJ + ANALYSIS + "-sig-differential.qza",
sig_table =  OUTPUTDIR + "/differential-expression/differentials.tsv",

# longitudinal
relative_frequency = OUTPUTDIR  + "/" + PROJ + "-relative_frequency_table.qza",
pw_diff = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "pairwise-differences.qzv",
pw_dist = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_pairwise-distances.qzv",
lme = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_linear-mixed-effects.qzv",
vol = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_volatility.qzv",
shannon_fd = OUTPUTDIR + "/longitudinal/" + PROJ + "_shannon-first-differences.qza",
first_diff = OUTPUTDIR + "/longitudinal/" + PROJ + "_first_distances.qza",
first_dist = OUTPUTDIR + "/longitudinal/" + PROJ + "_first_distances_LME.qzv",
baseline_fd = OUTPUTDIR + "/longitudinal/" + PROJ + "_first_distances_baseline_0.qza",
nmit = OUTPUTDIR + "/longitudinal/" + PROJ  + "_nmit-dm.qza",
nmit_viz = OUTPUTDIR + "/longitudinal/" + PROJ + ANALYSIS + "_nmit.qzv",
nmit_pcoa = OUTPUTDIR + "/longitudinal/" + PROJ  + "_nmit-pc.qza",
nmit_emp = OUTPUTDIR + "/longitudinal/" + PROJ  + "_nmit_emperor.qzv",
vol_table = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "-filtered_table.qza",
vol_est = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "-sample_estimator.qza",
vol_imp = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "-feature_importance.qza",
vol_acc = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "-accuracy_results.qzv",
vol_plot = OUTPUTDIR + "/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "-volatility_plot.qzv",
maturity_maz =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-maz_scores.qza",
maturity_emp =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-sample_estimator.qza",
maturity_imp =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-feature_importance.qza",
maturity_pred = OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-predictions.qza",
maturity_acc =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-accuracy_results.qzv",
maturity_vol =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-volatility_plots.qzv",
maturity_clus = OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-clustermap.qzv",
maturity_mod =  OUTPUTDIR + "/longitudinal/maturity/" + PROJ + ANALYSIS +  "-model_summary.qzv"

# List of non-longitudinal input files
rule_all_input_list = [# Preprocessing
CLASSIFIER,
raw_html,
raw_zip,
raw_multi_html,
trimmedData,
trim_html,
raw_qc,
trim_qc,
trim_multi_html,
q2_import,
q2_primerRM,
raw,
primer,
table,
cleaned_table,
filtered_rep,
filtered_table,
cleaned_metadata,
rep,
stats,
stats_viz,
sklearn,
filtered_sklearn,
biom,
tabulated_cleaned_metadata,
table_tsv,
table_tax,
seqs,
masked_align,relative_frequency,
rooted_tree,
align,
tree,
# Phylogeny
evenness_vector,
faith_pd_vector,
observed_asv,
rarefied_table,
id_filtered_table,shannon_vector,
bray_curtis,
bray_curtis_pcoa,
jaccard_distance,
jaccard_pcoa,
unweighted_unifrac_mat,
unweighted_unifrac_pcoa,
weighted_unifrac_mat,
weighted_unifrac_pcoa,
weighted_unifrac_viz,
unweighted_unifrac_viz,
bray_curtis_emperor,
jaccard_emperor,
unweighted_unifrac,
weighted_unifrac,
shannon_signif,
shannon_correl,
observed_asv_signif,
observed_asv_correl,
bray_curtis_signif,
barplots,
picrust2tree,
marker,
# Metagenomic (picrust2)
EC,
metagenome_contrib,
metagenome_unstrat,
seqtab_norm,
weighted_nsti,
marker_describe,
path_abun_unstrat,
path_abun_unstrat_describe,
differentials,
viz_differentials,
sig_diff,
sig_table]

# List of longitudinal input fils
rule_all_longitudinal_input = [pw_diff,
pw_dist,
lme,
vol,
shannon_fd,
first_diff,
first_dist,
baseline_fd,
nmit,
nmit_viz,
nmit_pcoa,
nmit_emp,
vol_table,
vol_est,
vol_imp,
vol_acc,
vol_plot,
maturity_maz,
maturity_emp,
maturity_imp,
maturity_pred,
maturity_acc,
maturity_vol,
maturity_clus,
maturity_mod]

if DROP_BLANKS == 'yes':

    print('dropping blanks')

else:
    print('no samples dropped')


if LONGITUDINAL == 'yes':

    rule_all_input_list.extend(rule_all_longitudinal_input)

    include: 'rules/longitudinal.smk'

    print("Will perform a longitudinal analysis")

else:
    print("no longitudinal analysis")


rule all:
    input:
        data = rule_all_input_list

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load rules #####

include: "rules/qc.smk"
include: "rules/dada2.smk"
include: "rules/phylogeny.smk"
include: "rules/picrust2.smk"
include: "rules/differential.smk"

