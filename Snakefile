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
INPUTDIR = config["raw_data"]
SCRATCH = config["scratch"]
OUTPUTDIR = config["outputDIR"]
HOME = config["home"]
LONGITUDINAL = config['perform_longitudinal']

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

LONG_RULE_ALL  = list()

if LONGITUDINAL == "yes":
    LONG_RULE_ALL.append("pw_diff = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "pairwise-differences.qzv"",
    "pw_dist = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_pairwise-distances.qzv"",
    "lme = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_linear-mixed-effects.qzv"",
    "vol = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_volatility.qzv"",
    "shannon_fd = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_shannon-first-differences.qza"",
    "first_diff = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances.qza"",
    "first_dist = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances_LME.qzv"",
    "baseline_fd = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_first_distances_baseline_0.qza"",
    "nmit = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit-dm.qza"",
    "nmit_viz = OUTPUTDIR + "/qiime2/asv/longitudinal/" + PROJ + ANALYSIS + "_nmit.qzv"",
    "nmit_pcoa = OUTPUTDIR + "/qiime2/asv/longitudinal/output/" + PROJ + ANALYSIS + "_nmit-pc.qza"",
    "nmit_emp = OUTPUTDIR + "/qiime2/asv/longitudinal/output/" + PROJ + ANALYSIS + "_nmit_emperor.qzv"",
    "vol_table = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "filtered_table.qza"",
    "vol_est = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "sample_estimator.qza"",
    "vol_imp = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "feature_importance.qza"",
    "vol_acc = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "accuracy_results.qzv"",
    "vol_plot = OUTPUTDIR + "/qiime2/asv/longitudinal/mld_feat_volatility/" + PROJ + ANALYSIS + "volatility_plot.qzv"",
    "maturity_maz =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "maz_scores.qza"",
    "maturity_emp =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "sample_estimator.qza"",
    "maturity_imp =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "feature_importance.qza"",
    "maturity_pred = OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "predictions.qza"",
    "maturity_acc =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "accuracy_results.qzv"",
    "maturity_vol =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "volatility_plots.qzv"",
    "maturity_clus = OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "clustermap.qzv"",
    "maturity_mod =  OUTPUTDIR + "/qiime2/asv/longitudinal/maturity/" + PROJ + ANALYSIS +  "model_summary.qzv"")
else:
    print('null')


rule all:
  input:
    # fastqc output before trimming
    raw_html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    raw_zip = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html",
    raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    # Trimmed data output
    trimmedData = expand("{scratch}/trimmed/{sample}_{num}_trim.fastq.gz", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS), 
    trim_html = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name
    trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt",
    # QIIME2 outputs
    q2_import = SCRATCH + "/qiime2/" + PROJ + "-PE-demux.qza",
    q2_primerRM = SCRATCH + "/qiime2/" + PROJ + "-PE-demux-noprimer.qza",
    #vizualization stats
    raw = OUTPUTDIR + "/qiime2/asv/viz/" + PROJ + "-PE-demux.qzv",
    primer = OUTPUTDIR + "/qiime2/asv/viz/" + PROJ + "-PE-demux-noprimer.qzv",
    #ASV outputs:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza",
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza",
    cleaned_metadata = HOME + "noblank-sample-metadata.tsv",
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    stats = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qza",
    stats_viz = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qzv",
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ +	"-tax_sklearn.qza",
    biom = OUTPUTDIR + "/qiime2/asv/table/feature-table.biom",
    table_tsv = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.tsv",
    table_tax = OUTPUTDIR + "/qiime2/asv/tax_assigned/taxonomy.tsv",
    seqs = OUTPUTDIR + "/qiime2/asv/picrust2/dna-sequences.fasta",
    masked_align = OUTPUTDIR + "/qiime2/asv/" + "mafft-fasttree-output/" + PROJ + "-masked_alignment.qza",
    rooted_tree = OUTPUTDIR + "/qiime2/asv/" + "mafft-fasttree-output/" + PROJ + "-rooted_tree.qza",
    align = OUTPUTDIR + "/qiime2/asv/" + "mafft-fasttree-output/" + PROJ + "-alignment.qza",
    tree = OUTPUTDIR + "/qiime2/asv/" + "mafft-fasttree-output/" + PROJ + "-tree.qza",
    # core metrics
    evenness_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-evenness_vector.qza",
    faith_pd_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-faith_pd_vector.qza",
    observed_asv = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
    rarefied_table = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-rarefied_table.qza",
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon_vector.qza",
    bray_curtis = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_distance_matrix.qza",   
    bray_curtis_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_pcoa_results.qza",
    jaccard_distance = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_distance_matrix.qza",
    jaccard_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_pcoa_results.qza",
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza",
    unweighted_unifrac_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_pcoa_results.qza",
    weighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_distance_matrix.qza",
    weighted_unifrac_pcoa = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_pcoa_results.qza",
    weighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-weighted-unifrac-group-site-significance.qzv",
    unweighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-unweighted-unifrac-group-site-significance.qzv",
    bray_curtis_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_emperor.qzv",
    jaccard_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_emperor.qzv",
    unweighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_emperor.qzv",
    weighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_emperor.qzv",   
    shannon_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance.qzv",
    shannon_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance-association.qzv",
    observed_asv_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus-significance.qzv",
    observed_asv_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed-otus-significance-association.qzv",
    bray_curtis_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + ANALYSIS + "-bray-curtis-group-significance.qzv",
    barplots = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-taxa-bar-plots.qzv",
    #picrust2
    picrust2tree = OUTPUTDIR + "/qiime2/asv/picrust2/out.tre",
    marker = OUTPUTDIR + "/qiime2/asv/picrust2/marker_predicted_and_nsti.tsv.gz",
    EC = OUTPUTDIR + "/qiime2/asv/picrust2/EC_predicted.tsv.gz",
    metagenome_contrib = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_contrib.tsv.gz",
    metagenome_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    seqtab_norm = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/seqtab_norm.tsv.gz",
    weighted_nsti =OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/weighted_nsti.tsv.gz",
  	marker_describe = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz",
    path_abun_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat.tsv.gz",
    path_abun_unstrat_describe = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz",
    # longitudinal
    LONG_RULE_ALL


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load rules #####

include: "rules/qc.smk"
include: "rules/dada2.smk"
include: "rules/phylogeny.smk"
include: "rules/picrust2.smk"

if LONGITUDINAL == 'yes':
    include: 'rules/longitudinal.smk'
    print("Will perform a longitudinal analysis")
else:
    print("no longitudinal analysis")