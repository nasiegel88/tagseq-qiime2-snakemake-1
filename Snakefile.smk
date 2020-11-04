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
INPUTDIR = config["raw_data"]
SCRATCH = config["scratch"]
OUTPUTDIR = config["outputDIR"]

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
DB = config["database"]
DB_classifier = config["database_classified"]
DB_tax = config["database_tax"]

# Phylogeny
META = config["metadata"]
METACATEGORY = config["metadata-category"]
ALPHASTATISTIC = config["alpha-div-p-method"]
BETASTATISTIC = config["beta-div-p-method"]
PERMNUMBER = config["permutations"]

#----DEFINE RULES----#

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
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    stats = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qza",
    stats_viz = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-stats-dada2.qzv",
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ +	"-tax_sklearn.qza",
    biom = OUTPUTDIR + "/qiime2/asv/table/feature-table.biom",
    table_tsv = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.tsv",
    table_tax = OUTPUTDIR + "/qiime2/asv/tax_assigned/taxonomy.tsv",
    seqs = OUTPUTDIR + "/qiime2/asv/picrust2/dna-sequences.fasta",

    # new files
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
    weighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted-unifrac-group-site-significance.qzv",
    unweighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted-unifrac-group-site-significance.qzv",
    bray_curtis_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_emperor.qzv",
    jaccard_emperor = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-jaccard_emperor.qzv",
    unweighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_emperor.qzv",
    weighted_unifrac = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_emperor.qzv",   
    shannon_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance.qzv",
    shannon_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance-association.qzv",
    observed_asv_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus-significance.qzv",
    observed_asv_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed-otus-significance-association.qzv",
    bray_curtis_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray-curtis-group-significance.qzv",
    barplots = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-taxa-bar-plots.qzv",
    # stop!
  
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
    path_abun_unstrat_describe = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz"


####
#QC#
####

rule fastqc:
  input:    
    INPUTDIR + "/{sample}_{num}" + SUF
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule trimmomatic_pe:
  input:
    r1 = INPUTDIR + "/{sample}_" + R1_SUF + SUF,
    r2 = INPUTDIR + "/{sample}_" + R2_SUF + SUF
  output:
    r1 = SCRATCH + "/trimmed/{sample}_" + R1_SUF + "_trim.fastq.gz",
    r2 = SCRATCH + "/trimmed/{sample}_" + R2_SUF + "_trim.fastq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired = SCRATCH + "/trimmed/{sample}_1.unpaired.fastq.gz",
    r2_unpaired = SCRATCH + "/trimmed/{sample}_2.unpaired.fastq.gz"
  log:
    SCRATCH + "/trimmed/logs/trimmomatic/{sample}.log"
  params:
    trimmer = ["LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:25"],
    extra = ""
  wrapper:
    "0.35.2/bio/trimmomatic/pe"

rule fastqc_trim:
  input:
    SCRATCH + "/trimmed/{sample}_{num}_trim.fastq.gz"
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}_trimmed.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule multiqc:
  input:
    raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS)
  output:
    raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html", 
    raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", 
    trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"

  conda:
   "envs/multiqc-env.yaml"
  shell: 
    """
    multiqc -n multiqc.html {input.raw_qc} #run multiqc
    mv multiqc.html {output.raw_multi_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.raw_multi_stats} #move and rename stats
    rm -rf multiqc_data #clean-up
    #repeat for trimmed data
    multiqc -n multiqc.html {input.trim_qc} #run multiqc
    mv multiqc.html {output.trim_multi_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.trim_multi_stats} #move and rename stats
    rm -rf multiqc_data	#clean-up
    """ 

rule import_qiime:
  input:
    MANIFEST
  output:
    q2_import = SCRATCH + "/qiime2/" + PROJ + "-PE-demux.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
   "qiime tools import \
       --type 'SampleData[PairedEndSequencesWithQuality]' \
       --input-path {input} \
       --output-path {output.q2_import} \
       --input-format PairedEndFastqManifestPhred33"

rule rm_primers:
  input:
    q2_import = SCRATCH + "/qiime2/" + PROJ + "-PE-demux.qza"
  output:
    q2_primerRM = SCRATCH + "/qiime2/" + PROJ + "-PE-demux-noprimer.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_primer_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime cutadapt trim-paired \
       --i-demultiplexed-sequences {input.q2_import} \
       --p-front-f {config[primerF]} \
       --p-front-r {config[primerR]} \
       --p-error-rate {config[primer_err]} \
       --p-overlap {config[primer_overlap]} \
       --o-trimmed-sequences {output.q2_primerRM}"

rule get_stats:
  input:
    q2_import = SCRATCH + "/qiime2/" + PROJ + "-PE-demux.qza",
    q2_primerRM = SCRATCH + "/qiime2/" + PROJ + "-PE-demux-noprimer.qza"
  output:
    raw = OUTPUTDIR + "/qiime2/asv/viz/" + PROJ + "-PE-demux.qzv",
    primer = OUTPUTDIR + "/qiime2/asv/viz/" + PROJ + "-PE-demux-noprimer.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_getviz_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    """
     qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
     qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    """

#####
#ASV#
#####
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
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime dada2 denoise-paired \
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
        --o-denoising-stats {output.stats}"


rule drop_blanks:
  input:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza"
  output:
    cleaned_table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-no_blanks-asv-table.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-remove-blanks.qza"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
     """
     var=(config["remove_blanks"])
     if [ "${{var}}" == 'yes' ]; then

         qiime feature-table filter-samples
			      –i-table {input.table}
			      –m-metadata-file {config[metadata]}
			      --p-exclude-ids TRUE
			      --p-where config["blanks"]
			      –o-filtered-table {output.cleaned_table}

     elif [ "${{var}}" == 'no' ]; then

         qiime feature-table filter-samples
			      –i-table {input.table}
			      –m-metadata-file {config[metadata]}
			      --p-exclude-ids FALSE
			      –o-filtered-table {output.cleaned_table}

     else
         print('Error')
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
    "envs/qiime2-2019.10.yaml"
  shell:
   "qiime metadata tabulate \
       --m-input-file {input.stats} \
       --o-visualization {output.stats_viz}"

rule assign_tax:
  input:
    rep = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-rep-seqs.qza",
    db_classified = DB_classifier
  output:
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_sklearn_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime feature-classifier classify-sklearn \
	     --i-classifier {input.db_classified} \
	     --i-reads {input.rep} \
	     --o-classification {output.sklearn}"

rule gen_table:
  input:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza"
  output:
    table_tsv = OUTPUTDIR + "/qiime2/asv/table/feature-table.biom"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_exportBIOM_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/table")
  shell:
    "qiime tools export --input-path {input.table} --output-path {params}"

rule convert:
  input:
    OUTPUTDIR + "/qiime2/asv/table/feature-table.biom"
  output:
    OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.tsv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-exportTSV_q2.log"
  conda:
    "envs/qiime2-2019.10.yaml"
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
    "envs/qiime2-2019.10.yaml"
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
    "envs/qiime2-2019.10.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/picrust2")
  shell:
    "qiime tools export --input-path {input.rep} --output-path {params}"

########### 
#Phylogeny#
###########  

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
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input.rep} \
        --o-alignment {output.align} \
        --o-masked-alignment {output.masked_align} \
        --o-tree {output.tree} \
        --o-rooted-tree {output.rooted_tree}"

rule core_metrics:
  input:
    rooted_tree = OUTPUTDIR + "/qiime2/asv/mafft-fasttree-output/" + PROJ + "-rooted_tree.qza",
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza",
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
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity core-metrics-phylogenetic \
        --i-phylogeny {input.rooted_tree} \
        --i-table {input.table} \
        --p-sampling-depth {config[sampling-depth]} \
        --m-metadata-file {config[metadata]} \
        --o-rarefied-table {output.rarefied_table} \
        --o-faith-pd-vector {output.faith_pd_vector} \
        --o-observed-otus-vector {output.observed_otus} \
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
        --o-bray-curtis-emperor {output.bray_curtis_emperor}"

rule richness:
  input:
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon_vector.qza",
  output:
    shannon_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ +  "-shannon-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity alpha-group-significance \
        --i-alpha-diversity {input.shannon_vector} \
        --m-metadata-file {config[metadata]} \
        --o-visualization {output.shannon_signif}"

rule richcorr:
  input:
    shannon_vector = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon_vector.qza",
  output:
    shannon_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-shannon-significance-association.qzv",
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-shannon-significance-association.qzv"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity alpha-correlation \
        --i-alpha-diversity {input.shannon_vector} \
        --m-metadata-file {config[metadata]} \
        --p-method {config[alpha-div-p-method]} \
        --o-visualization {output.shannon_correl}"

rule asv_signif:
  input:
    observed_asv = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
  output:
    observed_asv_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-observed_otus-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity alpha-group-significance \
        --i-alpha-diversity {input.observed_asv} \
        --m-metadata-file {config[metadata]} \
        --o-visualization {output.observed_asv_signif}"

rule asv_corr:
  input:
    observed_asv = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed_otus_vector.qza",
  output:
    observed_asv_correl = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-observed-otus-significance-association.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-observed_otus-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity alpha-correlation \
        --i-alpha-diversity {input.observed_asv} \
        --m-metadata-file {config[metadata]} \
        --p-method {config[alpha-div-p-method]} \
        --o-visualization {output.observed_asv_correl}"

rule evenness:
  input:
    bray_curtis = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray_curtis_distance_matrix.qza",
  output:
    bray_curtis_signif = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-bray-curtis-group-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-bray-curtis-group-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity beta-group-significance \
        --i-distance-matrix {input.bray_curtis} \
        --m-metadata-file {config[metadata]} \
        --m-metadata-column {config[metadata-category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.bray_curtis_signif} \
        --p-no-pairwise"

rule unifrac:
  input:
    unweighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted_unifrac_distance_matrix.qza",
  output:
    unweighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-unweighted-unifrac-group-site-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-unweighted-unifrac-group-site-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity beta-group-significance \
        --i-distance-matrix {input.unweighted_unifrac_mat} \
        --m-metadata-file {config[metadata]} \
        --m-metadata-column {config[metadata-category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.unweighted_unifrac_viz} \
        --p-no-pairwise"

rule weighted_unifrac:
  input:
    weighted_unifrac_mat = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted_unifrac_distance_matrix.qza",
  output:
    weighted_unifrac_viz = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-weighted-unifrac-group-site-significance.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-weighted-unifrac-group-site-significance.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime diversity beta-group-significance \
        --i-distance-matrix {input.weighted_unifrac_mat} \
        --m-metadata-file {config[metadata]} \
        --m-metadata-column {config[metadata-category]} \
        --p-method {config[beta-div-p-method]} \
        --p-permutations {config[permutations]} \
        --o-visualization {output.weighted_unifrac_viz} \
        --p-no-pairwise"

rule barplot:
  input:
    table = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-asv-table.qza",
    sklearn = OUTPUTDIR + "/qiime2/asv/" + PROJ + "-tax_sklearn.qza",
  output:
    barplots = OUTPUTDIR + "/qiime2/asv/core-metrics-results/" + PROJ + "-taxa-bar-plots.qzv"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "-taxa-bar-plots.log"
  conda:
    "envs/qiime2-2019.10.yaml"
  shell:
    "qiime taxa barplot \
        --i-table {input.table} \
        --i-taxonomy {input.sklearn} \
        --m-metadata-file {config[metadata]} \
        --o-visualization {output.barplots}"

#########
#picrust#
#########

rule asv_reftree:
	input:
		seqs = OUTPUTDIR + "/qiime2/asv/picrust2/dna-sequences.fasta"
	output:
		picrust2tree = OUTPUTDIR + "/qiime2/asv/picrust2/out.tre"
	log:
		SCRATCH + "/qiime2/logs/" + PROJ + "_exportTREE_picrust.log"
	conda:
		"envs/picrust2-env.yaml"
	params:
		directory(OUTPUTDIR + "/qiime2/asv/picrust2/intermediate/place_seqs")
	shell:
		"place_seqs.py -s {input.seqs} -o {output.picrust2tree} -p 1 --intermediate {params}"

rule hsp:
	input:
		picrust2tree = OUTPUTDIR + "/qiime2/asv/picrust2/out.tre"
	output:
		marker = OUTPUTDIR + "/qiime2/asv/picrust2/marker_predicted_and_nsti.tsv.gz",
		EC = OUTPUTDIR + "/qiime2/asv/picrust2/EC_predicted.tsv.gz"
	log:
		SCRATCH + "/qiime2/logs/" + PROJ + "_marker_EC_predictions_picrust.log"
	conda:
		"envs/picrust2-env.yaml"
	shell:
	    """
			 hsp.py -i 16S -t {input.picrust2tree} -o {output.marker} -p 1 -n
			 hsp.py -i EC -t {input.picrust2tree} -o {output.EC} -p 1
	    """


rule metagenome:
	input:
		table_tsv = OUTPUTDIR + "/qiime2/asv/table/feature-table.biom",
		marker = OUTPUTDIR + "/qiime2/asv/picrust2/marker_predicted_and_nsti.tsv.gz",
		EC = OUTPUTDIR + "/qiime2/asv/picrust2/EC_predicted.tsv.gz"
	output:
		metagenome_contrib = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_contrib.tsv.gz",
		metagenome_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
		seqtab_norm = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/seqtab_norm.tsv.gz",
		weighted_nsti =OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/weighted_nsti.tsv.gz"
	log:
		SCRATCH + "/qiime2/logs/" + PROJ + "_marker_EC_predictions_picrust.log"
	conda:
		"envs/picrust2-env.yaml"
	params: 
		directory(OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out")
	shell:
		"metagenome_pipeline.py -i {input.table_tsv} -m {input.marker} -f {input.EC} \
			                       -o {params} --strat_out"

# add strat_out and metagenome_contrib to config.yaml			                       

rule pl_infer:
  input:
    metagenome_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
  output:
    path_abun_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat.tsv.gz"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_pre_metagenome_picrust.log"
  conda:
    "envs/picrust2-env.yaml"
  params:
    directory(OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out")
  shell:
    "pathway_pipeline.py -i {input.metagenome_unstrat} \
                          -o {params} -p 1"


rule add_describe:
	input:
		metagenome_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
		path_abun_unstrat = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat.tsv.gz"
	output:
		marker_describe = OUTPUTDIR + "/qiime2/asv/picrust2/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz",
		path_abun_unstrat_describe = OUTPUTDIR + "/qiime2/asv/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz"
	log:
		SCRATCH + "/qiime2/logs/" + PROJ + "_metagenome_description_picrust.log"
	conda:
		"envs/picrust2-env.yaml"
	shell:
	    """
	      add_descriptions.py -i {input.metagenome_unstrat} -m EC \
	                            -o {output.marker_describe}
	      add_descriptions.py -i {input.path_abun_unstrat} -m METACYC \
	                            -o {output.path_abun_unstrat_describe}
	    """