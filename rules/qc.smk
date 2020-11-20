# Quality control
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

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
   "../envs/multiqc-env.yaml"
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
    rm -rf multiqc_data #clean-up
    """ 

rule import_qiime:
  input:
    MANIFEST
  output:
    q2_import = SCRATCH + "/qiime2/" + PROJ + "-PE-demux.qza"
  log:
    SCRATCH + "/qiime2/logs/" + PROJ + "_q2.log"
  conda:
    "../envs/qiime2-2019.10.yaml"
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
    "../envs/qiime2-2019.10.yaml"
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
    "../envs/qiime2-2019.10.yaml"
  shell:
    """
     qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
     qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    """