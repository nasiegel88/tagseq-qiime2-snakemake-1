# differential abundance testing
## Noah Siegel
## UC Davis
## nasiegel@ucdavis.edu

rule aldex2:
    input:
        filtered_table = OUTPUTDIR + PROJ + "-filtered-asv-table.qza",
        cleaned_metadata = HOME + "noblank-sample-metadata.tsv"
    output:
        differentials = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-differential.qza"
    conda:
        "../envs/qiime2-2019.7.yml"
    shell:
        """
        qiime aldex2 aldex2 \
            --i-table {input.filtered_table} \
            --m-metadata-file {input.cleaned_metadata} \
            --m-metadata-column {config[metadata_category]} \
            --o-differentials {output.differentials}
        """

rule aldex2_effect:
    input:
        differentials = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-differential.qza"
    output:
        viz_differentials = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-differential.qzv"
    conda:
        "../envs/qiime2-2019.7.yml"
    shell:
        """
        qiime aldex2 effect-plot \
            --i-table {input.differentials} \
            --o-visualization {output.viz_differentials}
        """

rule extract_sig:
    input:
        differentials = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-differential.qza"
    output:
        sig_diff = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-sig-differential.qza"
    conda:
        "../envs/qiime2-2019.7.yml"
    shell:
        """
        qiime aldex2 extract-differences \
            --i-table {input.differentials} \
            --o-differentials {output.sig_diff} \
            --p-sig-threshold 0.8 \
            --p-effect-threshold 0 \
            --p-difference-threshold 0
        """

rule convert_sig:
    input:
        sig_diff = OUTPUTDIR + "differential-expression/"  + PROJ + ANALYSIS + "-sig-differential.qza"
    output: 
        sig_table =  OUTPUTDIR + "differential-expression/differentials.tsv"
    conda:
        "../envs/qiime2-2019.7.yml"
    params:
        directory(OUTPUTDIR + "differential-expression")
    shell:
        "qiime tools export --input-path {input.sig_diff} --output-path {params}"
