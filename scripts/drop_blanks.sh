#!/bin/bash
var=(config["remove_blanks"])
if [ "${{var}}" == 'yes' ]; then

    qiime feature-table filter-samples \
        –i-table {input.table} \
        –m-metadata-file {config[metadata]} \
        --p-exclude-ids TRUE \
        --p-where config["blanks"] \ 
        –o-filtered-table {output.cleaned_table} 

elif [ "${{var}}" == 'no' ]; then

    qiime feature-table filter-samples \
        –i-table {input.table} \
        –m-metadata-file {config[metadata]} \
        --p-exclude-ids FALSE \ 
        –o-filtered-table {output.cleaned_table} 
else 

