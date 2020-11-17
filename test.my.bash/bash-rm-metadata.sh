#!/bin/bash

declare -a arr=("no")

if [ "${arr[@]}" == yes ]; then
    fgrep -v NS.Blank5 sample-metadata.tsv > outputfile.tsv
    echo "All blanks dropped"

elif [ "${arr[@]}" == 'no' ]; then
    cp sample-metadata.tsv outputfile.tsv
    echo "no blanks dropped"
fi