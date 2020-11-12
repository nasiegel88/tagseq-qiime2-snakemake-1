#!/bin/bash

declare -a arr=("yes")

if [ "${arr[@]}" == yes ]; then
  qiime feature-table filter-samples \
  --i-table prelim-mld-fecal-asv-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-exclude-ids TRUE  \
  --p-where "SampleID IN ('NS.Blank5')"  \
  --o-filtered-table clean-prelim-mld-fecal-asv-table.qza

elif [ "${arr[@]}" == 'no' ]; then
  qiime feature-table filter-samples \
  --i-table prelim-mld-fecal-asv-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-exclude-ids FALSE \
  --o-filtered-table clean-prelim-mld-fecal-asv-table.qza
fi
