#!/bin/bash
# Download classifier, reprentative sequences, and taxonomy
# As of 2020-08-12 the Silva database is reguarly curated
# Silva 138 SSURef NR99 full-length sequences-2019.10
#wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza" --no-hsts 
# Silva 138 SSURef NR99 full-length sequences-2020.8
curl -sL \
 "https://data.qiime2.org/2020.8/common/silva-138-99-nb-classifier.qza" > \
  "silva-138-99-nb-classifier_2020.8.qza"