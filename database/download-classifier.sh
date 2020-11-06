#!/bin/bash
# Download classifier, reprentative sequences, and taxonomy
# As of 2020-08-12 the Silva database is reguarly curated
# Silva 138 SSURef NR99 full-length sequences
wget -O "silva-138-99-seqs.qza" "https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza"
# Silva 138 SSURef NR99 full-length taxonomy
wget -O "silva-138-99-tax.qza" "https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza"



