#!/bin/bash
# Download classifier, reprentative sequences, and taxonomy
# As of 2020-08-12 the Silva database is reguarly curated

# Silva 138 SSURef NR99 full-length sequences
#curl -L https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza -o silva-138-99-seqs.qza

# Silva 138 SSURef NR99 full-length taxonomy
#curl -L https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza -o silva-138-99-tax.qza
# 2020-11-06

# Silva 138 SSURef NR99 full-length sequences
curl -L "https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza" -o "silva-138-99-seqs.qza"

# Silva 138 SSURef NR99 full-length taxonomy
curl -L "https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza" -o "silva-138-99-tax.qza"