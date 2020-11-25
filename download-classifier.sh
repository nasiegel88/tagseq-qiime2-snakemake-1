#!/bin/bash
# Download classifier, reprentative sequences, and taxonomy
# As of 2020-08-12 the Silva database is reguarly curated

# Silva 138 99% OTUs full-length sequences
curl -O -J -L https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza

# Silva 138 SSURef NR99 full-length sequences
curl -O -J -L https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza

# Silva 138 SSURef NR99 full-length taxonomy
curl -O -J -L https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza
# 2020-08-12