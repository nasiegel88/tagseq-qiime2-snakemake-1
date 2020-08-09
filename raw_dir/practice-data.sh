#!/bin/sh

## Generatate tar.gz file with all files in current direction
# tar cf file.tar.gz *.gz

# Preliminary sequencing run
curl -O -J -L https://osf.io/7uzjb/download
tar xf file.tar.gz
rm file.tar.gz