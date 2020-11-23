#!/bin/bash
## Generatate tar.gz file with all files in current directory

# tar cf file.tar.gz *.gz

# Preliminary sequencing run
curl -L https://osf.io/7uzjb/download -o file.tar.gz
tar xf file.tar.gz
rm file.tar.gz

# wget -c https://osf.io/7uzjb/download -O - | tar -xz