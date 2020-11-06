#!/bin/sh
## Generatate tar.gz file with all files in current directory

# tar cf file.tar.gz *.gz

# Preliminary sequencing run
#curl -O -J -L https://osf.io/7uzjb/download
#tar xf file.tar.gz
#rm file.tar.gz

 wget -c https://osf.io/7uzjb/download -O - | tar -xz