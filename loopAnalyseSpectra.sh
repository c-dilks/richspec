#!/bin/bash

if [ $# -ne 1 ]; then echo "usage: $0 [dataDir]"; exit; fi
datadir=$1

# analyse spectra
for file in ${datadir}/*.bin.hist.root; do
  root -b -q -l analyseSpectra.C'("'${file}'",1)'
done
echo "----------DONE SPECTRA ANALYSIS----------"

# concatenate pdfs
outpdf=${datadir}/allPlots_$(echo $datadir | sed 's/\/$//' | sed 's/^.*\///g').pdf
gs -q -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
  -sOutputFile=$outpdf ${datadir}/run*.plots.pdf
echo "---------DONE: PRODUCED PDFS"
