#!/bin/bash

# analyse spectra
#for file in datadir/*.bin.hist.root; do
  #root -b -q -l analyseSpectra.C'("'${file}'",1)'
#done
#echo "----------DONE SPECTRA ANALYSIS----------"

# concatenate pdfs
gs -q -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
  -sOutputFile=datadir/allPlots.pdf datadir/run*.plots.pdf
./renameFile.sh datadir/allPlots.pdf
echo "---------DONE: PRODUCED PDFS"

