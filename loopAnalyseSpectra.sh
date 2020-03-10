#!/bin/bash
# analyse spectra
for file in `find ../data -name "*.bin.hist.root"`; do
  root -b -q -l analyseSpectra.C'("'${file}'",1)'
done
echo "----------DONE SPECTRA ANALYSIS----------"

# concatenate pdfs
for dir in `ls -d ../data/*/`; do
  pushd $dir
  outpdf=allPlots_$(echo $dir | sed 's/\/$//' | sed 's/^.*\///g').pdf
  gs -q -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
  -sOutputFile=$outpdf run*.plots.pdf
  echo "PRODUCED ${dir}${outpdf}"
  popd
done
