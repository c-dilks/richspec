#!/bin/bash

job="spectra.bat";
log="logfiles"
mkdir -p $log

echo "generate batchfile: $job"
echo "Executable = $(which root)" > $job
echo "Universe = vanilla" >> $job
echo "notification = never" >> $job
echo "getenv = True" >> $job
echo "" >> $job


for file in datadir/*.root; do
  suffix=$(echo $file | sed 's/^.*\///g')
  echo "prepare to analyze $file"
  echo "Arguments = -b -q -l analyseSpectra.C(\\\"${file}\\\",1)" >> $job
  echo "Log    = ${log}/spectra.${suffix}.log" >> $job
  echo "Output = ${log}/spectra.${suffix}.out" >> $job
  echo "Error  = ${log}/spectra.${suffix}.err" >> $job
  echo "Queue" >> $job
  echo "" >> $job
done

echo "submitting $job..."
condor_submit $job


while [ 1 ]; do
  n=$(condor_q $(whoami) | tail -n1 | awk '{print $1}')
  if [ $n -gt 0 ]; then
    echo "$n condor jobs remain..."
    sleep 10
  else break
  fi
done
echo "----------DONE SPECTRA ANALYSIS----------"
echo "concatenate pdfs..."
sleep 3

gs -q -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
  -sOutputFile=datadir/allPlots.pdf datadir/run*.plots.pdf
./renameFile.sh datadir/allPlots.pdf
echo "---------DONE PRODUCING PDFS-------------"
