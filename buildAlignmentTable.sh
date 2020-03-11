#!/bin/bash
# builds table used for alignment studies

if [ $# -ne 1 ]; then echo "usage: $0 [dataDir]"; exit; fi
dataDir=$1

function getVal { 
  grep -w $1 $configfile | awk '{print $3}' | sed 's/;$//'
}

pushd $dataDir
> alignment.dat
for datfile in *.dat; do
  echo "read $datfile"
  configfile=$(echo $datfile|sed 's/table\.dat$/log/')
  runnum=$(getVal runID)
  x=$(getVal x)
  y=$(getVal y)
  while read line; do
    echo "$runnum $x $y $line" >> alignment.dat
  done < $datfile
done
popd
