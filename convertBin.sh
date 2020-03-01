#!/bin/bash
if [ $# -ne 1 ]; then echo "usage: $0 [dataDir]"; exit; fi
pushd $1
for file in `ls *.bin`; do 
  echo $file
  if [ ! -f "${file}.hist.root" ]; then echo bin2hist $file; fi
done
cd ..
dirName=`echo $1|sed 's/^.*\///g'`
echo $dirName
#tar czvf ${dirName}.tar.gz ${dirName}/*.{root,log}
popd
