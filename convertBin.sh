#!/bin/bash
if [ $# -ne 1 ]; then echo "usage: $0 [dataDir]"; exit; fi
pushd $1
for file in `ls *.bin`; do 
  echo "convert $file"
  if [ ! -f "${file}.hist.root" ]; then
    /home/drewkenjo/decode/bin2hist $file
  fi
done
cd ..
dirName=`echo $1|sed 's/^.*\///g'`
echo $dirName
tar czvf ${dirName}.tar.gz ${dirName}/*.{root,log}
popd
