#!/bin/bash
for file in `find ../data -name "*.bin.hist.root"`; do
  root -b -q -l findLevel.C'("'${file}'",1)'
done
