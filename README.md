# richspec
```
# rich
convertBin.sh [datadir]

# local
pushd ../data
[optional:] mv ./* ../dataOld/
scp [rich]:~/data/rich2/*.gz ./
for f in *.gz; do tar xzvf $f; done; rm *.gz
popd
rm datadir
ln -s ../data/[datadir] ./datadir
loopAnalyseSpectra.sh
tarPDFs.sh

# alignment
root -b -q alignmentAnalysis.C
```
