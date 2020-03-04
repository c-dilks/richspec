# richspec
```
# rich
convertBin.sh [datadir]

# local
pushd ../data
mv ./* ../dataOld/
scp [rich]:~/data/rich2/*.gz ./
for f in *.gz; do tar xzvf $f; done; rm *.gz
popd
loopAnalyseSpectra.sh
tarPDFs.sh
```


