# richspec

## Convert bin files
```
convertBin.sh [datadir]
```

# Analyse data locally
```
pushd ../data
[optional:] mv ./* ../dataOld/
scp [rich]:~/data/rich2/*.gz ./
for f in *.gz; do tar xzvf $f; done; rm *.gz
popd
rm datadir
ln -s ../data/[datadir] ./datadir
[condor]loopAnalyseSpectra.sh // use condor version, if you have condor
```

# alignment
```
root -b -q alignmentAnalysis.C
```
* if you want to lock the threshold for all runs:
  * `root -b -q analyseSpectra.C` with `MODE=1`, on your preferred run id file; this
    will store thresholds for each channel
  * `[condor]loopAnalyseSpectra.sh` with `MODE=2` setting in `analyseSpectra.C`; this
    will read thresholds which were stored in the previous step
  * revert the setting to `MODE=0` afterward
