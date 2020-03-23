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
  * `[condor]loopAnalyseSpectra.sh` with `MODE=0` setting in `analyseSpectra.C`; check
    to see if the thresholds move around much (they shouldn't); also determine the
    maximum value of mu observed, and adjust `muMaxPlot` in `analyseSpectra.C`
    accordingly; if using the setting `SIMPLE>0`, you may need to adjust the value of
    `SIMPLE`.
    * for `726_728_729` data, use `SIMPLE=30`
    * for `753_760_752` data, use `SIMPLE=50`
  * `root -b -q analyseSpectra.C` with `MODE=1`, on your preferred run id file; this
    will store thresholds for each channel in `datadir/thresholds.dat`
  * `[condor]loopAnalyseSpectra.sh` with `MODE=2` setting in `analyseSpectra.C`; this
    will read thresholds which were stored in the previous step
