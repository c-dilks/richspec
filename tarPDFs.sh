#!/bin/bash
pushd ..
dir=plots.$(date +%y.%m.%d)
mv data $dir
tar czvf ${dir}{.tar.gz,/*/*.pdf}
mv $dir data
popd
