#!/bin/bash

DIR=QCDSimPlots_10_20
FILE="inclusivePlots_QCD_01.root"
REFFILE="inclusivePlots_QCD_01.root"

test -d $DIR || mkdir $DIR

for M in allGenSim allGenSim1 genPions genPionsStrip; do
    ./inclusiveMuonPlots_step2.py $FILE $M -t titles.txt -n none --select-group mykin --out $DIR/${M}
done

for M in mergedMuons; do
    ./inclusiveMuonPlots_step2.py $FILE $M -t titles.txt -n none --select-group mykin,prodvtx --out $DIR/${M}
done

for M in mergedMuons; do
    ./inclusiveMuonPlots_step2.py $FILE $M -r $FILE --rd genPionsStrip -t titles.txt -n none --select-group mykin --out $DIR/${M}_prob -R 
done
