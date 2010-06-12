#!/bin/bash

DIR=PiRecoPlots_50_60
FILE="inclusiveMuonPlots_Pi_reco_02.root"
REFFILE="inclusivePlots_Pi_02.root"

COMP=" --composite='Ghost(8),Punch(1),Light(207),Heavy(67)' "

test -d $DIR || mkdir $DIR

for M in globalMuons globalselectMuons; do
    ./inclusiveMuonPlots_step2.py $FILE $M -t inclusiveMuonPlots_titles.txt -n none --select-group kin --out $DIR/${M}
done

for M in globalMuons globalselectMuons; do
    ./inclusiveMuonPlots_step2.py $FILE $M -r $FILE --rd $M $COMP -t inclusiveMuonPlots_titles.txt -n none --select-group kin --out $DIR/${M}_comp
done

for M in globalselectMuonsPunch globalselectMuonsLight; do
    ./inclusiveMuonPlots_step2.py $FILE $M -r $REFFILE --rd genPionsStrip -n none --select-group kin --out $DIR/${M}_ratio_genPi
done

for M in globalselectMuonsPunch globalselectMuonsLight; do
    ./inclusiveMuonPlots_step2.py $FILE $M -r $REFFILE --rd mergedMuons -n none --select-group kin --out $DIR/${M}_ratio_mergedMu
done