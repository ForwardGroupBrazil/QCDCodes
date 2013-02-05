#! /bin/csh

foreach SAMPLE ('MultiJet-Run2012A' 'BJetPlusX-Run2012B' 'BJetPlusX-Run2012C' 'BJetPlusX-Run2012D' 'Jet-Run2012A' 'JetMon-Run2012B' 'JetMon-Run2012C' 'JetMon-Run2012D')
	echo $SAMPLE
        echo 'cleaning' 
        eos rm -r /eos/cms/store/cmst3/user/kkousour/CMG/$SAMPLE
        eos rm -r /eos/cms/store/cmst3/user/kkousour/CMG/$SAMPLE
        rm -r log$SAMPLE 
        echo 'submitting'
        cmsBatch.py 20 flat-$SAMPLE-cfg.py -o log$SAMPLE -r /store/cmst3/user/kkousour/CMG/$SAMPLE -b 'bsub -q 8nh < ./batchScript.sh'
end
