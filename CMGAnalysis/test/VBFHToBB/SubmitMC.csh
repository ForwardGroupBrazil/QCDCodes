#! /bin/csh

foreach SAMPLE ('VBF-Powheg115' 'VBF-Powheg120' 'VBF-Powheg125' 'VBF-Powheg130' 'VBF-Powheg135' 'GluGlu-Powheg125')
	echo $SAMPLE
        echo 'cleaning' 
        eos rm -r /eos/cms/store/cmst3/user/kkousour/CMG/$SAMPLE
        eos rm -r /eos/cms/store/cmst3/user/kkousour/CMG/$SAMPLE
        rm -r log$SAMPLE 
        echo 'submitting'
        cmsBatch.py 10 flat-$SAMPLE-cfg.py -o log$SAMPLE -r /store/cmst3/user/kkousour/CMG/$SAMPLE -b 'bsub -q 8nh < ./batchScript.sh'
end
