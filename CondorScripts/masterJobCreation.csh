#!/usr/local/bin/tcsh

        setenv tag Mu131-test
        setenv particle_list  "13 -13"
        setenv pt_list  "10 50 100 500 1000"

        foreach particle ($particle_list)
            foreach pt ($pt_list)
                setenv Pt $pt
                setenv Particle $particle

		./makeCondorProject.sh $tag -1 -1 5 $CMSSW_BASE MuValidTemplate.cfg Mu${particle}Pt${pt}.plain
#               ./makeBatchProject.sh $tag -1 -1 5 $CMSSW_BASE MuValidTemplate.cfg Mu${particle}Pt${pt}.plain 8nh

            end
        end
exit
