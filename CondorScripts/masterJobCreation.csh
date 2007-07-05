#!/usr/local/bin/tcsh

        setenv reco_tag Mu131_HLT4-junk6
        setenv particle_list  "13 -13"
        setenv pt_list  "10 50 100 500 1000"

        foreach particle ($particle_list)
            foreach pt ($pt_list)
                setenv Pt $pt
                setenv Particle $particle

		./makeCondorProject.sh $reco_tag -1 -1 5 $CMSSW_BASE MuValidTemplate.cfg Mu${particle}Pt${pt}.plain

            end
        end
exit
