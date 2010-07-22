#!/usr/local/bin/tcsh

setenv subtype BSUB
setenv subcommand bsub

set tag="Validation"
set cfgTemplate="MuValidTemplate.cfg"
set cfgScriptName="makeCfg.sh"
set que=8nh

set firstFile=-1
set lastFile=-1
set filesPerJob=1

set particle_list="13 -13"
set pt_list="10 50 100 500 1000"

foreach particle ($particle_list)
    foreach pt ($pt_list)
	
	set COMMON_DIR=${HOME}/scratch0/${tag}/Mu${particle}Pt${pt}
	mkdir -p $COMMON_DIR
	setenv totalJobs `./${cfgScriptName} $tag $cfgTemplate Mu${particle}Pt${pt}.txt $firstFile $lastFile $filesPerJob $COMMON_DIR`
	echo "totalJobs = $totalJobs"

	./makeProject.sh $tag $firstFile $lastFile $filesPerJob $cfgTemplate Mu${particle}Pt${pt}.txt $que
    end
end

exit
