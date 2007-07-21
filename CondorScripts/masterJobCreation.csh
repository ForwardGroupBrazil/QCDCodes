#!/usr/local/bin/tcsh

setenv subtype BSUB
setenv subcommand bsub

set tag="160-junk"
set cfgTemplate="MuonRecoTemplate.cfg"
set que=8nm

set particle_list="13 -13"
set pt_list="10 50 100 500 1000"

foreach particle ($particle_list)
    foreach pt ($pt_list)
	set Particle=$particle
	set Pt=$pt
#	
#	COMMON_DIR=${HOME}/w0/${tag}/Mu${particle}Pt${pt}
#	mkdir $COMMON_DIR
#	./makeCfg.sh $tag $cfgTemplate Mu${particle}Pt${pt}.txt 1 1 1 $COMMON_DIR
#
	./makeProject.sh $tag 1 1 1 $cfgTemplate Mu${particle}Pt${pt}.plain $que
    end
end

exit
