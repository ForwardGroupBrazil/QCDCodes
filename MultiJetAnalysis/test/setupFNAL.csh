#! /bin/csh
echo "Setting gLite"
source /uscmst1/prod/grid/gLite_SL5_CRAB_27x.csh
echo "Setting scram"
source /uscmst1/prod/sw/cms/setup/cshrc prod
echo "Setting environment"
eval `scramv1 runtime -csh`
echo "Setting CRAB"
source /uscmst1/prod/grid/CRAB/crab.csh
#source /uscmst1/prod/grid/CRAB_2_7_3/crab.csh
