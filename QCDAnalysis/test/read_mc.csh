#!/bin/csh
setenv SCRAM_ARCH slc5_amd64_gcc434
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd /uscms/home/kkousour/work/dataAnalysis/7TeV/2011/CMSSW_4_1_2/src/KKousour/QCDAnalysis/test/
eval `scramv1 runtime -csh`
cmsRun ProcessedTreeProducer_mc_cfg.py
