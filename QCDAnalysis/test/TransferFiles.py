#! /usr/bin/env python
import os

#setenv LCG_GFAL_INFOSYS lcg-bdii.cern.ch

#-------------------- DATA ------------------------------
location    = "/uscms_data/d2/kkousour/7TeV/2011/Jets/data/July8th/"
filename    = "InclusiveJetsTree_data.root"

print "Transfer DATA file to FNAL:"
destination = "/11/store/user/kkousour/2011/Jets/production/July8th/"
command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+filename+"\""
print command
os.system(command)

print "Transfer DATA file to CERN:"
destination = "/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/July8th/"
command = "lcg-cp -v -n 10 \"file:///"+location+filename+"\" "+"\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN="+destination+filename+"\""
print command
os.system(command)

#-------------------- MC-- ------------------------------
#location    = "/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/"
#filename    = "Summer11Herwig23Flat_InclusiveJetsTree_mc.root"

#print "Transfer MC file to FNAL:"
#destination = "/11/store/user/kkousour/2011/Jets/production/mc/"
#command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+filename+"\""
#print command
#os.system(command)

#print "Transfer MC file to CERN:"
#destination = "/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/mc/"
#command = "lcg-cp -v -n 10 \"file:///"+location+filename+"\" "+"\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN="+destination+filename+"\""
#print command
#os.system(command)

