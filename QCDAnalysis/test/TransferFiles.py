#! /usr/bin/env python
import os

#-------------------- DATA ------------------------------
location    = "/uscms_data/d2/kkousour/7TeV/2011/Jets/data/June21st/"
filename    = "InclusiveJetsTree_data.root"

print "Transfer DATA file to FNAL:"
destination = "/11/store/user/kkousour/2011/Jets/production/June21st/"
command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+filename+"\""
print command

print "Transfer DATA file to CERN:"
destination = "/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/June21st/"
command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN="+destination+filename+"\""
print command

#-------------------- MC-- ------------------------------
location    = "/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/"
filename    = "InclusiveJetsTree_mc.root"

print "Transfer MC file to FNAL:"
destination = "/11/store/user/kkousour/2011/Jets/production/mc/"
command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+filename+"\""
print command

print "Transfer MC file to CERN:"
destination = "/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/mc/"
command = "srmcp \"file:///"+location+filename+"\" "+"\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN="+destination+filename+"\""
print command

