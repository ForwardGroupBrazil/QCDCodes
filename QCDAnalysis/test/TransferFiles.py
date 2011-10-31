#! /usr/bin/env python
import os

#setenv LCG_GFAL_INFOSYS lcg-bdii.cern.ch

#-------------------- DATA ------------------------------
location    = '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Oct11th/'
filename    = ['Jet_Run2011A_Aug05_ProcessedTree_data.root',
               'Jet_Run2011A_PromptV4_ProcessedTree_data.root',
               'Jet_Run2011A_May10_ProcessedTree_data.root',
               'Jet_Run2011A_PromptV6_ProcessedTree_data.root',
               'Jet_Run2011B_PromptV1_ProcessedTree_data.root'
              ] 
#location    = '/uscmst1b_scratch/lpc1/3DayLifetime/kkousour/'
#filename    = 'InclusiveJetsTree_data.root'

#print "Transfer DATA files to FNAL:"
#destination = '/11/store/user/kkousour/2011/Jets/production/Oct11th/'
#for ff in filename: 
#  command = "srmcp \"file:///"+location+ff+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+ff+"\""
#  print command
#  os.system(command)

print 'Transfer DATA files to CERN:'
destination = '/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Oct11th/'
for ff in filename:
  command = 'lcg-cp -v -n 10 \"file:///'+location+ff+'\" '+'\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN='+destination+ff+'\"'
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

