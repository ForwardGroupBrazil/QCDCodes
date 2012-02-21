#! /usr/bin/env python
import os

#setenv LCG_GFAL_INFOSYS lcg-bdii.cern.ch

#-------------------- DATA ------------------------------
location    = '/uscms_data/d2/kkousour/7TeV/2011/VBF/'
filename    = ['QCD_flatTree_weights.root'
               #'HT100_flatTree.root',
               #'HT250_flatTree.root',
               #'HT500_flatTree.root',
               #'HToBB_flatTree_weights.root'
              ] 
mc_location = '/uscmst1b_scratch/lpc1/3DayLifetime/kkousour/'
mc_filename = ['Pythia6Flat_ProcessedTree_mc.root',
               'HerwigFlat_ProcessedTree_mc.root'  
              ]

#print "Transfer DATA files to FNAL:"
#destination = '/11/store/user/kkousour/2011/VBF/'
#for ff in filename: 
#  command = "srmcp \"file:///"+location+ff+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+ff+"\""
#  print command
#  os.system(command)

print 'Transfer DATA files to CERN:'
destination = '/castor/cern.ch/user/k/kkousour/7TeV/vbf/'
for ff in filename:
  command = 'lcg-cp -v -n 10 \"file:///'+location+ff+'\" '+'\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN='+destination+ff+'\"'
  print command
  os.system(command)

#-------------------- MC-- ------------------------------
#print "Transfer MC file to FNAL:"
#destination = "/11/store/user/kkousour/2011/Jets/production/mc/"
#for ff in mc_filename: 
#  command = "srmcp \"file:///"+mc_location+ff+"\" "+"\"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN="+destination+ff+"\""
#  print command
#  os.system(command)

#print "Transfer MC file to CERN:"
#destination = "/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/mc/"
#for ff in mc_filename:
#  command = "lcg-cp -v -n 10 \"file:///"+mc_location+ff+"\" "+"\"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN="+destination+ff+"\""
#  print command
#  os.system(command)

