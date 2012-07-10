#! /usr/bin/env python
import os
import getopt
import sys

path        = "/store/cmst3/user/kkousour/ZJets_DY/"
prefix      = "root://eoscms//eos/cms/"
destination = "/afs/cern.ch/user/k/kkousour/scratch0/CMGTools/CMSSW_5_2_5/src/KKousour/ZJetAnalysis/test/"

ss = "flatTree"
NFILES = 33

ff = 0
command = "hadd -f "+destination+ss+".root "

while (ff < NFILES):
  command += prefix+path+"/"+ss+"_"+str(ff)+".root "
  ff+=1

print command
os.system(command)
