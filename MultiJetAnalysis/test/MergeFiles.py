#! /usr/bin/env python
import os
import getopt
import sys

path        = "/pnfs/cms/WAX/11/store/user/kkousour/2011/Multijets"
prefix      = "dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kkousour/2011/Multijets/"
destination = "/uscms/home/kkousour/work/dataAnalysis/7TeV/2011/CMSSW_4_2_5/src/KKousour/MultiJetAnalysis/test/"

DIR = os.listdir(path)

ss = "flatTree"

for dd in DIR:
  print "Reading directory "+dd
  ROOTFILES = os.listdir(path+"/"+dd)
  command = "hadd -f "+destination+"/"+dd+"_"+ss+".root "
  for ll in ROOTFILES:
    command += prefix+dd+"/"+ll+" "
  print command
  os.system(command)
