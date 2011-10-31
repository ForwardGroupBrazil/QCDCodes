#! /usr/bin/env python
import os
import getopt
import sys

path        = "/pnfs/cms/WAX/11/store/user/kkousour/2011/Jets/mc/processed/"
prefix      = "dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr/cms/WAX/11/store/user/kkousour/2011/Jets/mc/processed/"
destination = "/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/"

DIR = os.listdir(path)

ss = "InclusiveJetsTree_mc"

for dd in DIR:
  print "Reading directory "+dd
  ROOTFILES = os.listdir(path+"/"+dd)
  command = "hadd -f "+destination+"/"+dd+"_"+ss+".root "
  for ll in ROOTFILES:
    command += prefix+dd+"/"+ll+" "
  print command
  os.system(command)
