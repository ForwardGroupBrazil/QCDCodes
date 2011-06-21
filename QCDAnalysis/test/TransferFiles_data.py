#! /usr/bin/env python
import os
import getopt
import sys

path = "/pnfs/cms/WAX/11/store/user/kkousour/2011/Jets/Processed/"
destination = "/uscms_data/d2/kkousour/7TeV/2011/Jets/data/June16th/"

DIR = os.listdir(path)

ss = "ProcessedTree_data"

for dd in DIR:
  print "Reading directory "+dd
  ROOTFILES = os.listdir(path+"/"+dd)
  for ll in ROOTFILES:
    command = "dccp "+path+dd+"/"+ll+" "+destination+"/"+dd+"_"+ll
    print command
    os.system(command)
  command = "hadd -f "+destination+"/"+dd+"_"+ss+".root "+destination+"/"+dd+"_"+ss+"_*.root"
  os.system(command)
  command = "rm "+destination+"/"+dd+"_"+ss+"_*.root"
  os.system(command)
