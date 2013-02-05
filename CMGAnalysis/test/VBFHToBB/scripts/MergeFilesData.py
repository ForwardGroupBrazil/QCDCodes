#! /usr/bin/env python
import os
from CMGTools.Production.eostools import *

path = '/store/cmst3/user/kkousour/CMG/'

datasets = ['MultiJet-Run2012A','BJetPlusX-Run2012B','BJetPlusX-Run2012C','BJetPlusX-Run2012D','Jet-Run2012A','JetMon-Run2012B','JetMon-Run2012C','JetMon-Run2012D']

for ds in datasets:
  print 'Adding ROOT files from dataset: /'+ds
  files = ls_EOS(path+'/'+ds)
  command = 'hadd -f ~/workspace/private/data/'+ds+'.root '

  for ff in files:
    if (ff.find('.root') > 0) :
      pfn = lfnToPFN(ff)
      print '\"'+pfn+'\",'
      command += pfn + ' '

  os.system(command)
