#! /usr/bin/env python
import os
from CMGTools.Production.eostools import *

path = '/store/cmst3/user/kkousour/CMG/'

datasets = ['VBF-Powheg115','VBF-Powheg120','VBF-Powheg125','VBF-Powheg130','VBF-Powheg135','GluGlu-Powheg125']

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
