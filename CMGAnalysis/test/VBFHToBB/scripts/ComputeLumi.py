#! /usr/bin/env python
import os, csv, sys

DIR = ['logJSONA','logJSONB','logJSONC','logJSOND']
NAME = ['json_MultiJet-Run2012A','json_BJetPlusX-Run2012B','json_BJetPlusX-Run2012C','json_BJetPlusX-Run2012D']

counter = 0
for d in DIR:
  print 'Readin csv files from directory: /'+d
  total_lumi = 0
  subdir = os.listdir(d)
  for ss in subdir:
    if (ss.find('Job_') == 0):
      files = os.listdir(d+'/'+ss)
      csv_name = NAME[counter]+'.csv'
      if (csv_name in files):
        #print csv_name  
        lumi_list_csv = open(d+'/'+ss+'/'+csv_name,'rb')
        for i in range(1):
          lumi_list_csv.next()
        lumi_dict = csv.DictReader(lumi_list_csv,delimiter=',',fieldnames=['Run','DeliveredLS','Delivered','SelectedLS','Recorded'])
        lumi = {}
        for l in lumi_dict:
          if (l['Recorded'] != 'n/a'):
            total_lumi+=float(l['Recorded'])      
  print 'Total Lumi = %1.3f /pb' % (total_lumi/1000000)
  counter += 1
