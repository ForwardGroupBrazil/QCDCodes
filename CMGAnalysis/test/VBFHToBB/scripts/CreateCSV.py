#! /usr/bin/env python
import os, csv, sys

DIR = ['logJSONA','logJSONB','logJSONC','logJSOND']
NAME = ['json_MultiJet-Run2012A','json_BJetPlusX-Run2012B','json_BJetPlusX-Run2012C','json_BJetPlusX-Run2012D']

counter = 0
for d in DIR:
  print 'Creating csv files from directory: /'+d
  total_lumi = 0
  files = os.listdir(d)
  for ff in files:
    if (ff.find('Job_') == 0):
      csv_name = NAME[counter]+'.csv'
      command = 'lumiCalc2.py -i '+d+'/'+ff+'/'+NAME[counter]+'.txt -o '+d+'/'+ff+'/'+csv_name+' --nowarning overview'
      os.system(command)
      #lumi_list_csv = open(csv_name,'rb')
      #for i in range(1):
      #  lumi_list_csv.next()
      #lumi_dict = csv.DictReader(lumi_list_csv,delimiter=',',fieldnames=['Run','DeliveredLS','Delivered','SelectedLS','Recorded'])
      #lumi = {}
      #for l in lumi_dict:
      #  if (l['Recorded'] != 'n/a'):
      #    total_lumi+=float(l['Recorded'])
      #os.system('rm '+csv_name)      
  #print 'Total Lumi = %1.3f /pb' % (total_lumi/1000000)
  counter += 1
