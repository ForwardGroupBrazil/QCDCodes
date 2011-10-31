#! /usr/bin/env python
import os

FILES = [ 'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_10_2_BqW.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_1_1_sh0.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_2_1_LCM.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_4_1_t5G.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_5_1_8H9.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_6_1_M1v.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_7_1_1C1.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_8_1_MDv.root',
          'axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_9_1_Pa4.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_10_1_MWP.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_1_1_Obr.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_2_1_DcH.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_3_1_Mav.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_4_1_LgO.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_5_1_0G7.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_6_1_9Nt.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_8_1_MOz.root',
          'axiaxi900_hulubar_ussbar_hulbaru_ubarddbar_9_1_GUz.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_10_1_9Ad.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_1_1_wRs.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_2_1_3fO.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_3_1_1wv.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_4_1_Ewn.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_5_1_WOw.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_6_1_zCQ.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_7_1_vCm.root',
          'axiaxi900_scalarp_scalarm_uubar_uubar_8_1_Em3.root'
        ]

for ff in FILES:
  #----- FROM CERN --------------
  #command = 'lcg-cp -v -n 10 \"srm://srm-cms.cern.ch:8443/srm/managerv2?SFN=/castor/cern.ch/user/a/aferapon/Axigluons/'+ff+'\"'+' \"file:////uscmst1b_scratch/lpc1/3DayLifetime/kkousour/'+ff+'\"'
  #print command
  #os.system(command)
  #----- TO FNAL ----------------
  command = 'lcg-cp -v -n 10 \"file:////uscmst1b_scratch/lpc1/3DayLifetime/kkousour/'+ff+'\"'+' \"srm://cmssrm.fnal.gov:8443/srm/managerv1?SFN=/11/store/user/kkousour/2011/Jets/mc/Axigluon/'+ff+'\"'
  print command
  os.system(command)
  
