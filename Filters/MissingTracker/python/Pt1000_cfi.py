import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ( 
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/666A5021-C078-DD11-808F-001D09F23944.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/708E6A51-C078-DD11-8684-001D09F28D4A.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/7A4DD829-C078-DD11-A782-001D09F250AF.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/A649F93F-C078-DD11-A084-001D09F24047.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/BE2B016B-C078-DD11-9140-001D09F282F5.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/1E56386A-EB78-DD11-A81A-000423D6A6F4.root') );


secFiles.extend( (
               ) )
