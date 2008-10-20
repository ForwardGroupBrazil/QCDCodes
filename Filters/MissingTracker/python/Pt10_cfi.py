import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ( 
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/18797612-C078-DD11-9B06-00304879FBB2.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/624AB625-C078-DD11-AF9C-001D09F24DDF.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/7CBB3E6C-C078-DD11-829F-0030487A3232.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/DC3C4C7D-C078-DD11-81CA-001D09F23A3E.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/D0A06A6D-EB78-DD11-85A9-000423D98DB4.root'
       ) );


secFiles.extend( (
               ) )
