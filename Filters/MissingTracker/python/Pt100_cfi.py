import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ( 
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/7C85FE46-C078-DD11-9335-001D09F241B4.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/B261B730-C078-DD11-B188-001D09F2525D.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/C89D8F1A-C078-DD11-B936-001D09F24D4E.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0000/FA154F51-C078-DD11-B42F-001D09F2432B.root',
       '/store/relval/CMSSW_2_1_6/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/90CBF37C-EA78-DD11-B87D-001617DF785A.root') );


secFiles.extend( (
               ) )
