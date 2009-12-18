import FWCore.ParameterSet.Config as cms

process = cms.Process("RecoMuon")
# Messages
#process.load("RecoMuon.Configuration.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")

process.load("RecoMuon.Configuration.RecoMuon_cff")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/uscms_data/d2/aeverett/RealData/CMSSW_3_3_5/src/myFastReco.root',
    )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myJustMuonReco.root')
                               )

#process.p = cms.Path(process.muonrecoComplete)

#process.load("ISpy/Analyzers/ISpy_Producer_cff")
#process.add_(
#    cms.Service("ISpyService",
#    outputFileName = cms.untracked.string('myJustMuonRecoSpy.ig'),
#    outputMaxEvents = cms.untracked.int32(100),
#    )
#)
#process.p1 = cms.Path(process.iSpy_sequence)

process.load("FastAnalysis.GlobalMatchingAnalyser.globalmatchinganalyser_cfi")

process.analyser_step = cms.Path(process.globalMatchingAnalyser) 

process.this_is_the_end = cms.EndPath(process.out)

process.GlobalTag.globaltag = 'GR09_P_V7::All'
