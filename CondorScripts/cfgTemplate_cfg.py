import FWCore.ParameterSet.Config as cms

process = cms.Process("validReco")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Services_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")

process.load("Configuration.EventContent.EventContent_cff")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32($skipEvents),
    fileNames = cms.untracked.vstring()
)
#--------------
#inputFileBlock

#--------------
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32($nEvents)
)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *',
        'keep *_*_*_validReco'),
    fileName = cms.untracked.string('$outFileName.root')
)

process.MyOut = cms.OutputModule("PoolOutputModule",
    process.FEVTDEBUGHLTEventContent,
    fileName = cms.untracked.string('$outFileName_edm.root')
)

#------------------
#recoBlock

#------------------

#------------------
#analyzerBlock

#------------------

process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)
process.outpath = cms.EndPath(process.MyOut)
process.schedule = cms.Schedule(process.pRecoOffline_step,process.postreco_step,process.postmuon_step,process.muonValidation_step,process.MEtoEDMConverter_step,process.outpath)

process.MessageLogger.categories = []
process.MessageLogger.FrameworkJobReport.FwkJob.limit = 0
process.MessageLogger.cerr.FwkReport.limit = 0

