import FWCore.ParameterSet.Config as cms

process = cms.Process("Merge")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1500

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = cms.Source(
    "PoolSource",
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
    $inputFileNames
    #"file:file1.root",
    #"file:file2.root",
    )
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('out_$outputFileName'),
    )

process.outpath = cms.EndPath(process.out)
