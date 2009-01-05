import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames = cms.untracked.vstring(
        '/store/user/aeverett/noteReco2112/SingleMuPt0_500-step2/step2-SingleMuPt0_500_cfi_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_137.root',
    )
)

process.demo = cms.EDAnalyzer('ParamSwitch',
                              muonsTag = cms.InputTag("muons"),
                              selectionTag = cms.InputTag("AllGlobalMuons"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo.root"),
    )


process.p = cms.Path(process.demo)
