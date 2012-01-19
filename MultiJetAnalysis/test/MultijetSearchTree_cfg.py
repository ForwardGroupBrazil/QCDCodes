import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('/uscms_data/d2/kkousour/MultijetSearchTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.htUnprescaled = cms.EDAnalyzer('MultijetSearchTree',
    filenames         = cms.vstring('/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Nov13th/HT_Run2011A_May10_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Nov13th/HT_Run2011A_PromptV4_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Nov13th/HT_Run2011A_Aug05_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Nov13th/HT_Run2011A_PromptV6_ProcessedTree_data.root',
                                   '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Nov13th/HT_Run2011B_PromptV1_ProcessedTree_data.root'),
    treename          = cms.string('ProcessedTree'),
    dirname           = cms.string('ak5'),
    triggers          = cms.vstring('HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4',
                                    'HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8',
                                    'HLT_HT550_v9','HLT_HT550_v10','HLT_HT550_v11',
                                    'HLT_HT600_v1','HLT_HT600_v2','HLT_HT600_v3','HLT_HT600_v4',
                                    'HLT_HT650_v1','HLT_HT650_v2','HLT_HT650_v3','HLT_HT650_v4',
                                    'HLT_HT700_v1','HLT_HT700_v2',
                                    'HLT_HT750_v1','HLT_HT750_v2','HLT_HT750_v3'),
    jetID             = cms.int32(2),
    hcalNoiseFilter   = cms.int32(0),
    etaMAX            = cms.double(2.5),
    ptMIN             = cms.double(30),
    puTagMIN          = cms.double(0),
    isPreScaled       = cms.untracked.bool(False),
    nEvents           = cms.int32(-1)
)

process.ht450 = process.htUnprescaled.clone(
   triggers = cms.vstring('HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4','HLT_HT450_v5','HLT_HT450_v6','HLT_HT450_v7','HLT_HT450_v8','HLT_HT450_v9','HLT_HT450_v10','HLT_HT450_v11'), 
   isPreScaled = True
)

process.eightJet = process.htUnprescaled.clone(
   triggers = cms.vstring('HLT_EightJet35_v1','HLT_EightJet35_v2','HLT_EightJet35_v3','HLT_EightJet40_v1','HLT_EightJet40_v2','HLT_EightJet40_v3')
)

process.p = cms.Path(process.htUnprescaled * process.ht450 * process.eightJet)

