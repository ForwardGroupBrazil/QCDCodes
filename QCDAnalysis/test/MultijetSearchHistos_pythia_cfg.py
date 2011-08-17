import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('MultijetSearchHistos_pythia_Pt50.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.multijet    = cms.EDAnalyzer('MultijetSearchHistos',
    filename        = cms.string('/uscms_data/d2/kkousour/7TeV/2011/Jets/mc/Summer11PythiaZ2_InclusiveJetsTree_mc.root'),
    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak5'),
    triggers        = cms.vstring(),
    minPt           = cms.vdouble(50,50,50,50,50,50,50,50),
    maxEta          = cms.double(2.5),
    isMC            = cms.untracked.bool(True),
    rank            = cms.int32(8),
    jetID           = cms.int32(3),
    hcalNoiseFilter = cms.int32(0),
    nEvents         = cms.int32(-1),
    ptHatBnd        = cms.untracked.vdouble(5,15,30,50,80,120,170,300,470,600,800,1000,1400,1800,3500),
    ptHatLumi       = cms.untracked.vdouble(4.49e-5,1.32e-2,1.22e-1,1.04,8.40,5.32e+1,2.57e+2,5.51e+3,5.68e+4,2.73e+5,2.20e+6,6.30e+6,2.02e+8,8.20e+8)
)

process.p = cms.Path(process.multijet)

