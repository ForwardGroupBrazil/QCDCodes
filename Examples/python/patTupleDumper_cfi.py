import FWCore.ParameterSet.Config as cms

def makePatTupleDumper():
    return cms.PSet(
        
        )

patTupleDumper = cms.EDAnalyzer(
    "PatTupleDumper",
    #makePatTupleDumper(),
    muons     = cms.InputTag('muons'),
    selection = cms.string(""),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    outputFileName = cms.string('ntuple.root')
)

