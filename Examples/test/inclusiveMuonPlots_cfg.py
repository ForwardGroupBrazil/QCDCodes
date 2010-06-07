import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START3X_V26A::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:pat_PiPat.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        fileMode    = cms.untracked.string("NOMERGE") )

## Define some utilities to declare bins easily
def nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
def evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

process.TFileService = cms.Service('TFileService',
    fileName=cms.string('inclusiveMuonPlots_Data.test.root')
)

###
# InclusiveMuonPlots can be run on PAT or RECO
# MuonAnalysis has general kinematic plots
# UserCode adds more figures about hits, segments, etc.
#from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlots;
commonInputs = cms.PSet(
    muons     = cms.InputTag('patMuons'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlots(),
    commonInputs,
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
)
process.globalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon"
    )
process.standAloneMuons = process.trackerMuons.clone(
    selection = "isStandAloneMuon"
    )
process.exclusiveGlobalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon && !isTrackerMuon"
    )
process.nonGlobalMuons = process.trackerMuons.clone(
    selection = "isTrackerMuon && isStandAloneMuon && !isGlobalMuon"
    )

process.p = cms.Path(process.trackerMuons *
                     process.globalMuons *
                     process.standAloneMuons *
                     process.exclusiveGlobalMuons *
                     process.nonGlobalMuons
                     )

if True: ## use when starting from pat::Muons
    process.trackerMuons.muons = 'patMuons'
    process.globalMuons.muons  = 'patMuons'


