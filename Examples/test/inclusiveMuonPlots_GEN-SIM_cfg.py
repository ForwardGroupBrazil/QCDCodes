import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
from UserCode.Examples.inclusiveMuonPlotsGENSIM_cfi import makeInclusiveMuonPlots;
commonInputs = cms.PSet(
    muons     = cms.InputTag('patMuons'),
    particleSrc = cms.InputTag('mergedTruth'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.allGenSim = cms.EDAnalyzer(
    "InclusiveMuonPlotsGENSIM",
    makeInclusiveMuonPlots(),
    commonInputs,
    selection = cms.string(""),
)
process.allGenSim1 = process.allGenSim.clone(
    selection = "status == 1"
    )
process.genMuons = process.allGenSim.clone(
    selection = "status == 1 && abs(pdgId)==13"
    )
process.genPions = process.allGenSim.clone(
    selection = "status == 1 && abs(pdgId)==211"
    )
process.genPionsStrip = process.allGenSim.clone(
    selection = "status == 1 && abs(pdgId)==211"
    )
process.mergedMuons = process.allGenSim.clone(
    selection = "abs(pdgId)==13"
    )
process.mergedMuonsSi = process.allGenSim.clone(
    selection = "abs(pdgId) == 13 && (vertex.Rho < 120 && abs(vertex.Z) < 300)"
    )
process.mergedMuonsCal = process.allGenSim.clone(
    selection = "abs(pdgId) == 13 && (vertex.Rho > 120 || abs(vertex.Z) > 300)"
    )
process.mergedMuonsSi2 = process.allGenSim.clone(
    selection = "abs(pdgId) == 13 && collisionId == 2 && (vertex.Rho < 120 && abs(vertex.Z) < 300)"
    )
process.mergedMuonsCal2 = process.allGenSim.clone(
    selection = "abs(pdgId) == 13 && collisionId == 2 && (vertex.Rho > 120 || abs(vertex.Z) > 300)"
    )

process.motherFilter = cms.EDFilter(
    "CandViewRefSelector",
    filter = cms.bool(True),
    src = cms.InputTag("mergedTruth"),
    cut = cms.string('status == 1 && pt > 10 && pt < 20')
    )

process.p = cms.Path(process.allGenSim *
                     process.allGenSim1 *
                     process.genMuons *
                     process.genPions )

process.p2 = cms.Path(process.motherFilter +
                      (process.genPionsStrip *
                      process.mergedMuons *
                      process.mergedMuonsSi *
                      process.mergedMuonsCal *
                      process.mergedMuonsSi2 *
                      process.mergedMuonsCal2)
                      )

if False: ## use when starting from pat::Muons
    process.trackerMuons.muons = 'patMuons'
    process.globalMuons.muons  = 'patMuons'

process.load("MergedPiFiles")
process.TFileService.fileName = "inclusivePlots_Pi_01.root"
