import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonPlots")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'GR09_R_35X_V3::All'
process.GlobalTag.globaltag = 'START3X_V26A::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:pat_MuPat.root'
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


## process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlots",
##     muons     = cms.InputTag('muons'),
##     selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
##     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
##     normalization   = cms.InputTag("countCollisionEvents"), # optional!
##     # ---- Kinematics ----
##     ptBins = evenBins( 0, 15, 0.25),
##     pBins  = evenBins( 0, 30, 0.5),
##     etaBins = evenBins( -2.6, 2.6, 0.2),
##     phiBins = evenBins(-3.2,  3.2, 0.2),
##     chargeBins = cms.vdouble(-2,0,2),
##     # ---- Vertex ----
##     dxyFineBins = evenBins(-0.5,0.5, 0.005), #  50um
##     dzFineBins  = evenBins(-1.0,1.0, 0.010), # 100um
##     dxyCoarseBins = evenBins(-10,10,0.1), # 1mm
##     dzCoarseBins  = evenBins(-30,30,0.1), # 1mm
##     # ---- Tracks ----
##     pixelHitsBins   = nBins(8,0,8),
##     trackerHitsBins = nBins(33,0,33),
##     muonHitsBins    = nBins(25,0,50),
##     globalHitsBins  = nBins(40,0,80),
##     trackerChi2nBins = evenBins(0, 10, 0.2),
##     muonChi2nBins    = evenBins(0, 10, 0.2),
##     globalChi2nBins  = evenBins(0, 10, 0.2),
##     # ---- Isolation ----
##     isolationBins = evenBins(0, 5, .5),
## )
## process.globalMuons = process.trackerMuons.clone(
##     selection = "isGlobalMuon"
## )

from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import makeInclusiveMuonPlots;
commonInputs = cms.PSet(
    muons     = cms.InputTag('patMuons'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlots",
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


process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.bptxAnd   = hltLevel1GTSeed.clone(L1SeedsLogicalExpression = cms.string('0'))
process.bscFilter = hltLevel1GTSeed.clone(L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)'))
process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 15 && position.Rho <= 2"),
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
process.physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])
process.countCollisionEvents = cms.EDProducer("EventCountProducer")
process.preFilter = cms.Sequence(process.bptxAnd * process.physDecl * process.bscFilter * process.oneGoodVertexFilter * process.noScraping * process.countCollisionEvents)


process.p = cms.Path(process.preFilter * process.trackerMuons * process.globalMuons * process.standAloneMuons)

if True: ## use when starting from pat::Muons
    process.p.remove(process.preFilter)
    #process.source.fileNames = [ 'file:patMuons.root' ]
    process.trackerMuons.muons = 'patMuons'
    process.globalMuons.muons  = 'patMuons'
    process.TFileService.fileName = 'inclusiveMuonPlots_Data.pat.root'

