import FWCore.ParameterSet.Config as cms

## Define some utilities to declare bins easily
def _nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
def _evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

def makeInclusiveMuonPlots(rebinFactor=1):
    return cms.PSet(
        # ---- Kinematics ----
        ptBins = _evenBins( 0, 100, 2 * rebinFactor),
        pBins  = _evenBins( 0, 500, 2  * rebinFactor),
        etaBins = _evenBins( -2.6, 2.6, 0.2 * rebinFactor),
        phiBins = _evenBins(-3.2,  3.2, 0.2 * rebinFactor),
        chargeBins = cms.vdouble(-2,0,2),
        # ---- Vertex ----
        dxyFineBins = _evenBins(-0.2, 0.2, 0.005), #  50um
        dzFineBins  = _evenBins(-0.5, 0.5, 0.010), # 100um
        dxyCoarseBins = _evenBins( -4,  4, 0.1), # 1mm
        dzCoarseBins  = _evenBins(-10, 10, 0.1), # 1mm
        # ---- Tracks ----
        pixelHitsBins       = _nBins(8,0,8),
        pixelLayersBins     = _nBins(5,0,5),
        trackerHitsBins     = _nBins(33,0,33),
        trackerLostHitsBins = _nBins(10,0,10),
        muonHitsBins        = _nBins(50,0,50),
        muonStationHitsBins = _nBins(20,0,20),
        muonBadHitsBins     = _nBins(20,0,20),
        globalHitsBins      = _nBins(80,0,80),
        trackerChi2nBins = _evenBins(0, 10, 0.2 * rebinFactor),
        muonChi2nBins    = _evenBins(0, 10, 0.2 * rebinFactor),
        globalChi2nBins  = _evenBins(0, 10, 0.2 * rebinFactor),
        # ---- Isolation ----
        isolationBins = _evenBins(0,  5, .25  * rebinFactor),
        relIsoBins    = _evenBins(0, .5, .025 * rebinFactor),
        # ---- Muon ID ----
        muonStationsBins    = _nBins(5,0,5), 
        segmentMatchesBins = _nBins(12,0,12),
        segmentCompatBins  = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        caloCompatBins     = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        boolBins = _nBins(2,-0.5,1.5),
        zBins = _nBins(100,-500.,500.),
        rBins = _nBins(100,0.,500.),
        rzXBins = cms.uint32(1000),
        rzXRange = cms.vdouble(-500.,500.),
        rzYBins = cms.uint32(500),
        rzYRange = cms.vdouble(0.,500.),
        pdgBins = _nBins(1000,-500,500),
    )

inclusiveMuonPlots = cms.EDAnalyzer("InclusiveMuonPlotsGENSIM",
    makeInclusiveMuonPlots(),
    muons     = cms.InputTag('muons'),
    particleSrc  = cms.InputTag("genParticles"),

    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
    selectionReco = cms.string(""),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    weight = cms.untracked.double(1.0),
    #mother = cms.untracked.int32(32),
    #daughter = cms.untracked.int32(13)                            

    #normalization   = cms.InputTag("countCollisionEvents"), ## read normalization from output of cms.EDProducer("EventCountProducer")
)

