import FWCore.ParameterSet.Config as cms

## Define some utilities to declare bins easily
##--- n bins spaced evenly in range [min, max]
def _nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
##--- n bins of size 1 centered on 0, 1, ... (n-1)
def _iBins(n): return _nBins(n,-0.5,n-0.5)
##--- even bins of size delta from min to max
def _evenBins(min,max,delta): 
    ret = cms.vdouble(min)
    x = min
    while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
        x += delta
        ret.append(x)
    return ret 

def makeInclusiveMuonPlotsMRTU(rebinFactor=1):
    return cms.PSet(
        # ---- Tracks ----
        trackerChi2RelBins  = _evenBins(0, 4, 0.2 * rebinFactor),
        muonChi2RelBins  = _evenBins(0, 20, 0.2 * rebinFactor),
        chi2LocalPositionBins = _nBins(100,0.,0.01),
        chi2LocalMomentumBins = _nBins(100,0,100),
        localDistanceBins = _nBins(100,0,100),
        globalDeltaEtaPhiBins = _nBins(100,0,1.),
        glbTrackProbabilityBins = _nBins(100,0,100),
        
        # ---- Muon ID ----
        muonHitsBins        = _nBins(50,0,50),        
        muonStationHitsBins = _nBins(20,0,20),        
        muonStationsBins    = _nBins(5,0,5), 
        segmentMatchesBins = _nBins(12,0,12),
        segmentCompatBins  = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        caloCompatBins     = _evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
        # ---- Production Vertex ----
        zBins = _nBins(100,-500.,500.),
        rBins = _nBins(100,0.,500.),
        rzXBins = cms.uint32(1000),
        rzXRange = cms.vdouble(-500.,500.),
        rzYBins = cms.uint32(500),
        rzYRange = cms.vdouble(0.,500.),
        chi2ldXBins = cms.uint32(100),
        chi2ldXRange = cms.vdouble(0.,0.01),
        chi2ldYBins = cms.uint32(100),
        chi2ldYRange = cms.vdouble(0.,100.),
        chi2mldXBins = cms.uint32(100),
        chi2mldXRange = cms.vdouble(0.,500.),
        chi2mldYBins = cms.uint32(100),
        chi2mldYRange = cms.vdouble(0.,100.),
        # ---- ----
        boolBins = _nBins(2,-0.5,1.5),
        deltaPtBins = _evenBins( -50., 50., 2 * rebinFactor),
        deltaPtnBins = _evenBins(-5.,5.,0.1 * rebinFactor),
        muonHitCountsratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
        muonHitCountsrpcratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
        ratioBins = _evenBins(0.,1.2,0.1 * rebinFactor),
    )

inclusiveMuonPlotsMRTU = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlotsMRTU(),
    muons     = cms.InputTag('muons'),
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    #normalization   = cms.InputTag("countCollisionEvents"), ## read normalization from output of cms.EDProducer("EventCountProducer") 
)

