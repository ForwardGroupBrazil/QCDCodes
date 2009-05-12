import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeedOrig

def SwitchToBaseline(process):
    process.hltTrajSeedIOHit   = hltL3TrajectorySeedOrig.clone()
    
def SwitchToBaselinePP(process):
    process.hltL3TrajSeedIOHit.TSGFromCombinedHits.PSetNames.append('thirdTSG')
    
def SwitchToOIState(process):
    process.hltTrajSeedOIState = hltL3TrajectorySeedOrig.clone()
    #switch off a few things
    process.hltTrajSeedOIState.TrackerSeedCleaner = cms.PSet()
    process.hltTrajSeedOIState.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltTrajSeedOIState.TSGFromPropagation = cms.PSet()
    process.hltTrajSeedOIState.MuonTrackingRegionBuilder = cms.PSet()
    process.hltTrajSeedOIState.TSGFromMixedPairs = cms.PSet()
    process.hltTrajSeedOIState.TSGFromPixelPairs = cms.PSet()
    process.hltTrajSeedOIState.TSGFromPixelTriplets = cms.PSet()
    process.hltTrajSeedOIState.TSGFromCombinedHits = cms.PSet()

    #    process.hltTrajSeedOIState
    
    #switch on the OIstate
    process.hltTrajSeedOIState.tkSeedGenerator = "TSGForRoadSearchOI"

    
    from ptEtaRescale import ptRange,etaRange,diagTerm,offDiagTerm

    from RecoMuon.TrackingTools.MuonErrorMatrixValues_cff import MuonErrorMatrixValues

    process.hltTrajSeedOIState.TSGForRoadSearchOI = cms.PSet(
        ComponentName = cms.string('TSGForRoadSearch'),
        option = cms.uint32(3),
        propagatorCompatibleName = cms.string('SteppingHelixPropagatorAny'),
        propagatorName = cms.string('SteppingHelixPropagatorAlong'),
        maxChi2 = cms.double(40.0),
        manySeeds = cms.bool(False),
        copyMuonRecHit = cms.bool(False),
        errorMatrixPset = cms.PSet(
                   action = cms.string('use'),
                   atIP = cms.bool(True),
                   errorMatrixValuesPSet = cms.PSet(
                       xAxis = ptRange,
                       yAxis = etaRange,
                       zAxis = cms.vdouble(-3.14159, 3.14159),
                       
                       pf3_V11 = diagTerm,
                       pf3_V22 = diagTerm,
                       pf3_V33 = diagTerm,
                       pf3_V44 = diagTerm,
                       pf3_V55 = diagTerm,

                       pf3_V12 = offDiagTerm,
                       pf3_V13 = offDiagTerm,
                       pf3_V14 = offDiagTerm,
                       pf3_V15 = offDiagTerm,
                       pf3_V23 = offDiagTerm,
                       pf3_V24 = offDiagTerm,
                       pf3_V25 = offDiagTerm,
                       pf3_V34 = offDiagTerm,
                       pf3_V35 = offDiagTerm,
                       pf3_V45 = offDiagTerm
                       )
                   )
        )
    process.hltTrajSeedOIState.ServiceParameters.Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' , 'SteppingHelixPropagatorAlong' )
        
    #process.MuonCkfTrajectoryBuilder.useSeedLayer = True
    #process.MuonCkfTrajectoryBuilder.rescaleErrorIfFail = 3.0



def SwitchToOIHit(process):
    process.hltTrajSeedOIHit   = hltL3TrajectorySeedOrig.clone()    
    process.hltTrajSeedOIHit.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltTrajSeedOIHit.TSGForRoadSearchOI = cms.PSet()
    process.hltTrajSeedOIHit.MuonTrackingRegionBuilder = cms.PSet()
    process.hltTrajSeedOIHit.TSGFromMixedPairs = cms.PSet()
    process.hltTrajSeedOIHit.TSGFromPixelPairs = cms.PSet()
    process.hltTrajSeedOIHit.TSGFromPixelTriplets = cms.PSet()
    process.hltTrajSeedOIHit.TSGFromCombinedHits = cms.PSet()

    process.hltTrajSeedOIHit.tkSeedGenerator = "TSGFromPropagation"

    process.hltTrajSeedOIHit.TSGFromPropagation = cms.PSet(
        ErrorRescaling = cms.double(3.0),
        ComponentName = cms.string('TSGFromPropagation'),
        errorMatrixPset = cms.PSet(),
        UpdateState = cms.bool(True),
        UseSecondMeasurements = cms.bool(False),
        SelectState = cms.bool(True),
        MaxChi2 = cms.double(15.0),
        UseVertexState = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyOpposite'),
        ResetRescaling = cms.bool(True),
        SigmaZ = cms.double(25.0),       
        )
    process.hltTrajSeedOIHit.ServiceParameters.Propagators = cms.untracked.vstring( 'SmartPropagatorAnyOpposite' )       

def SwitchToCombo(process):
    process.hltL3TrajectorySeed = cms.EDFilter(
        "L3SeedCombiner",
        seedCollections = cms.VInputTag( 
        cms.InputTag("hltTrajSeedIOHit"),
        cms.InputTag("hltTrajSeedOIState"),
        cms.InputTag("hltTrajSeedOIHit")
        )
        )

    process.comboSeeds_seq = cms.Sequence(
        process.hltTrajSeedIOHit
        + process.hltTrajSeedOIState
        + process.hltTrajSeedOIHit
        )
    
