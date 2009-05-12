import FWCore.ParameterSet.Config as cms

def SwitchToComboSeeds(process):
    process.l3SeedCombination = cms.EDFilter(
        "L3MuonTrajectorySeedCombiner",
        seedCollections = cms.VInputTag(
          cms.InputTag("hltTrajSeedIOHit"),
          cms.InputTag("hltTrajSeedOIState"),
          cms.InputTag("hltTrajSeedOIHit")
          )
        )
    SwitchToBaselinePP(process)
    process.hltTrajSeedIOHit = process.hltL3TrajectorySeed.clone()
    SwitchToOIState(process)
    process.hltTrajSeedOIState = process.hltL3TrajectorySeed.clone()
    SwitchToOIHit(process)
    process.hltTrajSeedOIHit = process.hltL3TrajectorySeed.clone()
    
    process.HLTL3muonrecoSequence.replace(process.hltL3TrajectorySeed, process.hltTrajSeedIOHit*process.hltTrajSeedOIState*process.hltTrajSeedOIHit*process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed = process.l3SeedCombination

def SwitchToAllCombined(process):
    process.hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltL3TrajectorySeed.TSGForRoadSearchOI = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPropagation = cms.PSet()
    
    process.hltL3TrajectorySeed.TSGFromCombinedSeeds = cms.PSet(
        ComponentName = cms.string("CombinedTSG"),
        PSetNames = cms.vstring('baselinePP','oiState','oiHit'),
        baselinePP = makeBaselinePP(),
        oiState= makeOIState(),
        oiHit= makeOIHit()
        )
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedSeeds"
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' , 'SteppingHelixPropagatorAlong',  'SmartPropagatorAnyOpposite' ,'PropagatorWithMaterial')

def makeBaseline():
    return cms.PSet(
        ComponentName = cms.string('CombinedTSG'),
        PSetNames = cms.vstring('firstTSG','secondTSG'),

        firstTSG = cms.PSet(
           ComponentName = cms.string('TSGFromOrderedHits'),
           OrderedHitsFactoryPSet = cms.PSet(
             ComponentName = cms.string('StandardHitTripletGenerator'),
             SeedingLayers = cms.string('PixelLayerTriplets'),
             GeneratorPSet = cms.PSet(
              useBending = cms.bool(True),
              useFixedPreFiltering = cms.bool(False),
              phiPreFiltering = cms.double(0.3),
              extraHitRPhitolerance = cms.double(0.06),
              useMultScattering = cms.bool(True),
              ComponentName = cms.string('PixelTripletHLTGenerator'),
              extraHitRZtolerance = cms.double(0.06)
              )
             ),
           TTRHBuilder = cms.string('WithTrackAngle')
           ),

        secondTSG = cms.PSet(
            ComponentName = cms.string('TSGFromOrderedHits'),
            OrderedHitsFactoryPSet = cms.PSet(
             ComponentName = cms.string('StandardHitPairGenerator'),
             SeedingLayers = cms.string('PixelLayerPairs')
             ),
            TTRHBuilder = cms.string('WithTrackAngle')
            )
        )
    
def makeBaselinePP():
    pset=makeBaseline()
    pset.PSetNames.append('thirdTSG')
    pset.thirdTSG = cms.PSet(
           ComponentName = cms.string('DualByEtaTSG'),
           PSetNames = cms.vstring('endcapTSG','barrelTSG'),
           barrelTSG = cms.PSet(    ),
           endcapTSG = cms.PSet(
            ComponentName = cms.string('TSGFromOrderedHits'),
            OrderedHitsFactoryPSet = cms.PSet(
              ComponentName = cms.string('StandardHitPairGenerator'),
              SeedingLayers = cms.string('MixedLayerPairs')
              ),
            TTRHBuilder = cms.string('WithTrackAngle')
            ),
           etaSeparation = cms.double(2.0)
          )
    return pset

def SwitchToBaseline(process):
    print "baseline is by default in the menu."
    
def SwitchToBaselinePP(process):
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaselinePP()

def makeOIState():
    from UserCode.L3Switches.ptEtaRescale import ptRange,etaRange,diagTerm,offDiagTerm

    return  cms.PSet(
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
    
def SwitchToOIState(process):
    #switch off a few things
    process.hltL3TrajectorySeed.TrackerSeedCleaner = cms.PSet()
    process.hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPropagation = cms.PSet()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()

    #    process.hltL3TrajectorySeed
    
    #switch on the OIstate
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGForRoadSearchOI"
    process.hltL3TrajectorySeed.TSGForRoadSearchOI =makeOIState()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' , 'SteppingHelixPropagatorAlong' )
        


def makeOIHit():
    return cms.PSet(
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

def SwitchToOIHit(process):
    
    process.hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltL3TrajectorySeed.TSGForRoadSearchOI = cms.PSet()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()

    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromPropagation"

    process.hltL3TrajectorySeed.TSGFromPropagation = cms.PSet(
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
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring( 'SmartPropagatorAnyOpposite' ,'PropagatorWithMaterial' )       
