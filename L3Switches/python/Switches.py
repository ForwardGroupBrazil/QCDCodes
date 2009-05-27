import FWCore.ParameterSet.Config as cms

def adaptTo2_2(process):
    if (not hasattr(process,"hltL3TkTracksFromL2")):
        process.hltL3TkTracksFromL2 = cms.EDProducer( "TrackProducer",
                                                      TrajectoryInEvent = cms.bool( True ),
                                                      useHitsSplitting = cms.bool( False ),
                                                      clusterRemovalInfo = cms.InputTag( "" ),
                                                      alias = cms.untracked.string( "" ),
                                                      Fitter = cms.string( "hltKFFittingSmoother" ),
                                                      Propagator = cms.string( "PropagatorWithMaterial" ),
                                                      src = cms.InputTag( "hltL3TrackCandidateFromL2" ),
                                                      beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
                                                      TTRHBuilder = cms.string( "WithTrackAngle" ),
                                                      AlgorithmName = cms.string( "undefAlgorithm" )
                                                      )
        process.hltL3Muons.L3TrajBuilderParameters.tkTrajLabel = cms.InputTag( "hltL3TkTracksFromL2" )
        process.HLTL3muonrecoNocandSequence.replace(process.hltL3Muons, process.hltL3TkTracksFromL2+process.hltL3Muons)
        ##end the necessary ES objects

        process.hltKFFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
                                              ComponentName = cms.string( "hltKFFitter" ),
                                              Propagator = cms.string( "PropagatorWithMaterial" ),
                                              Updator = cms.string( "KFUpdator" ),
                                              Estimator = cms.string( "Chi2" ),
                                              minHits = cms.int32( 3 ),
                                              appendToDataLabel = cms.string( "" )
                                              )
        process.hltKFFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
                                                       ComponentName = cms.string( "hltKFFittingSmoother" ),
                                                       Fitter = cms.string( "hltKFFitter" ),
                                                       Smoother = cms.string( "hltKFSmoother" ),
                                                       EstimateCut = cms.double( -1.0 ),
                                                       MinNumberOfHits = cms.int32( 5 ),
                                                       RejectTracks = cms.bool( True ),
                                                       BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
                                                       NoInvalidHitsBeginEnd = cms.bool( False ),
                                                       appendToDataLabel = cms.string( "" )
                                                       )
        process.hltKFSmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
                                                ComponentName = cms.string( "hltKFSmoother" ),
                                                Propagator = cms.string( "PropagatorWithMaterial" ),
                                                Updator = cms.string( "KFUpdator" ),
                                                Estimator = cms.string( "Chi2" ),
                                                errorRescaling = cms.double( 100.0 ),
                                                minHits = cms.int32( 3 ),
                                                appendToDataLabel = cms.string( "" )
                                                )
        
        ##additionnal parameters
        process.hltL3Muons.L3TrajBuilderParameters.ScaleTECxFactor = cms.double( -1.0 )
        process.hltL3Muons.L3TrajBuilderParameters.ScaleTECyFactor = cms.double( -1.0 )

    

def PCut(process):
    process.hltL3TrajectorySeed.PCut = cms.double(2.5)
    

############ baseline ##############

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


def SwitchToBaseline(process):
    PCut(process)
    print "baseline is by default in the menu."

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


############# baseline plus mixed seeding ######
    
def SwitchToBaselinePP(process):
    PCut(process)
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaselinePP()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()


############### OI-state based #############
    
def makeOIState():
    from UserCode.L3Switches.ptEtaRescale import ptRange,etaRange,diagTerm,offDiagTerm

    return  cms.PSet(
        ComponentName = cms.string('TSGForRoadSearch'),
        option = cms.uint32(3),
        propagatorCompatibleName = cms.string('SteppingHelixPropagatorOpposite'),
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

def OIStatePropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorCompatibleName.value())
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.propagatorName.value())
        
def SwitchToOIState(process):
    PCut(process)
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
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGForRoadSearchOI)



########## OI hit-based ###########

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
        ResetMethod = cms.string("discrete"),
        SigmaZ = cms.double(25.0),
        )

def OIHitPropagators(process,pset):
    if (not hasattr(process.hltL3TrajectorySeed.ServiceParameters,"Propagators")):
        process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append('PropagatorWithMaterial')
    process.hltL3TrajectorySeed.ServiceParameters.Propagators.append(pset.Propagator.value())

def SwitchToOIHit(process):
    PCut(process)
    process.hltL3TrajectorySeed.TSGForRoadSearchIOpxl = cms.PSet()
    process.hltL3TrajectorySeed.TSGForRoadSearchOI = cms.PSet()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromMixedPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelPairs = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromPixelTriplets = cms.PSet()
    process.hltL3TrajectorySeed.TSGFromCombinedHits = cms.PSet()

    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromPropagation"

    process.hltL3TrajectorySeed.TSGFromPropagation = makeOIHit()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromPropagation)

################## all combined in one module


def SwitchToAllCombined(process):
    PCut(process)
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

    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiState)
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiHit)


############### all combined from different modules
    
def SwitchToComboSeeds(process):
    PCut(process)
    process.l3SeedCombination = cms.EDProducer(
        "L3MuonTrajectorySeedCombiner",
        labels = cms.VInputTag(
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

    process.hltL3TrajectorySeed = process.l3SeedCombination
    process.HLTL3muonrecoSequence.replace(process.hltL3TrajectorySeed, process.hltTrajSeedIOHit + process.hltTrajSeedOIState + process.hltTrajSeedOIHit + process.hltL3TrajectorySeed)
