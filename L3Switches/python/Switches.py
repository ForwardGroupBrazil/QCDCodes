import FWCore.ParameterSet.Config as cms

def PCut(process):
    process.hltL3TrajectorySeed.PCut = cms.double(2.5)

def cleanTSG(m):
    m.TSGForRoadSearchIOpxl = cms.PSet()
    m.TSGForRoadSearchOI = cms.PSet()
    m.TSGFromPropagation = cms.PSet()
    m.TSGFromMixedPairs = cms.PSet()
    m.TSGFromPixelPairs = cms.PSet()
    m.TSGFromPixelTriplets = cms.PSet()
    m.TSGFromCombinedHits = cms.PSet()


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


def regionBuilder():
    return cms.PSet(EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
                    Eta_fixed = cms.double( 0.2 ),
                    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
                    OnDemand = cms.double( -1.0 ),
                    Rescale_Dz = cms.double( 3.0 ),
                    Eta_min = cms.double( 0.1 ),
                    Rescale_phi = cms.double( 3.0 ),
                    PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
                    DeltaZ_Region = cms.double( 15.9 ),
                    Phi_min = cms.double( 0.1 ),
                    PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
                    vertexCollection = cms.InputTag( "pixelVertices" ),
                    Phi_fixed = cms.double( 0.2 ),
                    DeltaR = cms.double( 0.2 ),
                    EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
                    UseFixedRegion = cms.bool( False ),
                    Rescale_eta = cms.double( 3.0 ),
                    UseVertex = cms.bool( False ),
                    EscapePt = cms.double( 1.5 )
                    )

def seedCleaner():
    return cms.PSet(cleanerFromSharedHits = cms.bool( True ),
                    ptCleaner = cms.bool( True ),
                    TTRHBuilder = cms.string( "WithTrackAngle" ),
                    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
                    directionCleaner = cms.bool( True )
                    )
        
def SwitchToBaseline(process):
    PCut(process)
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaseline()
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedHits"
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()
    

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
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.TSGFromCombinedHits = makeBaselinePP()
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedHits"
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = regionBuilder()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = seedCleaner()


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
    cleanTSG(process.hltL3TrajectorySeed)

    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()
    process.hltL3TrajectorySeed.TrackerSeedCleaner = cms.PSet()

    #    process.hltL3TrajectorySeed
    
    #switch on the OIstate
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGForRoadSearchOI"
    process.hltL3TrajectorySeed.TSGForRoadSearchOI = makeOIState()
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
    cleanTSG(process.hltL3TrajectorySeed)
    process.hltL3TrajectorySeed.MuonTrackingRegionBuilder = cms.PSet()

    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromPropagation"

    process.hltL3TrajectorySeed.TSGFromPropagation = makeOIHit()
    process.hltL3TrajectorySeed.ServiceParameters.Propagators = cms.untracked.vstring()
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromPropagation)

################## all combined in one module


def SwitchToAllCombined(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
    
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

def SwitchToOICombined(process):
    cleanTSG(process.hltL3TrajectorySeed)
    PCut(process)
    
    process.hltL3TrajectorySeed.TSGFromCombinedSeeds = cms.PSet(
        ComponentName = cms.string("CombinedTSG"),
        PSetNames = cms.vstring('oiState','oiHit'),
        oiState= makeOIState(),
        oiHit= makeOIHit()
        )
    process.hltL3TrajectorySeed.tkSeedGenerator = "TSGFromCombinedSeeds"

    OIStatePropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiState)
    OIHitPropagators(process,process.hltL3TrajectorySeed.TSGFromCombinedSeeds.oiHit)

############### all combined from different modules
    
def SwitchToComboSeeds(process):
    cleanTSG(process.hltL3TrajectorySeed)
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
