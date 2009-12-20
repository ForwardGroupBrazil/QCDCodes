import FWCore.ParameterSet.Config as cms

#redo the DigiLinks too
from SimTracker.Configuration.SimTracker_cff import *
#from IOMC.RandomEngine.IOMC_cff import *
RandomNumberGeneratorService = cms.Service(
    "RandomNumberGeneratorService",
    restoreStateLabel = cms.untracked.string('randomEngineStateProducer'),
    simSiPixelDigis = cms.PSet(
    initialSeed = cms.untracked.uint32(1234567),
    engineName = cms.untracked.string('HepJamesRandom')
    ),
    simSiStripDigis = cms.PSet(
    initialSeed = cms.untracked.uint32(1234567),
    engineName = cms.untracked.string('HepJamesRandom')
    )
    )

from Configuration.StandardSequences.MixingNoPileUp_cff import *
from SimGeneral.TrackingAnalysis.trackingParticles_cfi import *

TPLink = cms.Sequence(mix*mergedtruth*trDigi) 
