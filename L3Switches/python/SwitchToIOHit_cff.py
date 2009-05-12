import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.HLT_2E30_cff import hltL3TrajectorySeedOrig

#from regularL3_cff import hltL3Muons
from regularL3_cff import *

hltL3TrajectorySeed = cms.EDFilter("L3SeedCombiner",
    seedCollections = cms.VInputTag( 
        cms.InputTag("hltL3TrajectorySeedOrig"),
    )
)

comboSeeds_seq = cms.Sequence(
    hltL3TrajectorySeedOrig
    )
