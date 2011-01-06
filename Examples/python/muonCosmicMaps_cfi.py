import FWCore.ParameterSet.Config as cms

cosmicCompatibility = cms.EDProducer(
    "MuonCosmicMaps",
    src=cms.InputTag("cosmicsVeto"),
    result = cms.string("cosmicCompatibility")
    )
timeCompatibility = cosmicCompatibility.clone(result = 'timeCompatibility')
backToBackCompatibility = cosmicCompatibility.clone(result = 'backToBackCompatibility')
overlapCompatibility = cosmicCompatibility.clone(result = 'overlapCompatibility')

cosmicCompatibilityLoader = cms.Sequence(
    cosmicCompatibility *
    timeCompatibility *
    backToBackCompatibility *
    overlapCompatibility
    )

def addUserData(patMuonProducer, labels=['cosmicCompatibility','timeCompatibility','backToBackCompatibility','overlapCompatibility']):
    for label in labels:
        patMuonProducer.userData.userFloats.src.append( cms.InputTag(label) )
