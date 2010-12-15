import FWCore.ParameterSet.Config as cms

muonCosmicMaps = cms.EDProducer("MuonCosmicMaps",
    src = cms.InputTag("muons"),
    useGlobalTrack = cms.bool(False), ## Use the standAlone track
    inputMuonCosmicCompatibilityValueMap = cms.InputTag("cosmicsVeto"),
)

## Helper function to add this info into a pat::Muon
def addUserData(patMuonProducer, label="muonCosmicMaps"):
    patMuonProducer.userData.userFloats.src += [
        cms.InputTag(label,""),
        cms.InputTag(label,"timeCompatibility"),
        cms.InputTag(label,"backToBackCompatibility"),
        cms.InputTag(label,"overlapCompatibility"),
        cms.InputTag(label,"ipCompatibility"),
        cms.InputTag(label,"vertexCompatibility"),
        ]
    
    
