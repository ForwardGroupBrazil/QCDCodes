import FWCore.ParameterSet.Config as cms

muonHitCounts = cms.EDProducer("MuonHitCounter",
    src = cms.InputTag("muons"),
    useGlobalTrack = cms.bool(False) ## Use the standAlone track
)

## Helper function to add this info into a pat::Muon
def addUserData(patMuonProducer, label="muonHitCounts"):
    patMuonProducer.userData.userInts.src += [
        cms.InputTag(label,""),
        cms.InputTag(label,"any"),
        cms.InputTag(label,"v1"),
        cms.InputTag(label,"v2"),
        cms.InputTag(label,"v3"),
        cms.InputTag(label,"v4"),
        cms.InputTag(label,"1any"),
        cms.InputTag(label,"2any"),
        cms.InputTag(label,"3any"),
        cms.InputTag(label,"4any"),
        cms.InputTag(label,"dt1"),
        cms.InputTag(label,"dt2"),
        cms.InputTag(label,"dt3"),
        cms.InputTag(label,"dt4"),
        cms.InputTag(label,"dt1any"),
        cms.InputTag(label,"dt2any"),
        cms.InputTag(label,"dt3any"),
        cms.InputTag(label,"dt4any"),
        cms.InputTag(label,"csc1"),
        cms.InputTag(label,"csc2"),
        cms.InputTag(label,"csc3"),
        cms.InputTag(label,"csc4"),
        cms.InputTag(label,"csc1any"),
        cms.InputTag(label,"csc2any"),
        cms.InputTag(label,"csc3any"),
        cms.InputTag(label,"csc4any"),
        cms.InputTag(label,"rpc1"),
        cms.InputTag(label,"rpc2"),
        cms.InputTag(label,"rpc3"),
        cms.InputTag(label,"rpc4"),
        cms.InputTag(label,"rpc1any"),
        cms.InputTag(label,"rpc2any"),
        cms.InputTag(label,"rpc3any"),
        cms.InputTag(label,"rpc4any"),
        ]
    
    
