import FWCore.ParameterSet.Config as cms

allHits             = cms.EDProducer("ShowerID",
                                     src=cms.InputTag("muonShower"),
                                     result = cms.string("allHits")
                                     )

showerSize = allHits.clone(result = 'showerSize')
showerDeltaR = allHits.clone(result = 'showerDeltaR')
correlatedHits = allHits.clone(result = 'correlatedHits')

muonShowerInfoLoader = cms.Sequence(
    allHits+
    correlatedHits+
    showerSize+
    showerDeltaR
    )

def addUserData(patMuonProducer, label="showerInformation"):
   patMuonProducer.userData.userFloats.src += [
      'allHits','correlatedHits','showerSize','showerDeltaR']
   
