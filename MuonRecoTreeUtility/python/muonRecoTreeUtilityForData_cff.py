import FWCore.ParameterSet.Config as cms

#migration to tfile service at some point
#TFileService = cms.Service("TFileService", fileName = cms.string('TFSReco.root') )

from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_reco_cfi import *

## #from Validation.RecoMuon.associators_cff import tpToL3MuonAssociation,tpToL2MuonAssociation

## #reDIGI_Path = cms.Paht( !tkSimDigiLinkAreThere + trdigi )

dummyModule = cms.EDFilter("DummyModule")
MRTU_Path = cms.Path(dummyModule)

#has to be in the endapth so that it can get the TriggerResults object
#TimerService = cms.Service("TimerService",useCPUtime = cms.untracked.bool(True))
#import HLTrigger.Timer.timer_cfi
#hltTimer = HLTrigger.Timer.timer_cfi.myTimer.clone()
MRTU_EndPath = cms.EndPath( recoMuonTreeMaker )

## #MRTUSchedule = cms.Schedule( reDIGI_Path + MRTU_Path )
## MRTUSchedule = cms.Schedule( MRTU_Path )
## #remember to extend with the EndPath

recoMuonTreeMaker.isRecoLevel = True

def insertMRTU(process):
    process.load('Workspace.MuonRecoTreeUtility.muonRecoTreeUtilityForData_cff')
    import FWCore.ParameterSet.SequenceTypes
    for p in process.schedule:
        if (p.__class__==FWCore.ParameterSet.SequenceTypes.EndPath):
            process.schedule.insert(process.schedule.index(p), process.MRTU_Path )
            break
    process.schedule.append( process.MRTU_EndPath )

    ##actually do the --no_output option
    if (hasattr(process,"out_step")):
        process.schedule.remove(process.out_step)
