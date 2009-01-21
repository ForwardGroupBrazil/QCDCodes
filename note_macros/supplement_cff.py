import FWCore.ParameterSet.Config as cms

# Stand alone muon track producer
#import RecoMuon.StandAloneMuonProducer.standAloneMuons_cff
import RecoMuon.StandAloneMuonProducer.standAloneMuons_cfi
standAloneMuonsNoRPC = RecoMuon.StandAloneMuonProducer.standAloneMuons_cfi.standAloneMuons.clone()

standAloneMuonsNoRPC.ServiceParameters.RPCLayers = False
standAloneMuonsNoRPC.STATrajBuilderParameters.FilterParameters.EnableRPCMeasurement = False
standAloneMuonsNoRPC.STATrajBuilderParameters.BWFilterParameters.EnableRPCMeasurement = False

import UserCode.MuonTrackSelector.muontrackselector_cfi
goodTracks =  UserCode.MuonTrackSelector.muontrackselector_cfi.muontrackselector.clone()

supplement_rec_seq = cms.Sequence(standAloneMuonsNoRPC*goodTracks)

#------------

from SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi import *

tpToTkmutkTrackAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByHits'),
    label_tp = cms.InputTag('mergedtruth','MergedTrackTruth'),
    label_tr = cms.InputTag('goodTracks:TrackerOnly'),
)

tpToStaTrackNoUpAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mergedtruth', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuons')
#    label_tr = cms.InputTag('muonSta')
)

tpToStaTrackNoRPCAssociation = cms.EDProducer('TrackAssociatorEDProducer',
    associator = cms.string('TrackAssociatorByDeltaR'),
    label_tp = cms.InputTag('mergedtruth', 'MergedTrackTruth'),
    label_tr = cms.InputTag('standAloneMuonsNoRPC:UpdatedAtVtx')
#    label_tr = cms.InputTag('muonSta')
)

import SimMuon.MCTruth.MuonAssociatorByHits_cfi

tpToTkMuontkAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaMuonNoUpAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
tpToStaMuonNoRPCAssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()

tpToTkMuontkAssociation.tpTag = 'mergedtruth:MergedTrackTruth'
tpToTkMuontkAssociation.tracksTag = 'goodTracks:TrackerOnly'
tpToTkMuontkAssociation.UseTracker = True
tpToTkMuontkAssociation.UseMuon = False
tpToTkMuontkAssociation.EfficiencyCut_track = 0.5
tpToTkMuontkAssociation.PurityCut_track = 0.75

tpToStaMuonNoUpAssociation.tpTag = 'mergedtruth:MergedTrackTruth'
tpToStaMuonNoUpAssociation.tracksTag = 'standAloneMuons'
tpToStaMuonNoUpAssociation.UseTracker = False
tpToStaMuonNoUpAssociation.UseMuon = True
tpToStaMuonNoUpAssociation.EfficiencyCut_muon = 0.5
tpToStaMuonNoUpAssociation.PurityCut_muon = 0.5

tpToStaMuonNoRPCAssociation.tpTag = 'mergedtruth:MergedTrackTruth'
tpToStaMuonNoRPCAssociation.tracksTag = 'standAloneMuonsNoRPC:UpdatedAtVtx'
tpToStaMuonNoRPCAssociation.UseTracker = False
tpToStaMuonNoRPCAssociation.UseMuon = True
tpToStaMuonNoRPCAssociation.EfficiencyCut_muon = 0.5
tpToStaMuonNoRPCAssociation.PurityCut_muon = 0.5

supplement_assoc_seq = cms.Sequence(
    tpToTkmutkTrackAssociation
    +tpToStaTrackNoUpAssociation
    +tpToStaTrackNoRPCAssociation
    +tpToTkMuontkAssociation
    +tpToStaMuonNoUpAssociation
    +tpToStaMuonNoRPCAssociation
    )

#-----

import Validation.RecoMuon.MultiTrackValidator_cfi

trkMuonTrackTkVTrackAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()

trkMuonTrackTkVTrackAssoc.associatormap = 'tpToTkmutkTrackAssociation'
trkMuonTrackTkVTrackAssoc.associators = 'TrackAssociatorByHits'
trkMuonTrackTkVTrackAssoc.label = ('goodTracks:TrackerOnly',)

staMuonNoUpTrackVTrackAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()

staMuonNoUpTrackVTrackAssoc.associatormap = 'tpToStaTrackNoUpAssociation'
staMuonNoUpTrackVTrackAssoc.associators = 'TrackAssociatorByDeltaR'
staMuonNoUpTrackVTrackAssoc.label = ('standAloneMuons',)

staMuonNoRPCTrackVTrackAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()

staMuonNoRPCTrackVTrackAssoc.associatormap = 'tpToStaTrackNoRPCAssociation'
staMuonNoRPCTrackVTrackAssoc.associators = 'TrackAssociatorByDeltaR'
staMuonNoRPCTrackVTrackAssoc.label = ('standAloneMuonsNoRPC:UpdatedAtVtx',)

staMuonNoUpTrackVMuonAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()

staMuonNoUpTrackVMuonAssoc.associatormap = 'tpToStaMuonNoUpAssociation'
staMuonNoUpTrackVMuonAssoc.associators = 'MuonAssociationByHits'
staMuonNoUpTrackVMuonAssoc.label = ('standAloneMuons',)

staMuonNoRPCTrackVMuonAssoc = Validation.RecoMuon.MultiTrackValidator_cfi.multiTrackValidator.clone()

staMuonNoRPCTrackVMuonAssoc.associatormap = 'tpToStaMuonNoRPCAssociation'
staMuonNoRPCTrackVMuonAssoc.associators = 'MuonAssociationByHits'
staMuonNoRPCTrackVMuonAssoc.label = ('standAloneMuonsNoRPC:UpdatedAtVtx',)



#-----------------

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Validation.RecoMuon.RecoMuonValidator_cfi import *

import Validation.RecoMuon.RecoMuonValidator_cfi

recoMuonNoRPCVMuAssoc = Validation.RecoMuon.RecoMuonValidator_cfi.recoMuonValidator.clone()

recoMuonNoRPCVMuAssoc.subDir = 'RecoMuonV/RecoMuon_NoRPCMuonAssoc'

recoMuonNoRPCVMuAssoc.trkMuLabel = 'goodTracks:TrackerOnly'
recoMuonNoRPCVMuAssoc.staMuLabel = 'standAloneMuonsNoRPC:UpdatedAtVtx'
recoMuonNoRPCVMuAssoc.glbMuLabel = 'globalMuons'

recoMuonNoRPCVMuAssoc.trkMuAssocLabel = 'tpToTkMuontkAssociation'
recoMuonNoRPCVMuAssoc.staMuAssocLabel = 'tpToStaMuonNoRPCAssociation'
recoMuonNoRPCVMuAssoc.glbMuAssocLabel = 'tpToGlbMuonAssociation'

recoMuonNoUpVMuAssoc = Validation.RecoMuon.RecoMuonValidator_cfi.recoMuonValidator.clone()

recoMuonNoUpVMuAssoc.subDir = 'RecoMuonV/RecoMuon_NoUpMuonAssoc'

recoMuonNoUpVMuAssoc.trkMuLabel = 'goodTracks:TrackerOnly'
recoMuonNoUpVMuAssoc.staMuLabel = 'standAloneMuons'
recoMuonNoUpVMuAssoc.glbMuLabel = 'globalMuons'

recoMuonNoUpVMuAssoc.trkMuAssocLabel = 'tpToTkMuontkAssociation'
recoMuonNoUpVMuAssoc.staMuAssocLabel = 'tpToStaMuonNoUpAssociation'
recoMuonNoUpVMuAssoc.glbMuAssocLabel = 'tpToGlbMuonAssociation'

recoMuonNoUpVTrackAssoc = Validation.RecoMuon.RecoMuonValidator_cfi.recoMuonValidator.clone()

recoMuonNoUpVTrackAssoc.subDir = 'RecoMuonV/RecoMuon_NoUpTrackAssoc'

recoMuonNoUpVTrackAssoc.trkMuLabel = 'goodTracks:TrackerOnly'
recoMuonNoUpVTrackAssoc.staMuLabel = 'standAloneMuons'
recoMuonNoUpVTrackAssoc.glbMuLabel = 'globalMuons'

recoMuonNoUpVTrackAssoc.trkMuAssocLabel = 'tpToTkmutkTrackAssociation'
recoMuonNoUpVTrackAssoc.staMuAssocLabel = 'tpToStaTrackNoUpAssociation'
recoMuonNoUpVTrackAssoc.glbMuAssocLabel = 'tpToGlbTrackAssociation'

recoMuonNoRPCVTrackAssoc = Validation.RecoMuon.RecoMuonValidator_cfi.recoMuonValidator.clone()

recoMuonNoRPCVTrackAssoc.subDir = 'RecoMuonV/RecoMuon_NoRPCTrackAssoc'

recoMuonNoRPCVTrackAssoc.trkMuLabel = 'goodTracks:TrackerOnly'
recoMuonNoRPCVTrackAssoc.staMuLabel = 'standAloneMuons'
recoMuonNoRPCVTrackAssoc.glbMuLabel = 'globalMuons'

recoMuonNoRPCVTrackAssoc.trkMuAssocLabel = 'tpToTkmutkTrackAssociation'
recoMuonNoRPCVTrackAssoc.staMuAssocLabel = 'tpToStaTrackNoRPCAssociation'
recoMuonNoRPCVTrackAssoc.glbMuAssocLabel = 'tpToGlbTrackAssociation'

supplement_Valid_seq = cms.Sequence(
    trkMuonTrackTkVTrackAssoc
    +staMuonNoUpTrackVTrackAssoc
    +staMuonNoRPCTrackVMuonAssoc
    +staMuonNoUpTrackVMuonAssoc
    +staMuonNoRPCTrackVMuonAssoc
    +recoMuonNoRPCVMuAssoc
    +recoMuonNoUpVMuAssoc
    +recoMuonNoRPCVTrackAssoc
    +recoMuonNoUpVTrackAssoc
    )
