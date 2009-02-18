import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
$inputFileNames
]);
secFiles.extend( [
]);
#import FWCore.ParameterSet.Config as cms

process = cms.Process("reVal")
process.load("FWCore.MessageService.MessageLogger_cfi")

### standard includes
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'STARTUP_V1::All'
process.GlobalTag.globaltag = 'IDEAL_V11::All'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = source

### validation-specific includes
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("Validation.RecoMuon.muonValidation_cff")


process.endjob_step = cms.Path(process.endOfProcess)

process.load("DQMServices.Components.EDMtoMEConverter_cff")

process.load("Validation.Configuration.postValidation_cff")

### configuration MultiTrackValidator ###
#process.multiTrackValidator.outputFile = 'mtv_'


process.cutsRecoTracks.algorithm = cms.string('')
process.cutsRecoTracks.quality = cms.string('')

process.multiTrackValidator.associators = ['TrackAssociatorByHits']

process.multiTrackValidator.label = ['generalTracks']
if (process.multiTrackValidator.label[0] == 'generalTracks'):
    process.multiTrackValidator.UseAssociators = cms.bool(True)
else:
    process.multiTrackValidator.UseAssociators = cms.bool(True)
######

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TFS_$outputFileName"),
                                   )

#process.load("UserCode.ParamSwitch.paramswitch_cfi")
import UserCode.ParamSwitch.paramswitch_cfi
process.paramswitch1 = UserCode.ParamSwitch.paramswitch_cfi.paramswitch.clone()
process.paramswitch1.ptFlag = False
process.paramswitch1.minPtCut = 82.
process.paramswitch1.maxPtCut = 97.
process.paramswitch1.chiFlag = False
process.paramswitch1.chi2Cut = 3.0
process.paramswitch1.hitFlag = False
process.paramswitch1.hitCut = 3
process.paramswitch1.muonsTag = "selMuons"

process.paramswitch2 = UserCode.ParamSwitch.paramswitch_cfi.paramswitch.clone()
process.paramswitch2.ptFlag = True
process.paramswitch2.minPtCut = 80.
process.paramswitch2.maxPtCut = 95.
process.paramswitch2.chiFlag = False
process.paramswitch2.chi2Cut = 3.0
process.paramswitch2.hitFlag = False
process.paramswitch2.hitCut = 3

process.load("UserCode.PtFilter.PtFilter_cfi")
process.ptFilter.MinPtCut = 90.0
process.ptFilter.MaxPtCut = 101.0
process.filter1 = cms.Path(process.ptFilter)

process.nStaFilter = cms.EDFilter(
    "TrackCountFilter",
    src = cms.InputTag('standAloneMuons:UpdatedAtVtx'),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(2)
    )

process.nTrkFilter = cms.EDFilter(
    "TrackCountFilter",
    src = cms.InputTag('generalTracks'),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(2)
    )

process.load("UserCode.MissingHitsFilter.MissingHitsFilter_cfi")
process.missingHitsFilter.hitCut = 2
process.filter2 = cms.Path(process.nStaFilter+process.nTrkFilter+process.missingHitsFilter)

process.load("UserCode.ChiFilter.chiFilter_cfi")
process.chiFilter.chi2Cut = 3.0
process.filter3 = cms.Path(process.chiFilter)

process.selMuons = cms.EDFilter(
    "MuonSelector",
    src = cms.InputTag("muons"),
#    cut = cms.string('isGlobalMuon = 1 & (innerTrack.hitPattern.numberOfValidTrackerHits - combinedMuon.hitPattern.numberOfValidTrackerHits > 4 | outerTrack.recHitsSize - combinedMuon.hitPattern.numberOfValidMuonHits > 4)')
    cut = cms.string('isGlobalMuon = 1 & (((innerTrack.hitPattern.numberOfValidTrackerHits - combinedMuon.hitPattern.numberOfValidTrackerHits) >= 2  ) || ((outerTrack.recHitsSize - combinedMuon.hitPattern.numberOfValidMuonHits) >= 2))')
    )

process.selMuonTracks = cms.EDProducer("TrackSelector",
    src = cms.InputTag("globalMuons"),
#    particleType = cms.string('mu+'),
    cut = cms.string('') #('pt > 95.0 & pt < 103.0 ')
)



#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring (
#    FILENAMES
#    )
#)
#process.extend("RelValTTbar_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.digi2track = cms.Sequence(process.siPixelDigis
                                  *process.SiStripRawToDigis
                                  *(process.trackerlocalreco
                                    +process.muonlocalreco)
                                  *(process.ckftracks
                                    *process.muonreco_plus_isolation)
                                  *process.cutsRecoTracks
                                  ##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                  *process.multiTrackValidator
                                  *process.recoMuonValidation
#                                  *process.endOfProcess
                                  )
#redo also tracking particles
process.digi2track_and_TP = cms.Sequence(process.mix
                                         *process.trackingParticles
                                         *process.siPixelDigis
                                         *process.SiStripRawToDigis
                                         *(process.trackerlocalreco
                                           +process.muonlocalreco)
                                         *(process.ckftracks
                                           *process.muonreco_plus_isolation)
                                         *process.cutsRecoTracks
                                         ##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                         *process.multiTrackValidator
                                         *process.recoMuonValidation
#                                         *process.endOfProcess
                                         )

process.re_tracking = cms.Sequence(process.siPixelRecHits
                                   *process.siStripMatchedRecHits
                                   *(process.ckftracks
                                     *process.muonreco_plus_isolation)
                                   *process.cutsRecoTracks
                                   ##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                   *process.multiTrackValidator
                                   *process.recoMuonValidation
#                                   *process.endOfProcess
                                   )

process.re_tracking_and_TP = cms.Sequence(process.mix
                                          *process.trackingParticles
                                          *process.siPixelRecHits
                                          *process.siStripMatchedRecHits
                                          *(process.ckftracks
                                            *process.muonreco_plus_isolation)
                                          *process.cutsRecoTracks
                                          ##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                          *process.multiTrackValidator
                                          *process.recoMuonValidation
#                                          *process.endOfProcess
                                          )

if (process.multiTrackValidator.label[0] == 'generalTracks'):
    process.only_validation = cms.Sequence(##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                           process.multiTrackValidator
                                           *process.recoMuonValidation
#                                           *process.endOfProcess
                                           )
else:
    process.only_validation = cms.Sequence(process.cutsRecoTracks
                                           ##process.cutsTPEffic*process.cutsTPFake* these modules are now useless
                                           *process.multiTrackValidator
                                           *process.recoMuonValidation
#                                           *process.endOfProcess
                                           )
    
if (process.multiTrackValidator.label[0] == 'generalTracks'):
    process.only_validation_and_TP = cms.Sequence(process.mix
                                                  *process.trackingParticles
                                                  *process.selMuonTracks
                                                  *process.selMuons
                                                  *(process.nStaFilter+process.nTrkFilter+process.missingHitsFilter)
                                                  *process.paramswitch1
                                                  *process.paramswitch2
                                                  *process.recoMuonValidation
                                                  *process.multiTrackValidator
#                                                  *process.endOfProcess
                                                  )
else:
    process.only_validation_and_TP = cms.Sequence(process.mix
                                                  *process.trackingParticles
                                                  *process.cutsRecoTracks
                                                  *process.selMuonTracks
                                                  *process.selMuons
                                                  *(process.nStaFilter+process.nTrkFilter+process.missingHitsFilter)
                                                  *process.paramswitch1
                                                  *process.paramswitch2
                                                  *process.recoMuonValidation
                                                  *process.multiTrackValidator
                                                  #                                                  *process.endOfProcess
                                                  )

### customized versoin of the OutputModule
### it save the mininal information which is necessary to perform tracking validation (tracks, tracking particles, 
### digiSimLink,etc..)

process.customEventContent = cms.PSet(
     outputCommands = cms.untracked.vstring('drop *')
 )

process.customEventContent.outputCommands.extend(process.RecoTrackerRECO.outputCommands)
process.customEventContent.outputCommands.extend(process.BeamSpotRECO.outputCommands)
process.customEventContent.outputCommands.extend(process.SimGeneralFEVTDEBUG.outputCommands)
process.customEventContent.outputCommands.extend(process.RecoLocalTrackerRECO.outputCommands)
process.customEventContent.outputCommands.append('keep *_simSiStripDigis_*_*')
process.customEventContent.outputCommands.append('keep *_simSiPixelDigis_*_*')
process.customEventContent.outputCommands.append('drop SiStripDigiedmDetSetVector_simSiStripDigis_*_*')
process.customEventContent.outputCommands.append('drop PixelDigiedmDetSetVector_simSiPixelDigis_*_*')



process.OUTPUT = cms.OutputModule("PoolOutputModule",
                                  process.customEventContent,
                                  fileName = cms.untracked.string('fullOutput_$outputFileName')
                                  )

process.VALOUTPUT = cms.OutputModule("PoolOutputModule",
                                     outputCommands = cms.untracked.vstring('keep *', "drop *_MEtoEDMConverter_*_*", "keep *_MEtoEDMConverter_*_reVal"),

                                     fileName = cms.untracked.string('$outputFileName'),
                                     SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('filter2')
    )

)

ValidationSequence="only_validation_and_TP"

if ValidationSequence=="harvesting":
    process.DQMStore.collateHistograms = False

    process.dqmSaver.convention = 'Offline'

    process.dqmSaver.saveByRun = cms.untracked.int32(-1)
    process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
    process.dqmSaver.forceRunNumber = cms.untracked.int32(1)


    process.dqmSaver.workflow = "/IDEAL_30X/RelValTTbar/Validation"
    process.DQMStore.verbose=3

    process.options = cms.untracked.PSet(
        fileMode = cms.untracked.string('FULLMERGE')
        )
    for filter in (getattr(process,f) for f in process.filters_()):
        if hasattr(filter,"outputFile"):
            filter.outputFile=""

process.harvesting= cms.Sequence(
    process.EDMtoMEConverter
    *process.postValidation
    *process.dqmSaver)




### final path and endPath
process.p = cms.Path(process.only_validation_and_TP)
if ValidationSequence!="harvesting":
    process.outpath = cms.EndPath(process.VALOUTPUT)

if ValidationSequence!="harvesting":
    process.schedule = cms.Schedule(
        process.p,
        process.filter2,
        process.endjob_step,
        process.outpath
        )
else:
    process.schedule = cms.Schedule(
        process.p
        )




process.tpToGlbTrackAssociation.label_tr = 'selMuonTracks'

process.glbMuonTrackVTrackAssoc.label = ('selMuonTracks',)
#process.glbMuonTrackVTrackAssoc.UseAssociators = True

process.recoMuonVTrackAssoc.glbMuLabel = 'selMuonTracks'
process.recoMuonVTrackAssoc.muonLabel = 'selMuons'

#process.tpToGlbMuonAssociation.tracksTag = ('selMuonTracks',)

#process.recoMuonVTrackAssoc.doAssoc = True
