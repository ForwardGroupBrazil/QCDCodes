# Auto generated configuration file
# using: 
# Revision: 1.99.2.8 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RAW2DIGI,RECO,POSTRECO,VALIDATION --mc -n100 --conditions FrontierConditions_GlobalTag,IDEAL_V11::All --filein file:step1.root --fileout step2_RECO_VALIDATION.root --eventcontent RECOSIM --datatier GEN-SIM-RECO --customise=Validation/RecoMuon/customise.py --python_filename=RelValSingleMuPt100_0.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SelStudy')

process.MessageLogger = cms.Service("MessageLogger",
    detailedInfo = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        extension = cms.untracked.string('.txt')
    ),
    debug = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        extension = cms.untracked.string('.txt'),
        SimLevelParent = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        MotherSearch = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        GlbSelectorStudy = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
        RecoMuonValidator = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)
        ),
    ),
    debugModules  = cms.untracked.vstring('glbSelStudy','recoMuonVTrackAssoc'),
    categories = cms.untracked.vstring('SimHitsAnlzrImproved',
        'GlbSelectorStudy',
        'MotherSearch',
        'SimLevelParent','RecoMuonValidator'),
    destinations = cms.untracked.vstring('detailedInfo','debug')
)


# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('step2 nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(25)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     '/store/relval/CMSSW_3_3_0/RelValInclusiveppMuX/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0002/889D6338-C4B7-DE11-BE10-001A928116FC.root',
#    '/store/user/aeverett/CMSSW_2_2_5//TTbar_Tauola//aeverett//TTbar_Tauola_CMSSW_2_2_5_IDEAL_step1//TTbar_Tauola_CMSSW_2_2_5_IDEAL_step1//95b7e266de590e364b9bb2c4d2a6f567//TTbar_Tauola_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37.root',
#    '/store/user/aeverett/SingleMuPt0_500_CMSSW_2_2_5_IDEAL_step1/SingleMuPt0_500_CMSSW_2_2_5_IDEAL_step1/3ea6209a8f9200d32c9696f6c7430454/SingleMuPt0_500_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_180.root',
#    '/store/user/aeverett/SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1/SingleKPt2_200_CMSSW_2_2_5_IDEAL_step1/abc92f2ca79c035bec8931c1df74b704/SingleKPt2_200_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19.root',
    )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('step2_RECO_VALIDATION.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)
process.output.outputCommands.append("keep *")
process.output.outputCommands.append("keep *_MEtoEDMConverter_*_*")

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V9::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.postreco_step = cms.Path(process.postreco_generator)
process.validation_step = cms.Path(process.validation)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

process.load("UserCode.GlbSelectorStudy.glbselectorstudy_cff")
process.glbSelStudy.doAssoc = False
process.glbSelStudy.trkMuAssocLabel = "tpToTkmuTrackAssociation"
process.glbSelStudy.staMuAssocLabel = "tpToStaTrackAssociation"
process.glbSelStudy.glbMuAssocLabel = "tpToGlbTrackAssociation"
process.glbSelStudy.trackProducer = "globalMuons"
process.glbSelStudy.trackAssociator = "TrackAssociatorByPosition"
process.p = cms.Path(process.muonAssociation_seq*process.glbSelStudy)
process.TFileService = cms.Service("TFileService", fileName = cms.string("TFS_selStudy.root"))

#process.muAnalyzer = cms.EDAnalyzer('MuAnalyzer',
#                                    simLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),)
#process.p2 = cms.Path(process.muAnalyzer)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.postreco_step,process.validation_step,process.p,process.endjob_step,process.out_step)




# Automatic addition of the customisation function

def customise(process):
    #process.load('sample')
    import Validation.RecoMuon.RelValCustoms as switch
    #switch.harvest_only(process)
    switch.validation_only(process)
    return(process)


# End of customisation function definition

process = customise(process)
