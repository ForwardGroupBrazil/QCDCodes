# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: MRTU -s L1,HLT:1E31 -n 100 --processName MHTU --conditions FrontierConditions_GlobalTag,MC_31X_V1::All --filein HLTDEBUG-file.root --no_output --customise Workspace/MuonRecoTreeUtility/muonRecoTreeUtility_customise.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('MRTU')

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
        extension = cms.untracked.string('.txt'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(00000)
        ),
        MuonRecoTreeUtility = cms.untracked.PSet(
            limit = cms.untracked.int32(000000)
        ),
        GlbSelectorStudy = cms.untracked.PSet(
     limit = cms.untracked.int32(0)
     ),
        SpecialBit = cms.untracked.PSet(
     limit = cms.untracked.int32(-1)
     )
    ),
    categories = cms.untracked.vstring(
     'MuonRecoTreeUtility','GlbSelectorStudy','SpecialBit'
     ),
    destinations = cms.untracked.vstring('detailedInfo','debug'),
    debugModules  = cms.untracked.vstring('hltMuonTreeMaker','glbSelStudy'),

)


# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/PostRecoGenerator_cff')
process.load('Configuration/StandardSequences/Validation_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.4 $'),
    annotation = cms.untracked.string('MHTU nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#     '/store/user/aeverett/CMSSW_3_2_1/SingleMuPt200//aeverett//SingleMuPt200_CMSSW_3_2_1_step1//SingleMuPt200_CMSSW_3_2_1_step1//3f9251f9631b9dac27d828bfadc5f286//step1_29.root',
#     "/store/user/aeverett//CMSSW_3_1_2//DYmumu_Mcut200-MC_31X_V3//PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_8.root",
#     "/store/user/aeverett//CMSSW_3_2_1//TTbar_Tauola//aeverett//TTbar_Tauola_CMSSW_3_2_1_step1//TTbar_Tauola_CMSSW_3_2_1_step1//b947059661a5e0b4111d8e6607110054//step1_8.root"
#     '/store/relval/CMSSW_3_3_0/RelValInclusiveppMuX/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_31X_V9-v1/0002/889D6338-C4B7-DE11-BE10-001A928116FC.root',
#     '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_15_20/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v1/0001/0847D5EB-5CB7-DE11-BCF1-0017312B5F3F.root',
#     '/store/relval/CMSSW_3_3_2/RelValProdTTbar/GEN-SIM-RAW/MC_31X_V9-v1/0004/0835072A-0DC7-DE11-A39D-001D09F2441B.root'
      '/store/user/aeverett/Commisioning/MC900//PurdueSkim2///goodGlobalMuons_24.root','/store/user/aeverett/Commisioning/MC900//PurdueSkim2///goodGlobalMuons_155.root','/store/user/aeverett/Commisioning/MC900//PurdueSkim2///goodGlobalMuons_126.root',

    )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MRTU_Reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'GR09_R_34X_V2::All'
process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)




# Schedule definition
#process.schedule = cms.Schedule()
process.schedule = cms.Schedule(process.L1simulation_step)

##
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
##
##//##
## process.load("UserCode.GlbSelectorStudy.glbselectorstudy_cff")
## process.glbSelStudy.doAssoc = False
## process.glbSelStudy.trkMuAssocLabel = "tpToTkmuTrackAssociation"
## process.glbSelStudy.staMuAssocLabel = "tpToStaTrackAssociation"
## process.glbSelStudy.glbMuAssocLabel = "tpToGlbTrackAssociation"
## process.glbSelStudy.trackProducer = "globalMuons"
## process.glbSelStudy.trackAssociator = "TrackAssociatorByDeltaR"
## process.pSelStudy = cms.Path(process.muonAssociation_seq*process.glbSelStudy)
## process.TFileService = cms.Service("TFileService", fileName = cms.string("TFS_selStudy.root"))
process.validation_step = cms.Path(process.validation)
##//##
process.schedule.extend([process.raw2digi_step,process.reconstruction_step])
process.schedule.extend([process.validation_step])
process.schedule.extend([process.endjob_step,process.out_step])


# Automatic addition of the customisation function
def customise(process):
     from Workspace.MuonRecoTreeUtility.muonRecoTreeUtility_cff import insertMRTU
     insertMRTU(process)
     return (process)


#process = customise(process)
