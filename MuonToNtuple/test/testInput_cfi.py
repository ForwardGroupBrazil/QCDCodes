import FWCore.ParameterSet.Config as cms

# from 

def testInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(200),
                   debugFlag = cms.untracked.bool(True),
                   
                   fileNames = cms.untracked.vstring(
     '/store/relval/CMSSW_3_1_0/RelValSingleMuPt100/GEN-SIM-RECO/MC_31X_V1-v1/0001/F671A054-6166-DE11-B0C7-001D09F2543D.root',
#     '/store/relval/CMSSW_2_2_9/RelValZMM/GEN-SIM-RECO/STARTUP_V11_v1/0001/CAB408F0-F831-DE11-8E1F-001617E30F50.root',
#     '/store/relval/CMSSW_2_2_9/RelValZMM/GEN-SIM-RECO/STARTUP_V11_v1/0001/80C561B5-FA31-DE11-9913-001617DC1F70.root',
#     '/store/mc/Summer08/Zmumu/GEN-SIM-RECO/IDEAL_V11_redigi_v2/0010/D4314D3D-0201-DE11-B0E9-00E08179177F.root',
		#	'/store/relval/CMSSW_2_1_7/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/F484E210-FA7D-DD11-A42A-000423D94534.root'
#     'file:/work/hdyoo/Work/DYanalysis/CMSSW_2_2_1/src/PYTHIA6_DYmumu_M6_40_filter_10TeV_cff_redigi_RECO.root'
       		   )
                 )
