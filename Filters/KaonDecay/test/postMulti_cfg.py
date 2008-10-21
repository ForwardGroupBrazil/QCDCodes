import FWCore.ParameterSet.Config as cms

process = cms.Process("postMulti")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")

#process.load("Configuration.EventContent.EventContent_cff")

process.load("DQMServices.Components.EDMtoMEConverter_cff")

process.load("DQMServices.Components.DQMEnvironment_cfi")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    'file:output_edm_1.root',
    'file:output_edm_2.root',
    'file:output_edm_3.root',
    'file:output_edm_4.root',
    'file:output_edm_5.root',
    'file:output_edm_6.root',
    'file:output_edm_7.root',
    'file:output_edm_8.root',
    'file:output_edm_9.root',
    'file:output_edm_10.root',
    'file:output_edm_11.root'
    )
                            )

process.PoolSource.fileNames = ['/store/user/hyxu/SingleKMinusPt2_200-newStep3/newStep3-SingleKMinusPt2_200-0062.root',]

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_*_*_validReco'),
    fileName = cms.untracked.string('postMulti.root')
)

#process.MyOut = cms.OutputModule("PoolOutputModule",
#    process.FEVTSIMHLTDEBUGEventContent,
#    fileName = cms.untracked.string('adam2.root')
#)

process.postProcessor = cms.EDFilter("PostProcessor",
    commands = cms.vstring("RecoMuonV/MultiTrack E v_Efficiency_eta \'Efficiency Plot vs #eta\' \'num_assoc(simToReco)_eta\' \'num_simul_eta\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_phi \'Efficiency Plot vs #eta\' \'num_assoc(simToReco)_phi\' \'num_simul_phi\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_pt \'Efficiency Plot vs pt\' \'num_assoc(simToReco)_pT\' \'num_simul_pT\'", 
        "RecoMuonV/MultiTrack E v_Efficiency_hit \'Efficiency Plot vs hits\' \'num_assoc(simToReco)_hit\' \'num_simul_hit\'", 
        "RecoMuonV/MultiTrack E v_Fake_eta \'Fake Plot vs #eta\' \'num_assoc(recoToSim)_eta\' \'num_reco_eta\'", 
        "RecoMuonV/MultiTrack E v_Fake_phi \'Fake Plot vs #eta\' \'num_assoc(recoToSim)_phi\' \'num_reco_phi\'", 
        "RecoMuonV/MultiTrack E v_Fake_pt \'Fake Plot vs pt\' \'num_assoc(recoToSim)_pT\' \'num_reco_pT\'", 
        "RecoMuonV/MultiTrack E v_Fake_hit \'Fake Plot vs hits\' \'num_assoc(recoToSim)_hit\' \'num_reco_hit\'", 
        "RecoMuonV/MultiTrack R v_dxy_Pull_eta \'Pull(dxy) vs #eta\' \'dxypull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dz_Pull_eta \'Pull(dz) vs #eta\' \'dzpull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_phi_Pull_eta \'Pull(#phi) vs #eta\' \'phipull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Pull_eta \'Pull(pt) vs #eta\' \'ptpull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_theta_Pull_eta \'Pull(#theta) vs #eta\' \'thetapull_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Pull_phi \'Pull(pt) vs #phi\' \'ptpull_vs_phi\'", 
        "RecoMuonV/MultiTrack R v_phi_Pull_phi \'Pull(#phi) vs #phi\' \'phipull_vs_phi\'", 
        "RecoMuonV/MultiTrack R v_theta_Pull_phi \'Pull(#theta) vs #phi\' \'thetapull_vs_phi\'", 
        "RecoMuonV/MultiTrack R v_dxy_Residual_eta \'Res(dxy) vs #eta\' \'dxyres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Resolution_eta \'Res(pt) vs #eta\' \'ptres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_dz_Residual_eta \'Res(dz) vs #eta\' \'dzres_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_phi_Residual_eta \'Res(#phi) vs #eta\' \'phires_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_cotTheta_Residual_eta \'Res(cot(#theta)) vs #eta\' \'cotThetares_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_eta_Residual_eta \'Res(#eta) vs #eta\' \'etares_vs_eta\'", 
        "RecoMuonV/MultiTrack R v_Pt_Resolution_phi \'Res(p_{t}) vs #phi\' \'ptres_vs_phi\'", 
        "RecoMuonV/MultiTrack R v_phi_Residual_phi \'Res(#phi) vs #phi\' \'phires_vs_phi\'", 
        "RecoMuonV/MultiTrack R v_dxy_Residual_pt \'Res(dxy vs p_{t}\' \'dxyres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_Pt_Resolution_pt \'Res(p_{t}) vs pt\' \'ptres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_dz_Residual_pt \'Res(dz) vs pt\' \'dzres_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_phi_Residual_pt \'Res(#phi) vs pt\' \'phires_vs_pt\'", 
        "RecoMuonV/MultiTrack R v_cotTheta_Residual_pt \'Res(cot(#theta)) vs pt\' \'cotThetares_vs_pt\'", 
        "RecoMuonV/MultiTrack C v_Efficiency_glb_tk \'Efficiency Plot\' \'globalMuons_tpToGlbAssociation/num_assoc(simToReco)_eta\' \'general_trackingParticleRecoAsssociation/num_assoc(simToReco)_eta\'", 
        "RecoMuonV/MultiTrack C v_Efficiency_glb_sta \'Efficiency Plot\' \'globalMuons_tpToGlbAssociation/num_assoc(simToReco)_eta\' \'standAloneMuons_UpdatedAtVtx_tpToStaAssociation/num_assoc(simToReco)_eta\'", 
        "RecoMuonV/MultiTrack C v_Efficiency_glb_tk \'Efficiency Plot\' \'globalMuons_tpToGlbAssociationPos/num_assoc(simToReco)_eta\' \'general_trackingParticleRecoAsssociationPos/num_assoc(simToReco)_eta\'", 
        "RecoMuonV/MultiTrack C v_Efficiency_glb_sta \'Efficiency Plot\' \'globalMuons_tpToGlbAssociationPos/num_assoc(simToReco)_eta\' \'standAloneMuons_UpdatedAtVtx_tpToStaAssociationPos/num_assoc(simToReco)_eta\'"),
    subDir = cms.untracked.string('RecoMuonV/MultiTrack'),
    outputFileName = cms.untracked.string('postMulti.root')
)

process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)
#process.outpath = cms.EndPath(process.MyOut)
process.p1 = cms.Path(process.EDMtoMEConverter*process.postProcessor)
process.schedule = cms.Schedule(process.p1)

process.MessageLogger.categories = []
process.MessageLogger.FrameworkJobReport.FwkJob.limit = 0
process.MessageLogger.cerr.FwkReport.limit = 0
process.EDMtoMEConverter.convertOnEndRun = True
process.DQMStore.referenceFileName = ''
process.dqmSaver.convention = 'RelVal'
process.dqmSaver.workflow = '/Validation/Test/RECO'


