import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source ('PoolSource',fileNames = readFiles)
readFiles.extend( (
'/store/user/kkousour/2011/Jets/JetAOD_84_1_VwO.root',
'/store/user/kkousour/2011/Jets/JetAOD_54_1_0B2.root',
'/store/user/kkousour/2011/Jets/JetAOD_50_1_oVM.root',
'/store/user/kkousour/2011/Jets/JetAOD_34_1_pba.root',
'/store/user/kkousour/2011/Jets/JetAOD_6_1_aI8.root',
'/store/user/kkousour/2011/Jets/JetAOD_56_1_zti.root',
'/store/user/kkousour/2011/Jets/JetAOD_74_1_K7s.root',
'/store/user/kkousour/2011/Jets/JetAOD_10_1_WAr.root',
'/store/user/kkousour/2011/Jets/JetAOD_14_1_cJ5.root',
'/store/user/kkousour/2011/Jets/JetAOD_68_1_1SN.root',
'/store/user/kkousour/2011/Jets/JetAOD_58_1_WPU.root',
'/store/user/kkousour/2011/Jets/JetAOD_25_1_cqo.root',
'/store/user/kkousour/2011/Jets/JetAOD_21_1_OfV.root',
'/store/user/kkousour/2011/Jets/JetAOD_48_1_Lac.root',
'/store/user/kkousour/2011/Jets/JetAOD_55_1_3ac.root',
'/store/user/kkousour/2011/Jets/JetAOD_47_1_tu5.root',
'/store/user/kkousour/2011/Jets/JetAOD_66_1_ypY.root',
'/store/user/kkousour/2011/Jets/JetAOD_62_1_2HF.root',
'/store/user/kkousour/2011/Jets/JetAOD_51_1_UP8.root',
'/store/user/kkousour/2011/Jets/JetAOD_30_1_AfL.root',
'/store/user/kkousour/2011/Jets/JetAOD_2_1_l1V.root',
'/store/user/kkousour/2011/Jets/JetAOD_90_1_yog.root',
'/store/user/kkousour/2011/Jets/JetAOD_88_1_d2r.root',
'/store/user/kkousour/2011/Jets/JetAOD_67_1_MWX.root',
'/store/user/kkousour/2011/Jets/JetAOD_20_1_7QQ.root',
'/store/user/kkousour/2011/Jets/JetAOD_27_1_l4o.root',
'/store/user/kkousour/2011/Jets/JetAOD_70_1_0VC.root',
'/store/user/kkousour/2011/Jets/JetAOD_4_1_N2v.root',
'/store/user/kkousour/2011/Jets/JetAOD_91_1_Jwt.root',
'/store/user/kkousour/2011/Jets/JetAOD_61_1_FSa.root',
'/store/user/kkousour/2011/Jets/JetAOD_71_1_wvU.root',
'/store/user/kkousour/2011/Jets/JetAOD_33_1_uw5.root',
'/store/user/kkousour/2011/Jets/JetAOD_11_1_Cw2.root',
'/store/user/kkousour/2011/Jets/JetAOD_85_1_TFw.root',
'/store/user/kkousour/2011/Jets/JetAOD_57_1_uKC.root',
'/store/user/kkousour/2011/Jets/JetAOD_35_1_fc2.root',
'/store/user/kkousour/2011/Jets/JetAOD_3_1_64R.root',
'/store/user/kkousour/2011/Jets/JetAOD_52_1_Kka.root',
'/store/user/kkousour/2011/Jets/JetAOD_45_1_E0b.root',
'/store/user/kkousour/2011/Jets/JetAOD_13_1_5UK.root',
'/store/user/kkousour/2011/Jets/JetAOD_78_1_1LW.root',
'/store/user/kkousour/2011/Jets/JetAOD_37_1_5iO.root',
'/store/user/kkousour/2011/Jets/JetAOD_32_1_MhV.root',
'/store/user/kkousour/2011/Jets/JetAOD_41_1_ND9.root',
'/store/user/kkousour/2011/Jets/JetAOD_59_1_ZGL.root',
'/store/user/kkousour/2011/Jets/JetAOD_8_1_1D9.root',
'/store/user/kkousour/2011/Jets/JetAOD_18_1_fEi.root',
'/store/user/kkousour/2011/Jets/JetAOD_82_1_aqq.root',
'/store/user/kkousour/2011/Jets/JetAOD_63_1_wWj.root',
'/store/user/kkousour/2011/Jets/JetAOD_15_1_hCO.root',
'/store/user/kkousour/2011/Jets/JetAOD_9_1_Bh3.root',
'/store/user/kkousour/2011/Jets/JetAOD_64_1_AwU.root',
'/store/user/kkousour/2011/Jets/JetAOD_5_1_I3y.root',
'/store/user/kkousour/2011/Jets/JetAOD_24_1_3dO.root',
'/store/user/kkousour/2011/Jets/JetAOD_22_1_UGn.root',
'/store/user/kkousour/2011/Jets/JetAOD_73_1_3jS.root',
'/store/user/kkousour/2011/Jets/JetAOD_44_1_qRa.root',
'/store/user/kkousour/2011/Jets/JetAOD_83_1_roi.root',
'/store/user/kkousour/2011/Jets/JetAOD_72_1_CtI.root',
'/store/user/kkousour/2011/Jets/JetAOD_79_1_JrU.root',
'/store/user/kkousour/2011/Jets/JetAOD_75_1_hT4.root',
'/store/user/kkousour/2011/Jets/JetAOD_60_1_30f.root',
'/store/user/kkousour/2011/Jets/JetAOD_87_1_acm.root',
'/store/user/kkousour/2011/Jets/JetAOD_53_1_Byn.root',
'/store/user/kkousour/2011/Jets/JetAOD_65_1_Tqu.root',
'/store/user/kkousour/2011/Jets/JetAOD_80_1_g0K.root',
'/store/user/kkousour/2011/Jets/JetAOD_81_1_XyQ.root',
'/store/user/kkousour/2011/Jets/JetAOD_93_1_1d1.root',
'/store/user/kkousour/2011/Jets/JetAOD_17_1_CGY.root',
'/store/user/kkousour/2011/Jets/JetAOD_1_1_Aqk.root',
'/store/user/kkousour/2011/Jets/JetAOD_69_1_tpJ.root',
'/store/user/kkousour/2011/Jets/JetAOD_76_1_gP7.root',
'/store/user/kkousour/2011/Jets/JetAOD_46_1_rga.root',
'/store/user/kkousour/2011/Jets/JetAOD_43_1_DRi.root',
'/store/user/kkousour/2011/Jets/JetAOD_36_1_hK3.root',
'/store/user/kkousour/2011/Jets/JetAOD_38_1_Pjj.root',
'/store/user/kkousour/2011/Jets/JetAOD_39_1_bxU.root',
'/store/user/kkousour/2011/Jets/JetAOD_26_1_K1V.root',
'/store/user/kkousour/2011/Jets/JetAOD_92_1_EU9.root',
'/store/user/kkousour/2011/Jets/JetAOD_49_1_BxO.root',
'/store/user/kkousour/2011/Jets/JetAOD_89_1_hYJ.root',
'/store/user/kkousour/2011/Jets/JetAOD_86_1_Evs.root',
'/store/user/kkousour/2011/Jets/JetAOD_7_1_Xep.root',
'/store/user/kkousour/2011/Jets/JetAOD_40_1_jck.root',
'/store/user/kkousour/2011/Jets/JetAOD_23_1_YvN.root',
'/store/user/kkousour/2011/Jets/JetAOD_12_1_lH6.root',
'/store/user/kkousour/2011/Jets/JetAOD_19_1_PqW.root',
'/store/user/kkousour/2011/Jets/JetAOD_31_1_h8x.root',
'/store/user/kkousour/2011/Jets/JetAOD_42_1_fDX.root',
'/store/user/kkousour/2011/Jets/JetAOD_16_1_dWk.root',
'/store/user/kkousour/2011/Jets/JetAOD_77_1_Ha8.root',
'/store/user/kkousour/2011/Jets/JetAOD_29_1_zGa.root',
'/store/user/kkousour/2011/Jets/JetAOD_28_1_xdu.root'
))