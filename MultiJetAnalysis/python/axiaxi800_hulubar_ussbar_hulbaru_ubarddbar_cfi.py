import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source ('PoolSource',fileNames = readFiles)
readFiles.extend((
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_10_2_BqW.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_1_1_sh0.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_2_1_LCM.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_4_1_t5G.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_5_1_8H9.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_6_1_M1v.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_7_1_1C1.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_8_1_MDv.root',
'/store/user/kkousour/2011/Jets/mc/Axigluon/axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_9_1_Pa4.root'
)
)
