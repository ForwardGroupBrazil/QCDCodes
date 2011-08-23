#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "KKousour/QCDAnalysis/plugins/MultijetFilter.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

using namespace reco;

typedef MultijetFilter<CaloJet> MultiCaloJetFilter;
DEFINE_FWK_MODULE(MultiCaloJetFilter);

typedef MultijetFilter<PFJet> MultiPFJetFilter;
DEFINE_FWK_MODULE(MultiPFJetFilter);

typedef MultijetFilter<GenJet> MultiGenJetFilter;
DEFINE_FWK_MODULE(MultiGenJetFilter);


