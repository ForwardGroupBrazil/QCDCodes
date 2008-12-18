#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "UserCode/GlobalTruncator/src/GlobalTruncator.h"
#include "UserCode/GlobalTruncator/src/TruncAnalyzer.h"


DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(GlobalTruncator);
DEFINE_FWK_MODULE(TruncAnalyzer);
