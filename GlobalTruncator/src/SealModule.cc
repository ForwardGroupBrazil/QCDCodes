#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "UserCode/GlobalTruncator/src/GlobalTruncator.h"


DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(GlobalTruncator);

