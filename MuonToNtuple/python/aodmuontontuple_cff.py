import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi import *
ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
ttrhbwr.ComputeCoarseLocalPositionFromDisk = True
#from TrackingTools.TrackRefitter.TracksToTrajectories_cff import *

from UserCode.MuonToNtuple.aodmuontontuple_cfi import *
