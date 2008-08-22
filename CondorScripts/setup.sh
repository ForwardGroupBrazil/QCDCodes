#!/bin/bash -x
CMSSW_VERSION=$1

source /opt/osg/setup.sh
#MCMSSW_BASE=/grp/cms/purdueCMS
#. ${MCMSSW_BASE}/cmsset_default.sh

#export SCRAM_ARCH=slc3_ia32_gcc323
export SCRAM_ARCH=slc4_ia32_gcc345

export MCMSSW_BASE=/apps/02/cmssoft/cms/
#MCMSSW_BASE=/grp/cms/purdueCMS
source ${MCMSSW_BASE}/cmsset_default.sh

source /grp/cms/crab/crab.sh

#. /grp/cms/crab/CRAB_1_5_3_pre3/crab.sh



voms-proxy-init -hours 192 -voms cms
#export X509_USER_PROXY=$(voms-proxy-info -path)
echo "here"
X509_USER_PROXY=/tmp/x509up_u$UID
export X509_USER_PROXY
echo "there $X509_USER_PROXY $CMSSW_VERSION"

#cd $RCAC_SCRATCH/${CMSSW_VERSION}/src
cd /scratch/scratch96/a/aeverett/${CMSSW_VERSION}/src
echo "here"
pwd
#cd /scratch/osg/cmssw/xu2/${CMSSW_VERSION}/src
eval `scramv1 runtime -sh`
export PATH=/scratch/scratch96/a/aeverett/reReco_2011:/grp/cms/tools/prod:/opt/d-cache/dcap/bin:$PATH
#cd ../..

EDG_WL_LOCATION=/opt/osg/glite
export EDG_WL_LOCATION

#LD_LIBRARY_PATH="${CMS_PATH}/${SCRAM_ARCH}/external/boost/1.33.1/lib/:${LD_LIBRARY_PATH}"
#LD_LIBRARY_PATH="${CMS_PATH}/${SCRAM_ARCH}/external/python/2.4.2/lib/:${LD_LIBRARY_PATH}"
#LD_LIBRARY_PATH="${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/${CMSSW_VERSION}/lib/${SCRAM_ARCH}:${LD_LIBRARY_PATH}"
#LD_LIBRARY_PATH="${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/${CMSSW_VERSION}/external/${SCRAM_ARCH}/lib:${LD_LIBRARY_PATH}"

#export LD_LIBRARY_PATH

#CMSSW_SEARCH_PATH="${CMS_PATH}/${SCRAM_ARCH}/cms/cmssw/${CMSSW_VERSION}/share:${CMSSW_SEARCH_PATH}"
#export CMSSW_SEARCH_PATH

