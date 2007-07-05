#!/bin/bash
#
# ---------------------------------------------------------
#
# Initialize arguments:
jobName=$1
declare -i firstFile=$2
declare -i lastFile=$3
declare -i filesPerJob=$4
CMSSW_HOME=$5
configTemplate=$6
configFileBlock=$7

let nFiles=$lastFile-$firstFile+1

if [ "$?" != "0" ]; then
  echo "Failed to initialize CMSSW environment with scram in $CMSSW_HOME."
  exit 1
fi


#
# environment setup
#
ORIGINAL_DIR=`pwd`
PATH=$PATH:$ORIGINAL_DIR
export PATH
cd $CMSSW_HOME/src
eval `scramv1 runtime -sh`

if [ "$?" != "0" ]; then
  echo "Failed to initialize CMSSW environment with scram in $CMSSW_HOME."
  exit 1
fi


#
# name project common directory
#
projectName=`basename $configTemplate Template.cfg`
blockName=`basename $configFileBlock .plain`

scratch_dir="/grp/cms/users/aeverett/scratch"

if ! [ -d $scratch_dir ]; then
  scratch_dir="/scratch"
fi
if ! [ -d $scratch_dir ]; then
  scratch_dir="/tmp"
fi

COMMON_HOME=${scratch_dir}/${USER}/${jobname}
echo "COMMON_HOME " $COMMON_HOME

COMMON_DIR=

COMMON_DIR=${COMMON_DIR:-${COMMON_HOME}/$jobName-$projectName/$blockName}
echo "COMMON_DIR " $COMMON_DIR


#
# make the project common directory
#
mkdir -p $COMMON_DIR
cd $COMMON_DIR || die "Failed to create directory $COMMON_DIR"
cd $ORIGINAL_DIR


#
# make the CFG files
#
totalJobs=`${ORIGINAL_DIR}/makeCfg.sh $jobName $configTemplate $configFileBlock $firstFile $lastFile $filesPerJob $COMMON_DIR`


#
# name the condor submit script
#
submitFile=$COMMON_DIR/submit

#
# add jobs to the condor submit script
#
${ORIGINAL_DIR}/makeCondorSubmit.sh $jobName $blockName $totalJobs $COMMON_DIR

#
# copy the run script to the common directory
#
cp condorRunScript.csh ${COMMON_DIR}/.
chmod +x ${COMMON_DIR}/condorRunScript.csh

#
# Submit the project
#
condor_submit $submitFile

echo -n "Jobs for $jobName are created in " $COMMON_DIR
echo ""
pwd
cd $ORIGINAL_DIR


