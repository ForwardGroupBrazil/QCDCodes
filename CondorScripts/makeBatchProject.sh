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
que=$8

let nFiles=$lastFile-$firstFile+1

if [ "$?" != "0" ]; then
  echo "Failed to initialize CMSSW environment with scram in $CMSSW_HOME."
  exit 1
fi


#
# environment setup
#
ORIGINAL_DIR=`pwd`
#PATH=$PATH:$ORIGINAL_DIR
#export PATH
#cd $CMSSW_HOME/src
#eval `scramv1 runtime -sh`
#
#if [ "$?" != "0" ]; then
#  echo "Failed to initialize CMSSW environment with scram in $CMSSW_HOME."
#  exit 1
#fi


#
# name project common directory
#
projectName=`basename $configTemplate Template.cfg`
blockName=`basename $configFileBlock .plain`

scratch_dir="${HOME}/w0"

if ! [ -d $scratch_dir ]; then
  scratch_dir="/scratch"
fi
if ! [ -d $scratch_dir ]; then
  scratch_dir="/tmp"
fi

COMMON_HOME=${scratch_dir}/${jobname}
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

echo "Done making CFG"

declare -i job=1
while [ $job -le $totalJobs ]; do
    
#
# name the batch run script
#
    submitFile=$COMMON_DIR/${blockName}_${job}
    rm $submitFile
#
#
#
    ${ORIGINAL_DIR}/makeBatchSubmit.sh $jobName $blockName $job $COMMON_DIR
    chmod +x $submitFile
    
#
# Submit the project
#
    bsub -q$que $submitFile
#-    rm $submitFile

    job=$(($job+1))

done

echo -n "Jobs for $jobName are created in " $COMMON_DIR
echo ""
pwd
cd $ORIGINAL_DIR


