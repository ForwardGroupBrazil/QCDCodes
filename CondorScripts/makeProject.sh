#!/bin/bash
#
# ---------------------------------------------------------
#
# Initialize arguments:
jobName=$1
declare -i firstFile=$2
declare -i lastFile=$3
declare -i filesPerJob=$4
configTemplate=$5
configFileBlock=$6
que=$7

let nFiles=$lastFile-$firstFile+1

echo "CMSSW_BASE $CMSSW_BASE"

#
# environment setup
#
ORIGINAL_DIR=`pwd`
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`

if [ "$?" != "0" ]; then
  echo "Failed to initialize CMSSW environment with scram in $CMSSW_BASE."
  exit 1
fi


#
# name project common directory
#
blockName=`basename $configFileBlock .txt`

scratch_dir="${HOME}/scratch0"

if ! [ -d $scratch_dir ]; then
  scratch_dir="/scratch"
fi
if ! [ -d $scratch_dir ]; then
  scratch_dir="/tmp"
fi

COMMON_DIR=${scratch_dir}/${jobName}/${blockName}
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

#
# Submit the project
#

if [ "$subtype" = "PBS" ] || [ "$subtype" = "BSUB" ]; then 
    declare -i job=1
    while [ $job -le $totalJobs ]; do	
	submitFile=batchRunScript.sh
	fileName=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d", $1, $2, $3)}'`
	$subcommand -q$que $submitFile $fileName $COMMON_DIR
	job=$(($job+1))	
    done
elif [ "$subtype" = "CONDOR" ]; then
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
    
    condor_submit $submitFile
else
    echo "Unknow submission type: " $subtype    
fi

echo -n "Jobs for $jobName are created in " $COMMON_DIR
echo ""
pwd
cd $ORIGINAL_DIR


