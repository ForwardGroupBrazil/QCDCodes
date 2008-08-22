#!/bin/bash 
#
# Setup:

#
# Usage:
#       farmoutAnalysisJobs <jobName> <CMSSW Version> <config file>
#
# The config file should refer to the following macros, which are automatically
# inserted by this script:
#
# $inputFileNames     ==>  Will be replaced by list of input files
# $outputFileName     ==>  Will be replaced by the $inputFileName-output.root
#
# Job parameters
#


# Initialize default settings:

# for storing output
SRM_SERVER=srm://dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=
PNFS_HOME=/store/user/aeverett/testCONDOR01

# for getting input
DCAP_SERVER=dcap://dcache.rcac.purdue.edu:22125
# path to directory containing lfns: /store/...
#PNFS_STORE=/store/user/hyxu
PNFS_STORE=/store/user/aeverett/testCONDOR01

SRM_HOME=${SRM_SERVER}/${PNFS_HOME}
DCAP_HOME=${DCAP_SERVER}${PNFS_STORE}

SITE_REQUIREMENTS='TARGET.HasAFS_OSG =?= True && TARGET.OSRedHatRelease =!= "Scientific Linux SL Release 3.0.4 (SL)"'
DISK_REQUIREMENTS=2000
MEMORY_REQUIREMENTS=1200
MIN_PROXY_HOURS=24

# special exit status to force job to leave the queue
FAIL_JOB=42

scratch_dir="/data"
if ! [ -d $scratch_dir ]; then
  scratch_dir="/scratch"
fi
if ! [ -d $scratch_dir ]; then
  scratch_dir="/tmp"
fi
SUBMIT_HOME=`pwd`

basename() {
  # This shell function is faster than calling /bin/basename
  path=$1
  suffix=$2
  path=${path##*/}  # get everything after the final '/'
  if [ ! -z $suffix ]; then
    path=${path%$suffix}
  fi
  echo $path
}

realpath() {
  readlink -f $1
}

logerror() {
  echo 2>&1 "$@"
}

die() {
  logerror
  logerror "$@"
  exit 1
}

outputFileExists() {
  fname=$1

  #Strip off srm://hostname:8443 to get raw pnfs path.
  local_fname=${fname#srm://*:8443}

  if [ -f "$local_fname" ]; then
    return 0
  fi
  return 1
}

check_proxy() {
  hours=$1
  proxy=$2
  if ! [ -f "$proxy" ]; then
    logerror
    logerror "NOTE: No grid proxy found.  (Expecting to find it here: $proxy.)"
    return 1
  fi

  #Issue a warning if less than this many seconds remain:
  min_proxy_lifetime=$((3600*$hours))

  seconds_left="`voms-proxy-info --timeleft --file=$proxy 2>/dev/null`"

  if [ "$seconds_left" = "" ]; then
    echo "WARNING: cannot find time remaining for grid proxy."
    voms-proxy-info -all -path $proxy
    return 0
  fi
  if [ "$seconds_left" -lt "$min_proxy_lifetime" ]; then
    logerror
    logerror "NOTE: grid proxy is about to expire:"
    logerror "voms-proxy-info"
    voms-proxy-info --file=$proxy
    return 1
  fi

}

PrintUsage() {
  echo "USAGE: farmoutAnalysisJobs [options] <jobName> <CMSSW Path> <config file>"
  echo ""
  echo "OPTIONS:"
  echo "  --skip-srmcp=0"
  echo "  --output-dir=${SRM_HOME}/<jobName>-<runName>"
  echo "  --input-dir=${DCAP_HOME}/<jobName>"
  echo "  --input-dbs-path=/dbs/dataset/path   (in place of --input-dir)"
  echo "  --submit-dir=${SUBMIT_HOME}/<jobName>-<runName>"
  echo "  --no-submit"
  echo "  --job-count=N            (limit the number of jobs that are created)"
  echo "  --input-files-per-job=1"
  echo "  --skip-existing-output   (do not create jobs if output file exists)"
  echo "  --skip-existing-jobs     (do not create jobs if job already created)"
  echo "  --match-input-files='*.root'"
  echo "  --exclude-input-files=pattern"
  echo "  --input-file-list=file"
  echo "  --memory-requirements=$MEMORY_REQUIREMENTS (megabytes)"
  echo "  --disk-requirements=$DISK_REQUIREMENTS   (megabytes)"
  echo "  --site-requirements=$SITE_REQUIREMENTS"
  echo "  --save-failed-datafiles  (save root file from failed cmsRun job)"
  echo "                           (in <output-dir>-cmsRun-failed)"
  echo "  --save-missing-input-file-list=file"
  echo "  --assume-input-files-exist (do not check file existence at submit time)"
  echo ""
  echo "Note that <runName> is taken from the name of the config file."
  exit 2
}

OPTS=`getopt -o "h" -l "help,output-dir:,input-dir:,submit-dir:,no-submit,job-count:,skip-existing-output,skip-srmcp,skip-existing-jobs,match-input-files:,exclude-input-files:,input-files-per-job:,disk-requirements:,memory-requirements:,input-file-list:,input-dbs-path:,save-failed-datafiles,save-missing-input-file-list:,assume-input-files-exist,site-requirements:" -- "$@"`
if [ $? -ne 0 ]; then PrintUsage; fi

eval set -- "$OPTS"

NO_SUBMIT=
JOB_LIMIT=
SKIP_EXISTING_OUTPUT=
SKIP_EXISTING_JOBS=
OUTPUT_DIR=
SKIP_SRMCP=0
INPUT_DIR=
SUBMIT_DIR=
MATCH_INPUT_FILES='*.root'
EXCLUDE_INPUT_FILES=
INPUT_FILES_PER_JOB=1
INPUT_FILE_LIST=
INPUT_DBS_PATH=
SAVE_FAILED_DATAFILES=
SAVE_MISSING_INPUT_FILE_LIST=
ASSUME_INPUT_FILES_EXIST=

while [ ! -z "$1" ]
do
  case "$1" in
    -h) PrintUsage;;
    --help) PrintUsage;;
    --no-submit) NO_SUBMIT=1;;
    --skip-srmcp) SKIP_SRMCP=1;;
    --job-count) shift; JOB_LIMIT=$1;;
    --skip-existing-output) SKIP_EXISTING_OUTPUT=1;;
    --skip-existing-jobs) SKIP_EXISTING_JOBS=1;;
    --output-dir) shift; OUTPUT_DIR=$1;;
    --input-dir) shift; INPUT_DIR=$1;;
    --submit-dir) shift; SUBMIT_DIR=$1;;
    --match-input-files) shift; MATCH_INPUT_FILES=$1;;
    --exclude-input-files) shift; EXCLUDE_INPUT_FILES=$1;;
    --input-files-per-job) shift; INPUT_FILES_PER_JOB=$1;;
    --disk-requirements) shift; DISK_REQUIREMENTS=$1;;
    --memory-requirements) shift; MEMORY_REQUIREMENTS=$1;;
    --input-file-list) shift; INPUT_FILE_LIST=$1;;
    --input-dbs-path) shift; INPUT_DBS_PATH=$1;;
    --save-failed-datafiles) SAVE_FAILED_DATAFILES=1;;
    --save-missing-input-file-list) shift; SAVE_MISSING_INPUT_FILE_LIST=$1;;
    --assume-input-files-exist) ASSUME_INPUT_FILES_EXIST=0;;
    --site-requirements) shift; SITE_REQUIREMENTS="$1";;
    --) shift; break;;
    *) die "Unexpected option $1";;
  esac
  shift
done

if [ "$#" -ne 3 ]; then PrintUsage; fi


# Check for some required utilities
for exe in scramv1 condor_submit cmsRun.sh readlink voms-proxy-info; do
  if ! which $exe >& /dev/null; then
    die "Cannot find $exe in PATH.  Your environment is not correctly set up."
  fi
done

# Additional command-line arguments

jobName=$1
#CMSSW_HOME=`realpath $2`
CMSSW_HOME=$2
configTemplate=`realpath $3`

# Now we have all the user's input.


runName=`basename $configTemplate .cfg`
echo "INPUT_DIR=" $INPUT_DIR
if [ "$INPUT_DIR" = "" ]; then
    if [ "$INPUT_FILE_LIST" != "" ] || [ "$INPUT_DBS_PATH" != "" ]; then
        INPUT_DIR=${DCAP_SERVER}${PNFS_STORE}
    fi
fi

OUTPUT_DIR=${OUTPUT_DIR:-${SRM_HOME}/$jobName-$runName}
INPUT_DIR=${INPUT_DIR:-${DCAP_HOME}/$jobName}
SUBMIT_DIR=${SUBMIT_DIR:-${SUBMIT_HOME}/$jobName-$runName}

#Strip off dcap://hostname:22125 to get raw pnfs path.
LOCAL_INPUT_DIR=${INPUT_DIR#dcap://*:22125}
echo $LOCAL_INPUT_DIR

#Get the part of INPUT_DIR that was stripped off.
INPUT_BASE=${INPUT_DIR%$LOCAL_INPUT_DIR}

if [ "$INPUT_FILE_LIST" != "" ]; then
    if ! [ -f "$INPUT_FILE_LIST" ]; then
        die "Error: No such file: $INPUT_FILE_LIST"
    fi
    INPUT_FILE_LIST=`realpath $INPUT_FILE_LIST`
#elif ! [ -d "$LOCAL_INPUT_DIR" ]; then
#  die "Error: No such input directory: $LOCAL_INPUT_DIR"
fi

if [ "$INPUT_DBS_PATH" != "" ]; then
    if [ "$DBSCMD_HOME" = "" ] || ! [ -d "$DBSCMD_HOME" ]; then
      die "DBS client is not correctly set up: DBSCMD_HOME is invalid"
    fi
fi

#if ! [ -d "$CMSSW_HOME" ]; then
#  die "Error: No such CMSSW directory: $CMSSW_HOME"
#fi

if [ -d "$SUBMIT_DIR" ] && [ "$SKIP_EXISTING_JOBS" != "1" ]; then
  logerror
  logerror "Error: Submit directory already exists: $SUBMIT_DIR"
  logerror
  logerror "You must either remove it, or specify --skip-existing-jobs, or"
  logerror "specify a different job name or submission directory with --submit-dir"
  exit 1
fi

proxy=${X509_USER_PROXY:-/tmp/x509up_u$UID}

if [ "$NO_SUBMIT" != 1 ] && ! check_proxy $MIN_PROXY_HOURS $proxy; then
  logerror
  logerror "Either rerun this command with --no-submit or create a new grid proxy"
  logerror "and rerun this command.  Example of how to create a grid proxy:"
  logerror
  logerror "voms-proxy-init --voms=cms --hours=48"
  exit 1
fi

# Check the config template

for macro in \$inputFileNames \$outputFileName; do
  if ! grep -F -q $macro $configTemplate; then
    die "$macro must appear on the configuration template.  I can't find it in $configTemplate"
  fi
done

#
# CMSSW environment setup
#
originalDir=`pwd`
PATH=$PATH:$originalDir
export PATH
cd $CMSSW_HOME || die "Failed to cd to $CMSSW_HOME."
eval `scramv1 runtime -sh`

if [ "$?" != "0" ]; then
  die "Failed to initialize CMSSW environment with scram in $CMSSW_HOME."
fi

runDir=$SUBMIT_DIR
submitFile=$runDir/submit

# Make sure submitFile name is unique in case we are continuing a previous
# submission.
if [ -f $submitFile ]; then
  num=1
  while [ -f $submitFile.$num ]; do
    num=$(($num+1))
  done
  submitFile=$submitFile.$num
fi


mkdir -p $runDir

cd $runDir || die "Failed to create directory $runDir"

#
# Job specification
#
echo `pwd`
PrePath=`pwd`
#Executable=`which cmsRun.sh`
Executable=/grp/cms/users/aeverett/farmTools/cmsRun.sh

#
# CMS Dashboard parameters
#
FARMOUT_DASHBOARD_REPORTER=`which farmout_dashboard.sh 2>/dev/null`
if [ "$FARMOUT_DASHBOARD_REPORTER" = "" ]; then
   echo "No farmout_dashboard.sh found, so no reporting to the CMS dashboard."
fi
if [ "$CMS_DASHBOARD_REPORTER" = "" ]; then
   echo "No CMS_DASHBOARD_REPORTER defined, so no reporting to the CMS dashboard."
fi
dboard="
dboard_taskId=${USER}-`hostname -f`-\$(Cluster)
dboard_jobId=\$(Process)
dboard_sid=\$\$([GlobalJobId])
dboard_application=`basename ${CMSSW_HOME}`
dboard_exe=cmsRun
dboard_tool=farmout
dboard_scheduler=local-condor
dboard_taskType=analysis
dboard_broker=local-condor-`hostname -f`
dboard_user=${USER}
dboard_SyncCE=${CMS_DASHBOARD_LOCAL_CE}
CMS_DASHBOARD_REPORTER=${CMS_DASHBOARD_REPORTER}
FARMOUT_DASHBOARD_REPORTER=${FARMOUT_DASHBOARD_REPORTER}
"
# convert newlines to spaces
dboard="`echo $dboard`"

if [ "$SAVE_FAILED_DATAFILES" != "" ]; then
  #cmsRun.sh checks for this in the environment
  save_failed_datafiles_env="SAVE_FAILED_DATAFILES=1"
fi

# First put all the submit file commands that are the same for all jobs.
    cat <<EOF > $submitFile
#Universe             = grid
Universe             = Vanilla
X509UserProxy        = ${proxy}
Executable           = $Executable
GetEnv               = true
gateway		     = osg.rcac.purdue.edu
jobmanager	     =jobmanager-condor
grid_resource        =gt2 \$(gateway)/\$(jobmanager)
globus_rsl           = (condor_submit=(requirements 'regexp(\"cms-*\",Machine)'))
WhenToTransferOutput = On_Exit
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || (ExitCode == ${FAIL_JOB} && JobRunCount>3)))
ImageSize            = $(($MEMORY_REQUIREMENTS*1000))
+DiskUsage           = $(($DISK_REQUIREMENTS*1000))
EOF


echo "Generating submit files in $runDir..."

#
# Loop over input files
#

EXCLUDE_ARG=
if [ "$EXCLUDE_INPUT_FILES" != "" ]; then
  EXCLUDE_ARG="-not -name $EXCLUDE_INPUT_FILES"
fi

#find_command="find $LOCAL_INPUT_DIR -name $MATCH_INPUT_FILES $EXCLUDE_ARG"
find_command="`srmls -2 srm://dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=$LOCAL_INPUT_DIR`|grep root"

dbs_query() {
    $DBSCMD_HOME/dbsCommandLine.py -c lsf "$@" | grep /store/
}
if [ "$INPUT_DBS_PATH" != "" ]; then
  find_command="dbs_query --path=$INPUT_DBS_PATH"
  check_input_file_existence=${ASSUME_INPUT_FILES_EXIST:-1}
  prepend_local_input_dir=1
fi

if [ "$INPUT_FILE_LIST" != "" ]; then
  find_command="cat $INPUT_FILE_LIST"
  check_input_file_existence=${ASSUME_INPUT_FILES_EXIST:-1}
  prepend_local_input_dir=1
fi

count=0
  srmls -2 srm://dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=$LOCAL_INPUT_DIR|grep root |
while read size nextInputFile
do
    inputFileNames=""
    i=$INPUT_FILES_PER_JOB
    while [ $i -gt 0 ]; do
        nextInputFileOrigName="${nextInputFile}"
        if [ "$prepend_local_input_dir" = 1 ]; then
            nextInputFile="${LOCAL_INPUT_DIR}${nextInputFile}"
        fi

        if [ "$check_input_file_existence" = "1" ]; then
            if ! [ -f "$nextInputFile" ]; then
                echo "$nextInputFile does not exist, skipping"

                if [ "$SAVE_MISSING_INPUT_FILE_LIST" != "" ]; then
                    echo "$nextInputFileOrigName" >> "$SAVE_MISSING_INPUT_FILE_LIST"
                fi

		read size nextInputFile || break
                continue
            fi
        fi

        if [ "$inputFileNames" = "" ]; then
            inputFileNames="\"$nextInputFile\""
            firstInputFile="${nextInputFile}"
        else
            inputFileNames="${inputFileNames},\"$nextInputFile\""
        fi

        i=$(($i-1))
        if [ $i -gt 0 ]; then
            read size nextInputFile || break
        fi
    done
    [ "$inputFileNames" = "" ] && break

#
# Name the files
#
    rootname=`basename $firstInputFile .root`
    jobtag=$runName-$rootname
    if [ ${#jobtag} -gt 245 ]; then
        # Condor (as of 6.9.4) cannot handle file names longer than 256
        jobtag=$rootname
    fi

    consub=$jobtag.sub
    conlog=$jobtag.log
    stdout=$jobtag.out
    stderr=$jobtag.err
    jobcfg=$jobtag.cfg
    outputFileName=$jobtag.root

#
# Create and set to the job subdirectory
#

    cd $runDir || die "Failed to cd to $runDir"


    if [ "$SKIP_EXISTING_JOBS" = "1" ] && [ -d $jobtag ]; then
      continue
    fi

    if [ "$SKIP_EXISTING_OUTPUT" = "1" ]; then
      # Check for existing output file
      if outputFileExists $OUTPUT_DIR/$outputFileName; then
        continue
      fi
    fi

    count=$(($count+1))
    if [ ! -z $JOB_LIMIT ] && [ $count -gt $JOB_LIMIT ]; then
        echo "Job limit $JOB_LIMIT reached.  Halting creation of jobs."

        # eat up rest of input to avoid broken pipe message
        while read junk; do junk=""; done

        break
    fi
    echo -n "."

    mkdir -p $jobtag || die "Failed to mkdir $jobtag"

#
# Prepare job configuration file
#

    sed < $configTemplate \
        "s|\\\$inputFileNames|$inputFileNames|g;
         s|\\\$outputFileName|$outputFileName|g" > $jobtag/$jobcfg

#
# Prepare condor submit file for the job
#
    cat >> $submitFile <<EOF

InitialDir           = $PrePath/$jobtag
Arguments            = $jobcfg `basename $outputFileName` $OUTPUT_DIR $SKIP_SRMCP
Transfer_Input_Files = $jobcfg
output               = $stdout
error                = $stderr
Log                  = $conlog
# aaa
# transfer_output_files = $outputFileName
+CMSJob              = True
Queue
EOF
done || exit 1

echo ""


#
# Submit the job
#
if ! grep -q ^Queue $submitFile; then
  echo "No jobs were created, so there is nothing to do."
  #rm $submitFile
  exit 0
elif [ -z "$NO_SUBMIT" ]; then
  # The job is messed up if X509_USER_PROXY is defined, because then
  # Condor doesn't override this value to point to the actual proxy
  # location on the execution node.
  unset X509_USER_PROXY

  condor_submit $submitFile || die "Failed to submit $submitFile"
else
  echo "Submit file $submitFile has been created but not submitted."
fi

echo -n "Jobs for $jobName are created in "
