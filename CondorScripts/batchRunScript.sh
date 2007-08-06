#!/bin/csh -f

if ( $subtype != "PBS" ) then
    set fileName=$1
    set COMMON_DIR=$2
endif

# LSF parameters
#
#$subtype -q cmsprs
#$subtype -R "mem > 256"
#$subtype -oo ${COMMON_DIR}/${fileName}.out
#$subtype -eo ${COMMON_DIR}/${fileName}.err
#
# ---------------------------------------------------------
#

set nonomatch

pwd
echo "fileName: $fileName"
if ( $subtype == "PBS" ) then
    echo ------------------------------------------------------
    echo "PBS: qsub is running on $PBS_O_HOST"
    echo "PBS: originating queue is $PBS_O_QUEUE"
    echo "PBS: executing queue is $PBS_QUEUE"
    echo "PBS: working directory is $PBS_O_WORKDIR"
    echo "PBS: execution mode is $PBS_ENVIRONMENT"
    echo "PBS: job identifier is $PBS_JOBID"
    echo "PBS: job name is $PBS_JOBNAME"
    echo "PBS: current home directory is $PBS_O_HOME"
    #echo "PBS: PATH = \$PBS_O_PATH"
    echo ------------------------------------------------------
endif

echo ""
echo "Job is running on `uname -a`"
if ( `uname -s` == "Linux" ) then
  set processor = `sort /proc/cpuinfo | uniq | gawk -F: '(substr($1,1,10)=="model name"){print $2}'`
  set rate = `sort /proc/cpuinfo | uniq | gawk -F: '(substr($1,1,7)=="cpu MHz"){print substr($2,1,6)}'`
  echo "Processor info : " $processor $rate "MHz"
endif
echo "Job started on `date`"
echo ""

echo " "
if($subtype == "PBS") then
    echo "Running PBS"
    set WORKDIR=/tmp/PBS_${PBS_JOBID}
    mkdir $WORKDIR
endif
echo "Workdir is : " $WORKDIR
echo " "
#echo "Job has been submitted from : " $LS_SUBCWD

if ( `uname -s` != "Solaris" && `uname -s` != "Linux") exit

echo -----------------------------------------
cd $WORKDIR
ls -al
echo -----------------------------------------
echo Copying $COMMON_DIR/${fileName}.cfg  to $WORKDIR
cp $COMMON_DIR/${fileName}.cfg $WORKDIR/.
ls -al
echo -----------------------------------------

cd ${CMSSW_BASE}/src
eval `scramv1 run -csh`
cd $WORKDIR
setenv CMSSW_SEARCH_PATH `pwd`:$CMSSW_SEARCH_PATH
pwd

rm ${fileName}.log
cmsRun -p ${fileName}.cfg >& ${fileName}.log

echo -----------------------------------------
pwd
ls -al
echo -----------------------------------------
cp ${WORKDIR}/*.root ${COMMON_DIR}/.
#cp ${WORKDIR}/*.log  ${COMMON_DIR}/.  

