#!/bin/csh -f
set nonomatch

set CONDOR_CLUSTER = $1
set CONDOR_PROCESS = $2
set PARAMETER_SET  = $3
set COMMON_DIR     = $4

echo "cluster " ${CONDOR_CLUSTER}
echo "process " ${CONDOR_PROCESS}

echo ""
echo "Job is running on `uname -a`"
echo ""
if ( ${OSTYPE} == "linux" ) then
  set processor = `sort /proc/cpuinfo | uniq | gawk -F: '(substr($1,1,10)=="model name"){print $2}'`
  set rate = `sort /proc/cpuinfo | uniq | gawk -F: '(substr($1,1,7)=="cpu MHz"){print substr($2,1,6)}'`
  echo "Processor info : " $processor $rate "MHz"
endif
set start = `date`
echo "Job started on `date`"
echo ""

pwd
ls

setenv WORKDIR `pwd`

echo "cluster " ${CONDOR_CLUSTER}
echo "process " ${CONDOR_PROCESS}

#
#----------------------------------------------------------
# s e t   t h e   r u n t i m e   e n v i r o n m e n t
#----------------------------------------------------------
#
#-- Is this necessary if GetEnv=true ?
setenv PATH /usr/local/bin:/usr/bin:/bin:/opt/osg/srmclient/bin
setenv CMS_PATH /apps/02/cmssoft/cms/
source ${CMS_PATH}/cmsset_default.csh

setenv CMSSW_RELEASE CMSSW_2_0_11
setenv SCRAM_ARCH slc4_ia32_gcc345 

setenv SRCDIR /scratch/scratch96/a/aeverett/${CMSSW_RELEASE}

echo " "
echo "Output directory is : " $COMMON_DIR
echo " "

#
#----------------------------------------------------------
# s e t   t h e   r u n t i m e   e n v i r o n m e n t
#----------------------------------------------------------
#

cd ${SRCDIR}/src
#- Is this necessary of GetEnv=true ?
eval `scramv1 runtime -csh` 

#
#----------------------------------------------------------
# c o p y   e x e   a n d   c o n f i g    f i l e s
#----------------------------------------------------------
#
cd ${WORKDIR}
echo "Workdir " ${WORKDIR}
 

#
#----------------------------------------------------------
# e x e c u t e   j o b
#----------------------------------------------------------
#

 echo "Processing Event Collection :  " 

 /usr/bin/time -f "%e %U %S %x" -o x cmsRun -p ${PARAMETER_SET}

 set rtime = `tail -1 x | cut -f 1 -d " "`
 set utime = `tail -1 x | cut -f 2 -d " "`
 set stime = `tail -1 x | cut -f 3 -d " "`
 set stat  = `tail -1 x | cut -f 4 -d " "`

#
#----------------------------------------------------------
# c o p y   o u t p u t
#----------------------------------------------------------
#
#  if ( \${status} != 0 ) cp \${WORKDIR}/*.root \${COMMON_DIR}/H_1.root
#    cp \${WORKDIR}/*.root \${COMMON_DIR}/${basename}.root
    cp ${WORKDIR}/*.root ${COMMON_DIR}/.
    cp ${WORKDIR}/*.log  ${COMMON_DIR}/.  
#
# show what is being left behind...
#
  echo ""
  echo "Current directory:"
  pwd
  ls -lrtAFh
#
  set end = `date`
  echo ""
  echo "Job end `date`"
  echo ""
#
#
exit ${status}
    
    
    
