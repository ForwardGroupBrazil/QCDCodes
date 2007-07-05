#!/bin/bash

jobName=$1
blockName=$2
job=$3
COMMON_DIR=$4

echo "makeBatchSubmit $jobName $blockName $job $COMMON_DIR"

submitFile=$COMMON_DIR/${blockName}_${job}

fileName=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d", $1, $2, $3)}'`

echo "fileName $fileName"

# First put all the submit file commands that are the same for all jobs.
cat << EOF > $submitFile
#!/bin/csh -f
# LSF parameters
#
#BSUB -q cmsprs
#BSUB -R "mem > 256"
#BSUB -u aeverett@purdue.edu 
#BSUB -o ${COMMON_DIR}/${blockName}_${job}.log
#BSUB -J ${blockName}_${job}
#
# ---------------------------------------------------------
#
set nonomatch

echo ""
echo "Job is running on \`uname -a\`"
if ( \${OS} == "Linux" ) then
  set processor = \`sort /proc/cpuinfo | uniq | gawk -F: '(substr(\$1,1,10)=="model name"){print \$2}'\`
  set rate = \`sort /proc/cpuinfo | uniq | gawk -F: '(substr(\$1,1,7)=="cpu MHz"){print substr(\$2,1,6)}'\`
  echo "Processor info : " \$processor \$rate "MHz"
endif
echo "Job started on \`date\`"
echo ""

echo " "
echo "Workdir is : " \$WORKDIR
echo " "
echo "Job has been submitted from : " \$LS_SUBCWD

if ( \${OS} != "Solaris" && \${OS} != "Linux" ) exit


echo -----------------------------------------
ls -al
echo -----------------------------------------
echo Copying $COMMON_DIR/${fileName}.cfg  to \$WORKDIR
cp $COMMON_DIR/${fileName}.cfg \$WORKDIR/.
ls -al
echo -----------------------------------------

cd \${CMSSW_BASE}/src
eval \`scramv1 run -csh\`
cd -
setenv CMSSW_SEARCH_PATH \`pwd\`:$CMSSW_SEARCH_PATH
pwd

rm ${fileName}.log
cmsRun -p ${fileName}.cfg >& ${fileName}.log

echo -----------------------------------------
pwd
ls -al
echo -----------------------------------------
    cp \${WORKDIR}/*.root ${COMMON_DIR}/.
    cp \${WORKDIR}/*.log  ${COMMON_DIR}/.  

EOF
