#!/bin/bash

jobName=$1
blockName=$2
totalJobs=$3
COMMON_DIR=$4

#
# Name the condor submit script
#
submitFile=$COMMON_DIR/submit

# Make sure submitFile name is unique in case we are continuing a previous
# submission.
if [ -f $submitFile ]; then
    num=1
    while [ -f $submitFile.$num ]; do
	num=$(($num+1))
    done
    submitFile=$submitFile.$num
fi

# First put all the submit file commands that are the same for all jobs.
cat << EOF > $submitFile
Universe             = vanilla
Executable           = condorRunScript.csh
Copy_To_Spool        = false
Notification         = never
GetEnv               = true
#should_transfer_files = YES
WhenToTransferOutput = On_Exit
on_exit_remove       = (ExitBySignal == FALSE && ExitStatus == 0)
+IsFastQueueJob      = True
requirements = (Arch == "X86_64")&&regexp("cms",Name)
+CMSJob = True
EOF

declare -i job=1
while [ $job -le $totalJobs ]; do
    echo "adding to subscript for $job of $totalJobs"
#
# Name the files
#    
    conlog=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d.log", $1, $2,$3)}'`
    stdout=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d.out", $1, $2,$3)}'`
    stderr=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d.err", $1, $2,$3)}'`
    
    jobcfg=`echo $jobName $blockName $job | awk '{printf("%s-%s-%4.4d.cfg", $1, $2,$3)}'`
    
    let job=$job+1    
    
#
# Prepare condor submit file for the job
#
    cat >> $submitFile <<EOF

InitialDir           = $COMMON_DIR
Arguments            = \$(CLUSTER) \$(PROCESS) $jobcfg $COMMON_DIR
Transfer_Input_Files = $jobcfg
output               = $stdout
error                = $stderr
log                  = $conlog

Queue
EOF
    
done
