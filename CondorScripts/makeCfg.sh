#!/bin/bash

fileOfSamples=MW.txt
COMMON_DIR=${HOME}/data/SLC4/CMSSW_1_7_1/src/UserCode/CondorScripts/tmpdir
subtype=CONDOR
jobName=name1
configTemplate="cfgTemplate.cfg"
recoFragment="offlineReco.cff"
analyzerFragment="MuValidTemplate.cfg"

# setup runtime environment
ORIGINAL_DIR=`pwd`
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd -
if [ "$?" != "0" ]; then
  echo "Failed to initialize CMSSW environment with scram in $CMSSW_BASE."
  exit 1
fi

# make the common directory for this job
mkdir -p ${COMMON_DIR}
echo "COMMON_DIR = ${COMMON_DIR}"

# read in the fileOfSamples and parameters
count=0
while read -a file;
  do
  count=$(($count+1))
  files[count]=${file[0]}
  tag[count]=${file[1]}
  nEvents[count]=${file[2]}
  nEventsJob[count]=${file[3]}
done < ${fileOfSamples}

# loop over all samples
fileName=""
submitFile=""
sample=0
while [ $sample -lt $count ]; do
    sample=$(($sample+1))

    # make a directory for each sample
    RUN_DIR=${COMMON_DIR}/${tag[sample]}
    mkdir -p ${RUN_DIR}
    
    declare -i M=0
    declare -i N=0
    
    if [ "$subtype" = "CRAB" ]; then
	# set the number of subSamples for CRAB to 1
	let M=1
	let N=1
    fi
    if [ "$subtype" = "CONDOR" ]; then
	# set the number of subSamples for CONDOR to depend on the
        # numberOfEventsPerJob and the totalNumberOfEvents
	let M=nEvents[sample]
	let N=nEventsJob[sample]

	# each subSample must have a submit file
	submitFile="${RUN_DIR}/submit"
	if [ -f $submitFile ]; then
	    num=1
	    while [ -f $submitFile.$num ]; do
		num=$(($num+1))
	    done
	    submitFile=$submitFile.$num
	fi	
    fi

    # loop over all subSamples
    skipEvents=0
    subSample=1
    while [ $skipEvents -lt $M ]; do
	
	if [ "$subtype" = "CRAB" ]; then
	    fileName=`echo $jobName ${tag[sample]} | awk '{printf("%s-%s",$1, $2)}'`
	    
	    # customize the crab configuration file
	    rm -f tmpCrab.txt
	    cat << EOF >> tmpCrab.txt
datasetpath = ${files[$sample]}
total_number_of_events = ${nEvents[$sample]}
number_of_jobs = ${nEventsJob[$sample]}
EOF
	    
	    sed -e '/#CMSSWBlock/ r 'tmpCrab.txt'' \
		-e "s/\\\$configFile/${fileName}.cfg/" \
		-e "s/\\\$outFileName/${fileName}/g" < crabTemplate.cfg > ${RUN_DIR}/crab.${fileName}.cfg
	    
	    rm -f tmpCrab.txt
	    
	fi
	if [ "$subtype" = "CONDOR" ];then
	    fileName=`echo $jobName ${tag[sample]} $subSample | awk '{printf("%s-%s-%3.3d",$1, $2, $3)}'`

	    # add common info to the condor submit file
            # First put all the submit file commands that 
            # are the same for all jobs.
	    if [ $subSample -eq 1 ]; then

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
	    fi
	    
	    # name the files for this subSample
	    conlog=${fileName}.log
	    stdout=${fileName}.out
	    stderr=${fileName}.err
	    jobcfg=${fileName}.cfg
	    
	    # prepare condor submit file for the job
	    cat >> $submitFile <<EOF

InitialDir           = $RUN_DIR
Arguments            = \$(CLUSTER) \$(PROCESS) $jobcfg $RUN_DIR
Transfer_Input_Files = $jobcfg
output               = $stdout
error                = $stderr
log                  = $conlog

Queue
EOF
	    
	    cp condorRunScript.csh ${RUN_DIR}/condorRunScript.csh
	    chmod +x ${RUN_DIR}/condorRunScript.csh	    

	fi

	# prepare the configuration file
	rm -f tmp.cfg
	sed -e '/#inputFileBlock/ r '${files[sample]}'' \
	    -e '/#recoBlock/ r '${recoFragment}'' \
	    -e '/#analyzerBlock/ r '${analyzerFragment}'' < ${configTemplate} > tmp.cfg
	
	sed -e "s/\\\$skipEvents/${skipEvents}/" \
	    -e "s/\\\$outFileName/${fileName}/" \
	    -e "s/\\\$nEvents/${N}/" < tmp.cfg > ${RUN_DIR}/${fileName}.cfg
	
	subSample=$(($subSample+1))
	skipEvents=$(($skipEvents+N))
	
	rm -f tmp.cfg
	
    done

    if [ "$subtype" = "CRAB" ]; then
	cd ${RUN_DIR}
#	crab -cfg crab.${fileName}.cfg -create -submit
	cd -
	echo "crab -cfg crab.${fileName}.cfg -create -submit"
    fi

    if [ "$subtype" = "CONDOR" ]; then
#	condor_submit $submitFile
	echo "condor_submit $submitFile"
    fi
    
done

