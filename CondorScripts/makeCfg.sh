#!/bin/bash

jobName=$1
configTemplate=$2
configFileBlock=$3
declare -i firstFile=$4
declare -i lastFile=$5
declare -i filesPerJob=$6
COMMON_DIR=$7

blockName=`basename $configFileBlock .plain`

lines=`wc -l $configFileBlock | awk '{printf("%d",$1)}'`

Pt=`basename $configFileBlock .plain | sed -e "s/Mu13/YYY/g" \
    -e "s/Mu-13/ZZZ/g" | cut -c4-`

if [ $firstFile -lt 0 ]; then
    let firstFile=1
    let lastFile=$lines
fi

let nFiles=$lastFile-$firstFile+1

#
# Get an array of file names
#
count=1
while IFS= read -r file
do
  files[count]=$file
  count=$(($count+1))
done < $configFileBlock

while read  -a array;
  do
  count=$(($count+1))
  if [ ${array[0]} == $Pt  ]; then
      let minPt=${array[1]}
      let maxPt=${array[2]}
  fi
done < MuValidSubstitutions.txt


#
# Starting values for the job loop
#
declare -i nFilesSubmitted=0
declare -i totalJobs=1
declare -i fileIndex=$firstFile

#
# Make the CFG files
#
while [ $nFiles -gt $nFilesSubmitted ]; do
    rm tmp.txt
    let indexCount=1
    while [ $indexCount -le $filesPerJob ]; do
	if [ $fileIndex -le $lastFile ]; then     
	    cat << EOFFF >> tmp.txt
            replace PoolSource.fileNames +=  '${files[$fileIndex]}'
EOFFF
	fi
	let indexCount=$indexCount+1
	let fileIndex=$fileIndex+1
    done

    fileName=`echo $jobName $blockName $totalJobs | awk '{printf("%s-%s-%4.4d", $1, $2, $3)}'`

    sed -e '/#inputFileBlock/ r 'tmp.txt'' \
	-e "s/\\\$minPt/$minPt/" \
	-e "s/\\\$maxPt/$maxPt/" \
	-e "s/\\\$outputFileName/$fileName/g" < ${configTemplate} > ${fileName}.cfg

    mv ${fileName}.cfg ${COMMON_DIR}/${fileName}.cfg

#
# Prepare for the next job
#
    let nFilesSubmitted=$nFilesSubmitted+$filesPerJob
    if (( $nFilesSubmitted < $nFiles )); then
	let totalJobs=$totalJobs+1
#    else
#	echo "$nJobsSubmitted >= $totalJobs done!"
    fi
done

rm tmp.txt

echo $totalJobs



