#!/bin/bash

castorDir=$1
castorDir2=$2

echo ${castorDir}
echo ${castorDir2}

#read -r -a FileList <<< $(rfdir /castor/cern.ch/user/a/aeverett/note-reco209-01/minusPt1000)
read -r -a FileList <<< $(rfdir ${castorDir}/${castorDir2})

element_count=${#FileList[@]}
let "end_count = $element_count - 1"

index=8

#   for i in "${FileList[@]}"
#   do
#     echo "$i"
#   done

cat << EOF >> ${castorDir2}.txt
replace PoolSource.fileNames = {
EOF

while [ "$index" -lt "$element_count" ]
  do    # List all the elements in the array.
  if [ "$index" -ne "$end_count" ]; then
      cat << EOF >> ${castorDir2}.txt
   'rfio:${castorDir}/${FileList[$index]}',
EOF
  fi
  if [ "$index" -eq "$end_count" ]; then
      cat << EOF >> ${castorDir2}.txt
  'rfio:${castorDir}/${FileList[$index]}'
EOF
  fi
  let "index = $index + 9"
  # Or:
  #    index+=1
  # if running Bash, version 3.1 or later.
done

cat << EOF >> ${castorDir2}.txt
}
EOF

#echo ${FileList[17]}


