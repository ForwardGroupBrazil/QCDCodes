#!/bin/bash

castorDir=$1
castorDir2=$2

echo ${castorDir}
echo ${castorDir2}

read -r -a FileList <<< $(rfdir ${castorDir}/${castorDir2})

element_count=${#FileList[@]}
let "end_count = $element_count - 1"

index=8

cat << EOF >> ${castorDir2}.txt
replace PoolSource.fileNames = {
EOF
match="edm"
first=1
while [ "$index" -lt "$element_count" ]
  do    # List all the elements in the array.
#  if [ "$index" -ne "$end_count" ]; then
      if [[ "${FileList[$index]}" =~ "${match}" ]]; then
	  if [ $first -eq 1 ]; then
	      let "first = 0"
	      cat << EOF >> ${castorDir2}.txt
   'rfio:${castorDir}/${castorDir2}/${FileList[$index]}'
EOF
	  fi
	  cat << EOF >> ${castorDir2}.txt
   ,'rfio:${castorDir}/${castorDir2}/${FileList[$index]}'
EOF
#      fi
  fi
  let "index = $index + 9"
done

cat << EOF >> ${castorDir2}.txt
}
EOF
