#!/bin/bash

castorDir=$1
castorDir2=$2

echo ${castorDir}
echo ${castorDir2}

read -r -a FileList <<< $(srmls -2 srm://dcache.rcac.purdue.edu:8443/srm/manager
v2?SFN=${castorDir}/${castorDir2})
#read -r -a FileList <<< $(ls ${castorDir}/${castorDir2})

element_count=${#FileList[@]}
let "end_count = $element_count - 1"
echo "element_count $element_count"

echo "${FileList[0]} ${FileList[1]} ${FileList[2]} ${FileList[3]}"

index=3
#index=0

match="root"
while [ "$index" -lt "$element_count" ]
  do    # List all the elements in the array.
  if [[ "${FileList[$index]}" =~ "${match}" ]]; then
      cat << EOF >> ${castorDir2}_manip.txt
   srmrm -2 srm://dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=${FileList[$inde
x]}
EOF

  fi
  let "index = $index + 2"
#  let "index = $index + 1"
done

