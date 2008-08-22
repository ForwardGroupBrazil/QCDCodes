#!/bin/sh
#
RUNTIME_AREA=`pwd`
jobcfg=$1
datafile=$2
SRM_OUTPUT_DIR=$3
SKIP_SRMCP=$4
SRM_OUTPUT_FILE="$SRM_OUTPUT_DIR/$datafile"
SRM_FAILED_OUTPUT_FILE="${SRM_OUTPUT_DIR}-cmsRun-failed/${datafile}"
OSG_SETUP=/opt/osg/setup.sh

# special exit status to force job to leave the queue
FAIL_JOB=42

# core files have a nasty habit of filling up disks and causing trouble,
# so disable them.
ulimit -c 0

# load environment for using srmcp
source $OSG_SETUP
#WORKING_DIR=`/bin/mktemp -d /tmp/cms_XXXXXXXXXXXX`
#   if [ ! $? == 0 ] ; then
#       echo "$WORKING_DIR could not be created on WN `hostname`"
#       exit 1
#   fi
#   echo ">>> Created working directory: $WORKING_DIR"
#   cd $WORKING_DIR
#   echo ">>> current directory (WORKING_DIR): $WORKING_DIR"
#
#   export SCRAM_ARCH=slc4_ia32_gcc345
#   source /apps/02/cmssoft/cms/cmsset_default.sh
#   #source /grp/cms/purdueCMS/cmsset_default.sh
#   scramv1 project CMSSW CMSSW_2_0_11
#status=$?
#if [ $status != 0 ] ; then
#    echo "ERROR CMSSW"
#    rm -rf $WORKING_DIR
#    exit 1
#fi
#echo ">> create CMSSW_2_0_11"
#cd CMSSW_2_0_11/src
#SOFTWARE_DIR=`pwd`
#echo $SOFTWARE_DIR
#eval `scramv1 runtime -sh`
#which cmsRun
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/src/Configuration .
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/src/Validation .
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/lib ../lib/.
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/bin ../bin/.
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/share ../share/.
#cp -r /usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/module ../module/.
##export PATH=/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/bin/slc4_ia32_gcc345:$PATH
##export CMSSW_SEARCH_PATH=/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/src:/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/share:$CMSSW_SEARCH_PATH
##export LD_LIBRARY_PATH=/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/lib/slc4_ia32_gcc345:/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/external/slc4_ia32_gcc345:$LD_LIBRARY_PATH
##export PYTHONPATH=/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/lib/slc4_ia32_gcc345:$PYTHONPATH
##export SEAL_PLUGINS=/usr/rmt_share/scratch96/x/xu2/CMSSW_2_0_11/module/slc4_ia32_gcc345:$SEAL_PLUGINS
#hostname
dashboard_completion() {
  export dboard_ExeExitCode=$1

  if [ "${FARMOUT_DASHBOARD_REPORTER}" != "" ]; then
      ${FARMOUT_DASHBOARD_REPORTER} completion
  fi
}

outputFileExists() {
  srm_fname="$1"

  info="`srm-get-metadata -retry_num=0 $srm_fname 2>/dev/null`"
  if [ "$?" != "0" ]; then
    return 1
  fi
  size_line="`echo \"$info\" | grep \"size :\"`"
  if [ "$?" != "0" ]; then
    return 1
  fi
  IFS=$' :\t'
  size="`echo \"$size_line\" | ( read label size; echo $size )`"
  unset IFS
  if [ "$size" != "0" ] && [ "$size" != "" ]; then
    return 0
  fi
  return 1
}

DoSrmcp() {
    src="$1"
    dest="$2"

    if [ "$SKIP_SRMCP" = "1" ]; then
	echo "SKIP DoSrmcp ------------"
	return 0
    fi
    echo "Not SKIPing DoSrmcp ------------"

    # The internal retries in srmcp are often useless, because a broken file
    # may be left behind (as of dCache 1.7.0-38), so further attempts to
    # copy the file will all fail with the error "file exists".  Therefore,
    # we have our own retry loop, including deletion of broken files.

    tries=0
    while [ "$tries" -lt 3 ]; do
      if [ "$tries" -gt 1 ]; then
        echo "Trying again at `date`: srmcp $src $dest"
      fi

      srmcp -2 -debug=true -retry_num=0 "$src" "$dest"
      rc=$?

      if [ "$rc" = "0" ]; then
        return 0
      fi

      echo
      echo
      echo "srmcp exited with non-zero status $rc at `date`."
      echo "This happened when copying $src to $dest."
      echo

      if outputFileExists "$dest"; then
         echo "Cleaning up failed destination file $dest."
         srm-advisory-delete -debug=true -retry_num=0 "$dest"

         rc=$?
         if [ "$rc" != "0" ]; then
           echo "srm-advisory-delete failed with exit status $rc at `date`."
         fi
      fi

      tries=$(($tries+1))
    done

    echo "Giving up after $tries attempts to copy $src to $dest."
    return 1
}


if outputFileExists $SRM_OUTPUT_FILE; then
  echo "File already exists: $SRM_OUTPUT_FILE; exiting as though successful."
  exit 0
fi

if [ "${FARMOUT_DASHBOARD_REPORTER}" != "" ]; then
    ${FARMOUT_DASHBOARD_REPORTER} submission
    ${FARMOUT_DASHBOARD_REPORTER} execution
fi

start_time=`date "+%s"`
echo $start_time
cp $RUNTIME_AREA/$jobcfg .
cmsRun $jobcfg
cmsRun_rc=$?

export dboard_ExeTime=$((`date "+%s"` -  $start_time))

echo "ls -ltr"
ls -ltr
echo "End of ls output"

if [ "$cmsRun_rc" != "0" ]; then
  echo "cmsRun exited with status $cmsRun_rc"
# aaa section 1
  if [ -f $datafile ] && [ "$SAVE_FAILED_DATAFILES" = "1" ]; then
      echo "----- Section 1 -----"
      if ! DoSrmcp "file:///`pwd`/$datafile" "$SRM_FAILED_OUTPUT_FILE"; then
	  echo "Failed to save datafile from failed run."
      fi
  fi
# aaa end
  rm -f $datafile

  dashboard_completion $cmsRun_rc

  # Do not try to run this job again.  Would be nice to detect transient
  # errors (e.g. dCache down) and retry in those cases.
  exit $FAIL_JOB
fi

if ! [ -f $datafile ]; then
  echo "cmsRun did not produce expected datafile $datafile"
  exit $FAIL_JOB
fi

# aaa section 2
if ! DoSrmcp "file:///`pwd`/$datafile" "$SRM_OUTPUT_FILE"; then
    echo "----- Section 2 -----"
    dashboard_completion 60307
    exit 1
fi
# aaa end

# aaa section 3
if [ "$SKIP_SRMCP" = "0" ]; then
   rm $datafile
fi
##rm -rf $WORKING_DIR
# Copy all other root files in the directory also

# aaa section 4
if [ "$SKIP_SRMCP" = "0" ]; then
    for file in `ls -1 *.root 2>/dev/null`; do
	if ! DoSrmcp file:///`pwd`/$file $SRM_OUTPUT_DIR/$file; then
	    dashboard_completion 60307
	    exit 1
	fi
 # aaa section 5
	rm $file
    done
fi
# aaa end

dashboard_completion 0

exit 0
