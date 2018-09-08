#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# Auxiliary script for CoCoALib configuration process.
# Expects env variable CXX to be set (to compiler's name).

# Script to see whether compiler is clang, and then link with special flags.
# If no warning is produced, the script prints -fPIC; otherwise it prints nothing.

if [ $# -ne 0 ]
then
  echo "ERROR: expected no args   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Check environment variable CXX
if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Create tmp directory, put test prog in it, compile and run.
umask 22
TODAY=`date "+%Y-%m-%d"`
TIME=`date "+%H:%M:%S"`
TMP_DIR=/tmp/CoCoALib-config-$USER-$TODAY/fpic-ldflag-$TIME-$$
/bin/rm -rf $TMP_DIR  &&  /bin/mkdir -p $TMP_DIR
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

cd $TMP_DIR

# test if it is clang:  .... a bit harsh, maybe...
/bin/cat > TestProg.C <<EOF
int main()
{
#ifdef __clang__
  exit(1);
#endif
}
EOF

FPIC_FLAG=-fPIC

"$CXX" -o TestProg TestProg.C > LogFile  2>& 1  &&  ./TestProg >> LogFile  2>&1
if [ $? -ne 0 ]
then
  FPIC_LDFLAG="-Wl,-no_pie";
fi

# Clean up TMP_DIR
cd
/bin/rm -rf "$TMP_DIR"
echo $FPIC_LDFLAG
