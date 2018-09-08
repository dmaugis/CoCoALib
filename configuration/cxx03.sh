#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# Auxiliary script for CoCoALib configuration process.
# Script expects the env variables CXX and CXXFLAGS to be set.

# Script to see whether the -std=c++03 compiler flag is recognised.
# Output is either "-std=c++03" (if the flag is recognized, otherwise the empty string.


if [ $# -ne 0 ]
then
  echo "ERROR: expected no args.   $SCRIPT_NAME"  > /dev/stderr
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
TMP_DIR=/tmp/CoCoALib-config-$USER-$TODAY/cxx03-$TIME-$$
/bin/rm -rf $TMP_DIR  &&  /bin/mkdir -p $TMP_DIR
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\".   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

cd $TMP_DIR

/bin/cat > test.C <<EOF
int main()
{}
EOF

CXX03="-std=c++03"
"$CXX" $CXX03 test.C  >> LogFile  2>& 1 
if [ $? -ne 0 ]
then
  CXX03=
fi

# Clean up TMP_DIR
cd
/bin/rm -rf "$TMP_DIR"
echo $CXX03
