#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# Auxiliary script for CoCoALib configuration process.
# Script expects the env variables CXX and CXXFAGS to be set.

# Script to see whether the -fPIC flag produces annoying compiler warnings.
# If no warning is produced, the script prints "-fPIC"; otherwise it prints nothing.

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


FPIC_FLAG=-fPIC

# Create tmp directory, put test prog in it, compile and run.
umask 22
TODAY=`date "+%Y-%m-%d"`
TIME=`date "+%H:%M:%S"`
TMP_DIR=/tmp/CoCoALib-config-$USER-$TODAY/fpic-flag-$TIME-$$
/bin/rm -rf $TMP_DIR  &&  /bin/mkdir -p $TMP_DIR
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

cd $TMP_DIR

/bin/cat > test.C <<EOF
int f(int x)
{
  return (x+1)*x+41;
}
EOF


COMPILER_MESG=`"$CXX" $FPIC_FLAG -c -o test.o test.C 2>& 1`
if [ $? -ne 0 ]
then
  echo "ERROR: test compilation failed   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

# Clean up TMP_DIR
cd
/bin/rm -rf "$TMP_DIR"
if [ -z "$COMPILER_MESG" ]
then
  echo $FPIC_FLAG
fi
