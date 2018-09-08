#! /bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# This script expects the env variables CXX and CXXFLAGS to be set.

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


# This script tries to check that CXX is a working C++ compiler,
# and that CXXFLAGS are suitable flags for it.  The check entails
# compiling a very simple source file (with and without CXXFLAGS).

# Check that CXX is an executable file:
FULLCXX=`which "$CXX" 2>/dev/null`
if [ $? -ne 0 -o \! \( -x "$FULLCXX" -a -r "$FULLCXX" -a -f "$FULLCXX" \) ]
then
  echo "Specified compiler \"$CXX\" is not an executable!"  > /dev/stderr
  exit 1
fi


# Create tmp directory, put test prog in it, compile and run.
umask 22
TODAY=`date "+%Y-%m-%d"`
TIME=`date "+%H:%M:%S"`
TMP_DIR=/tmp/CoCoALib-config-$USER-$TODAY/verify-compiler-$TIME-$$
/bin/rm -rf $TMP_DIR  &&  /bin/mkdir -p $TMP_DIR
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

cd $TMP_DIR

# Here is the simple source code we shall use to test the compiler:
/bin/cat > TestProg.C <<EOF
#include <iostream>
using namespace std;
int main()
{
#ifdef __GNUC__
  cout << "gnu";
#else
  cout << "not gnu";
#endif
}
EOF

# Try plain compiler (without CXXFLAGS):
$CXX TestProg.C -o TestProg  > LogFile  2>&1
if [ $? -ne 0 -o \! -f TestProg -o \! -x TestProg ]
then
  echo "ERROR: Are you sure \"$CXX\" is a C++ compiler?   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi
/bin/rm TestProg  # not necessary, just being tidy :-)

# Try compiler with CXXFLAGS:
$CXX $CXXFLAGS TestProg.C -o TestProg  > LogFile  2>&1
if [ $? -ne 0 -o \! -f TestProg -o \! -x TestProg ]
then
  echo "ERROR: Compilation flags \"$CXXFLAGS\" seem to be unsuitable for \"$CXX\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

COMPILER_TYPE=`./TestProg`

# Clean up TMP_DIR
cd
/bin/rm -rf $TMP_DIR

echo "$COMPILER_TYPE"
