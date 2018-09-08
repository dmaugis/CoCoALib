#! /bin/sh

SCRIPT_NAME=[[`basename "$0"`]]

# Script assumes that CXX and CXXFLAGS are set.
# This script tests whether CXX can compile with readline by
# just adding -lreadline as flag.
# The check entails compiling a very simple source file.

# Exit code is 0 if "-lreadline" works.
# Exit code is 1 if "-lreadline" does not work.
# Exit code is 2 if there was bad input or could not create tmp dir.
# (for debugging, an error message is printed if exit code is 2)

if [ $# -ne 0 ]
then
  echo "ERROR: expected no args.   $SCRIPT_NAME"  > /dev/stderr
  exit 2
fi
   

if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 2
fi


# Create tmp directory, put test prog in it, compile and run.
umask 22
TODAY=`date "+%Y-%m-%d"`
TIME=`date "+%H:%M:%S"`
TMP_DIR=/tmp/CoCoALib-config-$USER-$TODAY/readline-try-default-$TIME-$$
/bin/rm -rf $TMP_DIR  &&  /bin/mkdir -p $TMP_DIR
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"  >/dev/stderr
  exit 2
fi

cd $TMP_DIR

# Here is the simple source code we shall use to test for readline:
/bin/cat > test-readline.C <<EOF
#include "stdlib.h"
#include "stdio.h"
#include "unistd.h"
#include "readline/readline.h"
#include "readline/history.h"

#include <string>

int main()
{
  char *line;
  line = readline("");
  std::string str(line);
  free(line);
  // We shall test on an input of 3 chars...
  if (str.size() != 3) return 1;
  return 0;
}
EOF

"$CXX" $CXXFLAGS test-readline.C -lreadline -o test-readline  > LogFile 2>&1
if [ $? -ne 0 ]
then
    exit 1
fi

# Successful, so clean up TMP_DIR
cd
/bin/rm -rf $TMP_DIR
exit 0

