#!/bin/bash
#  choose one: --> with 2 args, second arg must be CoCoALib dir
#  src/CoCoA-5/release-files/release-linux.sh .
#  release-files/release-linux.sh ../..
#  cd release-files/; release-linux.sh
#--------------------------------------

# WEBSITE="www.dima.unige.it:/Volumes/WWW/webcocoa"
# WEBSITE="130.251.60.18:/Volumes/WWW/webcocoa"
WEBSITE="WWW/webcocoa"

#--------------------------------------
# CVS directory with latest version (either arg or default value)
if [ $# = 0 ]
then
  COCOALIBDIRWithCVS=`cd ../../..; pwd`
else
  COCOALIBDIRWithCVS=`cd "$1"; pwd`
fi

#--------------------------------------
source "$COCOALIBDIRWithCVS/src/CoCoA-5/release-files/release-common.sh"

#--------------------------------------
# names for release directories
COCOA_TEXT=cocoa-$VER

FULLPATH_RELEASE_TEXT_DIR=$RELEASE_DIR/$COCOA_TEXT

#------------------------------------------------------------
cd $RELEASE_DIR/   # <------------ always assume we are here
#------------------------------------------------------------
echo " --======-- CoCoA for linux 32/64 --======--"
echo " --REMEMBER to copy the executables in ~/bin: --vvvvvvvvvvvvvvvvv"
echo "cp $COCOALIBDIRWithCVS/src/CoCoA-5/CoCoAInterpreter ~/bin/CoCoAInterpreter-64"
echo " --REMEMBER! --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

###check-file "$HOME/bin/CoCoAInterpreter-32"
check-file "$HOME/bin/CoCoAInterpreter-64"

cd ..
make texdoc
make htmldoc
cd $RELEASE_DIR/   # <------------ always assume we are here

#------------------------------------------------------------
echo " --vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv--"
echo " GENERATING RELEASE from sources in"
echo "   $COCOALIBDIRWithCVS"
echo " HAVE YOU CHECKED OUT AND COMPILED?  If not, do it with either:"
echo "   cd $COCOALIBDIRWithCVS; ./configure --again; make"
echo "   cd $COCOALIBDIRWithCVS; make"
echo " RELEASE DIR(s) will be:"
echo "   $FULLPATH_RELEASE_TEXT_DIR"
echo " --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--"
#make -j3 library; make -j3

#------------------------------------------------------------
rename-existing-release "$FULLPATH_RELEASE_TEXT_DIR"

#------------------------------------------------------------
copy-packages    $FULLPATH_RELEASE_TEXT_DIR
copy-CoCoAManual $FULLPATH_RELEASE_TEXT_DIR
copy-emacs       $FULLPATH_RELEASE_TEXT_DIR

mkdir -p $FULLPATH_RELEASE_TEXT_DIR/bin
###/bin/cp ~/bin/CoCoAInterpreter-32  $FULLPATH_RELEASE_TEXT_DIR/bin/.
/bin/cp ~/bin/CoCoAInterpreter-64  $FULLPATH_RELEASE_TEXT_DIR/bin/.
/bin/cp cocoa5-linux  $FULLPATH_RELEASE_TEXT_DIR/cocoa5
/bin/cp ConfigEmacs-linux.sh  $FULLPATH_RELEASE_TEXT_DIR/emacs/ConfigEmacs.sh
chmod +x  $FULLPATH_RELEASE_TEXT_DIR/cocoa5
chmod +x  $FULLPATH_RELEASE_TEXT_DIR/emacs/ConfigEmacs.sh
MakeTGZ "$COCOA_TEXT" "linux"
echo "...done"

#------------------------------------------------------------
echo " --======-- suggest-sftp --======--"
echo "sftp storage1.dima.unige.it"
suggest-sftp $FULLPATH_RELEASE_TEXT_DIR "linux" tgz
echo "exit"
echo "touch $WEBSITE/download/download5.shtml"
echo " --======-- end --======--"

echo " --======-- CoCoA for linux 32/64 --======-- --REMEMBER! --"
