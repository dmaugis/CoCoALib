#!/bin/bash
#  cd release-files/; release-win.sh
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

FULLPATH_RELEASE_TEXT_DIR=/Applications/$COCOA_TEXT.$MINORMINOR-win/$COCOA_TEXT

#------------------------------------------------------------
cd $RELEASE_DIR/   # <------------ always assume we are here
#------------------------------------------------------------
echo " --======-- CoCoA for windows --======--"
echo " --REMEMBER to copy the executable in ToBeCopied: --vvvvvvvvvvvvvvvvv"
echo " --REMEMBER! --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

mkdir -p $FULLPATH_RELEASE_TEXT_DIR
check-file "$FULLPATH_RELEASE_TEXT_DIR/../ToBeCopied/CoCoAInterpreter.exe"

cd ..
make texdoc
make htmldoc
cd $RELEASE_DIR/   # <------------ always assume we are here

#------------------------------------------------------------
echo " --vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv--"
echo " GENERATING RELEASE from sources in"
echo "   $COCOALIBDIRWithCVS"
echo " HAVE YOU CHECKED OUT AND COMPILED on WINDOWS?  If not, do it:"
echo " RELEASE DIR(s) will be:"
echo "   $FULLPATH_RELEASE_TEXT_DIR"
echo " --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--"

#------------------------------------------------------------
rename-existing-release "$FULLPATH_RELEASE_TEXT_DIR"

#------------------------------------------------------------
echo " --======-- text CoCoA for Windows  --======--"
copy-packages    $FULLPATH_RELEASE_TEXT_DIR
copy-CoCoAManual $FULLPATH_RELEASE_TEXT_DIR
copy-emacs-el    $FULLPATH_RELEASE_TEXT_DIR

cd $FULLPATH_RELEASE_TEXT_DIR/..

/bin/cp ToBeCopied/CoCoAInterpreter.exe $FULLPATH_RELEASE_TEXT_DIR/
/bin/cp ToBeCopied/cyg* $FULLPATH_RELEASE_TEXT_DIR/
/bin/cp ToBeCopied/*emacs $FULLPATH_RELEASE_TEXT_DIR/emacs

MakeZIP "$COCOA_TEXT" "text-win"
echo "...done"

#------------------------------------------------------------
echo " --======-- suggest-sftp --======--"
echo "sftp storage1.dima.unige.it"
suggest-sftp $FULLPATH_RELEASE_TEXT_DIR "text-win" zip
echo "exit"
echo "touch $WEBSITE/download/download5.shtml"
echo " --======-- end --======--"

