#!/bin/bash
#  cd release-files/; release-win-fromMac.sh
#--------------------------------------

# WEBSITE="www.dima.unige.it:/Volumes/WWW/webcocoa"
# WEBSITE="130.251.60.18:/Volumes/WWW/webcocoa"
WEBSITE="WWW/webcocoa"

#--------------------------------------
# CVS directory with latest version
COCOALIBDIRWithCVS=`cd ../../..; pwd`

#--------------------------------------
source "$COCOALIBDIRWithCVS/src/CoCoA-5/release-files/release-common.sh"

source "$COCOALIBDIRWithCVS/configuration/version"
COCOALIB=CoCoALib-$VER_MAJ.$VER_MIN$VER_PATCH
#--------------------------------------
# names for release directories
COCOA_TEXT=cocoa-$VER

FULLPATH_RELEASE_TEXT_DIR=/Applications/$COCOA_TEXT.$MINORMINOR-win/$COCOA_TEXT

BINARY_DIR=$COCOALIBDIRWithCVS/../../CoCoABinaries/Windows-$VER.$MINORMINOR

#------------------------------------------------------------
cd $RELEASE_DIR/   # <------------ always assume we are here
#------------------------------------------------------------
echo " --======================--CoCoA-for-WINDOWS--======================--"
echo " --send-$COCOALIB.tgz:"
echo "cp ~/tmp/$COCOALIB.tgz /Volumes/SharedFolder/"
echo " --compile: ~/shell-scripts/cocoa5-win"
echo " --check emacs files in $FULLPATH_RELEASE_TEXT_DIR/../ToBeCopied"
echo " --copy-CoCoAInterpreter.exe:"
echo "cp /Volumes/SharedFolder/CoCoAInterpreter.exe $BINARY_DIR/"
echo " --REMEMBER! --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

check-file "$BINARY_DIR/CoCoAInterpreter.exe"

cd ..
make doc
cd $RELEASE_DIR/   # <------------ always assume we are here

#------------------------------------------------------------
echo " --generating_release from $COCOALIBDIRWithCVS"
echo " --release_dir: $FULLPATH_RELEASE_TEXT_DIR"

#------------------------------------------------------------
rename-existing-release "$FULLPATH_RELEASE_TEXT_DIR"

#------------------------------------------------------------
echo " --======-- text CoCoA for Windows  --======--"
copy-packages    $FULLPATH_RELEASE_TEXT_DIR
copy-CoCoAManual $FULLPATH_RELEASE_TEXT_DIR
copy-emacs-el    $FULLPATH_RELEASE_TEXT_DIR

cd $FULLPATH_RELEASE_TEXT_DIR/..

/bin/cp $BINARY_DIR/CoCoAInterpreter.exe $FULLPATH_RELEASE_TEXT_DIR/
/bin/cp ToBeCopied/cyg* $FULLPATH_RELEASE_TEXT_DIR/
/bin/cp ToBeCopied/*emacs $FULLPATH_RELEASE_TEXT_DIR/emacs

MakeZIP "$COCOA_TEXT" "text-win"
echo "...done"

#------------------------------------------------------------
echo " --======-- suggest-sftp --======--"
echo "sftp storage1.dima.unige.it"
suggest-sftp $FULLPATH_RELEASE_TEXT_DIR "text-win" "zip"
echo "exit"
echo "touch /Volumes/$WEBSITE/download/download5.shtml"
echo " --======-- end --======--"
