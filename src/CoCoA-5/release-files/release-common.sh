#----------------------------------------------
source version
VER=$COCOA5_VER_MAJ.$COCOA5_VER_MIN
MINORMINOR=$COCOA5_VER_MINMIN
#----------------------------------------------
REL=$VER.$MINORMINOR
#----------------------------------------------
DATE=`date "+%Y%m%d-%H%M"`
YEAR=`date "+%Y"`
#----------------------------------------------
# Convert the file name into an absolute path (in case it was not)
THIS_DIR=`dirname "$$0"`
RELEASE_DIR=`cd "$THIS_DIR"; pwd`

#----------------------------------------------

check-file()
{
  if [ \! -x "$1" ]
  then
      echo "$0: ERROR cannot find $1"
      exit 1
  fi
  if [ "$RELEASE_DIR/release-common.sh" -nt "$1" ]
  then
      echo "$0: ERROR file too old $1"
      echo " -----> if you know what you are doing...  ---->"
      echo " touch $1"
      exit 1
  fi
}

rename-existing-release()
{
  if [ -d "$@" ]
  then
      echo "  --renaming --> `basename $@`-$DATE"
      /bin/mv "$@" "$@-$DATE"
      echo "  rm -rf $@-$YEAR*"
  fi
}

#----------------------------------------------

copy-packages()
{
  echo "Copying packages in    ..${1:(-40)}/";  # last chars
  # rm -rf $1/packages;
  mkdir -p $1/packages;
  /bin/cp -R  ../packages/*  $1/packages/;
  /bin/rm -rf $1/packages/CVS/;
}

copy-CoCoAManual()
{
  echo "Copying CoCoAManual in ..${1:(-40)}/";  # last chars
  # rm -rf $1/CoCoAManual/;
  mkdir -p $1/CoCoAManual;
  /bin/cp -R  ../CoCoAManual/*  $1/CoCoAManual/;
  /bin/mv $1/CoCoAManual/tex/CoCoAManual.pdf $1/CoCoAManual/.;
  /bin/cp ../../../doc/CoCoATranslationTable.html $1/CoCoAManual/.;
  /bin/rm -rf $1/CoCoAManual/CVS/;
  /bin/rm -rf $1/CoCoAManual/aux-files/;
  /bin/rm -rf $1/CoCoAManual/tex/;
  /bin/rm -rf $1/CoCoAManual/html/CVS;
}

copy-emacs()
{
  echo "Copying emacs in       ..${1:(-40)}/";  # last chars
  mkdir -p $1/emacs;
  /bin/cp ../emacs/cocoa5.emacs  $1/emacs/;
  /bin/cp ../emacs/cocoa5.el     $1/emacs/;
}

copy-emacs-el()
{
  echo "Copying emacs.el in    ..${1:(-40)}/";  # last chars
  mkdir -p $1/emacs;
  /bin/cp ../emacs/cocoa5.el     $1/emacs/;
}

MakeTGZ()
{
  CURR_DIR=`pwd` # just for printing
  echo "compressing  ..${CURR_DIR:(-40)}/$1";
  tar -czf $1-$2.tgz $1;
}

MakeZIP()
{
  CURR_DIR=`pwd` # just for printing
  echo "compressing  ..${CURR_DIR:(-40)}/$1";
  zip -rq $1-$2.zip $1;
}

# suggest-scp()
# {
#   echo "scp $1-$2.tgz $WEBSITE/download/bin/cocoa-$REL-$2.tgz"
# }

suggest-sftp()
{
  echo "put $1-$2.$3 $WEBSITE/download/bin/cocoa-$REL-$2.$3"
}


