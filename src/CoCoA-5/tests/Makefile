# Makefile for CoCoALib/src/CoCoA-5/tests/ directory

COCOA_ROOT=../../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/CoCoA-5/tests/

TESTS=exsegv.cocoa5 \
      exbugs.cocoa5 \
      bug-EvalTwice.cocoa5 \
      ErrMesg.cocoa5 \
      AnonFunc.cocoa5 \
      demo-GeMiTo2011.cocoa5 \
      demo-Osaka2015.cocoa5 \
      test-manual.cocoa5 \
      whatiscocoa.cocoa5 \
      primary.cocoa5 \
      lecture-HF1.cocoa5 \
      lecture-HF2.cocoa5 \
      lecture-HF3.cocoa5 \
      lecture-HF4.cocoa5 \
      tut-CoCoLA1.cocoa5 \
      tut-CoCoLA2.cocoa5 \
      tut-CoCoLA3.cocoa5 \
      tut-CoCoLA4.cocoa5 \
      RealRoots.cocoa5 \
      ExtLibNORMALIZ.cocoa5 \
      test-ApproxSolve.cocoa5 \
      test-FactorAlgExt.cocoa5 \
      test-HomomorphismOps.cocoa5 \
      test-PrimaryDecomposition0.cocoa5 \
      test-implicit.cocoa5 \
      test-toric.cocoa5 \
      SourceAnna.cocoa5 \
      radical.cocoa5


.PHONY: all
all: ../CoCoAInterpreter check

.PHONY: lib
lib: library

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)

# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make RunTests" in the CoCoALib root
# directory.  Most of the work is done by the RunTests.sh script.
.PHONY: check
check: ../CoCoAInterpreter
	@./RunTests.sh $(TESTS) 2> /dev/null  # discard err mesgs from the shell script

.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  ./*.found  ./*.cerr  ./*~  ./.\#*

.PHONY: veryclean
veryclean: clean-local
	@echo "Verycleaned CoCoALib/$(CWD)"



$(EXECS): $(COCOA_LIB)


#############################################################################
# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/tests/Makefile,v 1.29 2018/03/15 14:19:58 bigatti Exp $
# $Log: Makefile,v $
# Revision 1.29  2018/03/15 14:19:58  bigatti
# -- test-radical now is last
#
# Revision 1.28  2017/12/18 22:09:31  bigatti
# -- some sorting (test-..)
#
# Revision 1.27  2017/12/18 22:02:02  bigatti
# -- renamed HomomorphismFns.cocoa5 into test-HomomorphismOps.cocoa5
#
# Revision 1.26  2017/05/17 06:29:36  bigatti
# -- renamed bug-evalIdeal into bug-EvalTwice
#
# Revision 1.25  2017/05/16 16:27:46  bigatti
# -- added bug-evalIdeal.cocoa5
#
# Revision 1.24  2017/05/11 10:27:06  bigatti
# -- added tests fr toric and ImplicitHypersurface
#
# Revision 1.23  2017/02/15 14:07:51  bigatti
# -- added tests for PrimaryDecomposition0 and FactorAlgExt (and MinPoly)
#
# Revision 1.22  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.21  2016/07/21 08:38:47  bigatti
# -- added test for system solving
#
# Revision 1.20  2016/07/20 11:32:02  bigatti
# -- added test for manual
#
# Revision 1.19  2016/05/18 12:28:55  abbott
# Summary: Added new target "lib": synonym for "library"
#
# Revision 1.18  2016/02/18 08:11:31  bigatti
# -- added tutorials for cocoa school in India
#
# Revision 1.17  2015/07/27 16:03:16  bigatti
# -- added demo Bigatti Osaka
#
# Revision 1.16  2014/07/28 14:58:56  abbott
# Summary: Improved targets (very)clean-local
# Author: JAA
#
# Revision 1.15  2014/07/22 12:17:00  bigatti
# -- added radical
#
# Revision 1.14  2014/07/04 14:42:14  bigatti
# -- added test for ker, IsSurjective, IsInjective...
#
# Revision 1.13  2014/05/02 15:38:15  abbott
# Summary: Added new test-suite exbugs.cocoa5
# Author: JAA
#
# Revision 1.12  2014/03/25 15:45:14  abbott
# Summary: Renamed lambdas.cocoa5 to AnonFunc.cocoa5
# Author: JAA
#
# Revision 1.11  2012/10/02 07:48:30  bigatti
# -- added test for Normaliz
#
# Revision 1.10  2012/05/24 14:08:35  bigatti
# -- added ErrMesg.cocoa5,out  for testing expected error messages
#
# Revision 1.9  2012/05/20 11:38:11  abbott
# Added new test for RealRoots (after unearthing a bug).
#
# Revision 1.8  2012/01/19 10:37:13  bigatti
# -- added primary
#
# Revision 1.7  2011/09/16 14:05:59  bigatti
# -- added whatiscocoa.cocoa5
#
# Revision 1.6  2011/05/26 14:54:06  bigatti
# -- added first 3 lectures for cocoaschool 2011
#
# Revision 1.5  2011/05/24 08:19:41  bigatti
# -- added lecture-HF4.cocoa5
#
# Revision 1.4  2011/05/13 17:04:31  abbott
# Removed EXECS by changing interface to RunTests.sh.
#
# Revision 1.3  2011/05/11 11:41:36  bigatti
# -- added demo-GeMiTo2011
#
# Revision 1.2  2011/05/11 11:12:20  bigatti
# -- added lambdas
#
# Revision 1.1  2011/05/10 16:45:02  bigatti
# -- first import
#
