# Makefile configuration for CoCoALib.
# Created automatically by the configure script.
# Created on  2018-08-09  at time  10:16:36
# Command was: 
# ./configure  "--with-libgmp=/Users/bigatti/0.99/gmp-6.1.2/.libs/libgmp.a"

PLATFORM="Darwin 15.6.0 x86_64"

##################################################################
# Definitions common to all Makefiles in CoCoALib.

# Version number of CoCoALib we shall build.
VERSION=0.99600

INSTALL_CMD=install
COCOALIB_INSTALL_DIR=/usr/local
COCOA5_INSTALL_DIR=/usr/local/bin

EXTLIBS=$(COCOA_ROOT)/configuration/ExternalLibs

# Compilation flags:
CXXFLAGS=-std=c++03 -Wall -pedantic -fPIC  -O2  


######################################################
# These variables were set by the configure script.

CXX=g++

# We use the following GMP installation:
GMP_VERSION=6.1.2
GMP_LIB=$(EXTLIBS)/lib/libgmp-symlink.a
GMP_LDLIB=-lgmp-symlink
GMPXX_LIB=$(EXTLIBS)/lib/libgmpxx-symlink.a
GMPXX_LDLIB=-lgmpxx-symlink

HAVE_QMAKE=no

# BOOST settings:
HAVE_BOOST=yes
BOOST_LDLIBS=

# READLINE settings:
HAVE_READLINE=yes
READLINE_LDLIBS=-lreadline

# CDD settings:
HAVE_CDD=no

# FROBBY settings:
HAVE_FROBBY=no

# GFAN settings:
HAVE_GFAN=no

# GSL settings:
HAVE_GSL=no

# MathSAT settings:
HAVE_MATHSAT =no

# Normaliz settings:
HAVE_NORMALIZ=no

LDLIBS=-Wl,-no_pie $(COCOA_LIB) -L$(EXTLIBS)/lib  $(FROBBY_LDLIBS)  $(GFAN_LDLIBS)  $(CDD_LDLIBS)  $(GSL_LDLIBS)  $(MATHSAT_LDLIBS)  $(NORMALIZ_LDLIBS)  $(GMPXX_LDLIB)  $(GMP_LDLIB)

COCOA5_CXX_DEFINES=-DCoCoA_WITH_READLINE
COCOA5_LDLIBS=$(SOCKET_LIB)  $(BOOST_LDLIBS)  $(READLINE_LDLIBS)


##################################################################
# This is the second fixed part of the common Makefile definitions.
# The configure script will copy it to autoconf.mk.
# NOTE: COCOA_ROOT is defined as a relative path in each individual Makefile.

COCOA_HDR=$(COCOA_ROOT)/include/CoCoA/library.H

INCLUDE=-I$(COCOA_ROOT)/include  -isystem $(EXTLIBS)/include
COMPILE=$(CXX)  $(CXXFLAGS)  $(INCLUDE)
COCOA_LIB=$(COCOA_ROOT)/lib/libcocoa.a



# Rule for compiling C++ code in *.C files into *.o object files
%.o: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE) -c -o $@ $<

# Rule for compiling and linking C++ code in *.C files
%: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE) -o $@ $< $(LDLIBS)
	@AppleDir="$@.dSYM" ; \
	echo " " $(CXXFLAGS) " " | fgrep " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


# Rule for compiling C++ code in *.cpp files into *.o object files
%.o: %.cpp
	@echo "Compiling `basename $@`"
	$(COMPILE) -c -o $@ $<

# Rule for compiling and linking C++ code in *.cpp files
%: %.cpp
	@echo "Compiling `basename $@`"
	$(COMPILE) -o $@ $< $(LDLIBS)
	@AppleDir="$@.dSYM" ; \
	echo " " $(CXXFLAGS) " " | fgrep " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


# The following are derived from the conventions for Makefiles for code
# which should become part of the GNU project.  It seems reasonable to
# adopt them here as well.  I found the recommendations in the online
# help for make (e.g. run the command "info make" on a GNU/Linux system).
SHELL=/bin/bash
.SUFFIXES:
.SUFFIXES: .C .o
