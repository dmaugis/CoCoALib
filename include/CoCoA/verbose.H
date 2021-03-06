#ifndef CoCoA_verbose_H
#define CoCoA_verbose_H

//   Copyright (c)  2016  John Abbott

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.
//   -------------------------------------------------------------------------

#include "CoCoA/config.H"

#include <iosfwd>
//using std::ostream;
#include <string>
//using std::string;

namespace CoCoA
{

  class VerboseLog
  {
  public:
    VerboseLog(const char* const FnName);
    ~VerboseLog();
    std::ostream& operator()(long level);
  private:
    VerboseLog(const VerboseLog&); // NEVER DEFINED
  private:
    const std::string myFnName;
    static long ourNestingDepth;
    static long ourVerbosityLevel;
    //    friend long SetVerbosityLevel(long);
    //    friend bool IsVerbosityLevel(long level);
    friend void SetVerbosityLevel(long);
    friend long VerbosityLevel();
  };

//   long SetVerbosityLevel(long NewLevel); // only log messages with level <= NewLevel will be printed
//   bool IsVerbosityLevel(long level); // true iff the (global) verbosity level is <= level
  void SetVerbosityLevel(long NewLevel); // only log messages with level <= NewLevel will be printed
  long VerbosityLevel(); // return current verbosity level

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/verbose.H,v 1.2 2017/03/02 10:04:22 bigatti Exp $
// $Log: verbose.H,v $
// Revision 1.2  2017/03/02 10:04:22  bigatti
// -- modified interface for VerbosityLevel
//
// Revision 1.1  2016/11/11 13:24:08  abbott
// Summary: new file for "verbose" capabilities
//
//

#endif
