#ifndef CoCoA_ThreadsafeCounter_H
#define CoCoA_ThreadsafeCounter_H

//   Copyright (c)  2012, 2017  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#include <iosfwd>
//using std::ostream;

#if __cplusplus > 199711L  &&  defined(CoCoA_THREADSAFE)
#include <atomic>
#endif

namespace CoCoA
{

  // IMPORTANT: NOT THREADSAFE IF your compiler is pre C++11
  class ThreadsafeCounter
  {
  public:
    ThreadsafeCounter(): myCount(0) {}
    long myAdvance(long n);  // returns previous value of counter
    friend std::ostream& operator<<(std::ostream& out, const ThreadsafeCounter& c);
  private: // data members
#if __cplusplus > 199711L  &&  defined(CoCoA_THREADSAFE)
    std::atomic<long> myCount;
#else
    long myCount;
#endif
  };

  std::ostream& operator<<(std::ostream& out, const ThreadsafeCounter& c);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ThreadsafeCounter.H,v 1.6 2017/10/04 10:33:28 abbott Exp $
// $Log: ThreadsafeCounter.H,v $
// Revision 1.6  2017/10/04 10:33:28  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2017/10/04 10:23:19  abbott
// Summary: Removed dependency on BOOST libs
//
// Revision 1.4  2017/03/13 12:17:09  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.3  2015/11/04 10:29:49  abbott
// Summary: Added some C++11 threadsafe bits
//
// Revision 1.2  2012/07/19 17:12:48  abbott
// Added default **NON-THREADSAFE** version if BOOST is not available.
//
// Revision 1.1  2012/05/29 14:54:00  abbott
// Made ThreadsafeCounter separate from symbol.
//
//

#endif
