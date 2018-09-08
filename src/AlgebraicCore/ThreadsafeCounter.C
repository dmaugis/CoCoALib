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


#include "CoCoA/ThreadsafeCounter.H"
#include "CoCoA/error.H"

#include <iostream>
// using std::ostream;


namespace CoCoA
{

  long ThreadsafeCounter::myAdvance(long skip)
  {
    if (skip <= 0) CoCoA_ERROR(ERR::NotPositive, "ThreadsafeCounter::myAdvance");
    const long OrigCount = myCount;
    myCount += skip;
    return OrigCount;
  }


  std::ostream& operator<<(std::ostream& out, const ThreadsafeCounter& c)
  {
    if (!out) return out;  // short-cut for bad ostreams
    const long count = c.myCount;
    return out << "ThreadsafeCounter(" << count << ")";
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ThreadsafeCounter.C,v 1.5 2017/10/04 10:28:54 abbott Exp $
// $Log: ThreadsafeCounter.C,v $
// Revision 1.5  2017/10/04 10:28:54  abbott
// Summary: Removed dependency on BOOST
//
// Revision 1.4  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2015/11/04 10:29:13  abbott
// Summary: Added some C++11 threadsafe bits
//
// Revision 1.2  2012/07/19 17:12:48  abbott
// Added default **NON-THREADSAFE** version if BOOST is not available.
//
// Revision 1.1  2012/05/29 14:54:00  abbott
// Made ThreadsafeCounter separate from symbol.
//
//
