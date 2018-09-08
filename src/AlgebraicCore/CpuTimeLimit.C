//   Copyright (c)  2017,2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/error.H"
#include "CoCoA/time.H"

#include <iostream>
using namespace std;

namespace CoCoA
{

  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  TimeoutException::~TimeoutException()
  {}


  std::ostream& TimeoutException::myOutputSelf(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CoCoA::TimeoutException(context=\"" << myContext << "\")";
    return out;
  }
    
  //------------------------------------------------------------------
  
  CpuTimeLimit::CpuTimeLimit(double interval):
      myCountSinceLastCheck(0),
      myDoCheckAtThisCount(1),
      myLastCheckTime(CpuTime()),
      myTargetCheckInterval(interval/10),
      myTriggerTime(myLastCheckTime + interval),
      myTotalCount(0)
  {
    if (interval < 0) CoCoA_ERROR(ERR::NotNonNegative, "CpuTimeLimit ctor");
    if (interval > 1000000) CoCoA_ERROR(ERR::ArgTooBig, "CpuTimeLimit ctor");
  }


  bool CpuTimeLimit::IamTimedOut() const
  {
//    if (IamUnlimited()) return false;  // Should never happen!
    myTotalCount += myCountSinceLastCheck;
    myCountSinceLastCheck = 0;
    const double t = CpuTime();
    if (t > myTriggerTime) return true;
    const double ratio = (t-myLastCheckTime)/myTargetCheckInterval;
    myLastCheckTime = t;
    if (ratio >= 2 && myDoCheckAtThisCount > 1) { myDoCheckAtThisCount /= 2; }
    if (ratio < 0.5 && myDoCheckAtThisCount < 1000000) { myDoCheckAtThisCount *= 2; }
    return false;
  }


  // Quick makeshift impl.
  std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& TimeLimit)
  {
    if (!out) return out;  // short-cut for bad ostreams
    
    if (IsUnlimited(TimeLimit)) return out << "CpuTimeLimit(UNLIMITED)";

    out << "CpuTimeLimit(TriggerTime=" << TimeLimit.myTriggerTime
        << ", CurrTime=" << CpuTime()
        << ",  CountSinceLastCheck=" << TimeLimit.myCountSinceLastCheck
        << ", DoCheckAtThisCount=" << TimeLimit.myDoCheckAtThisCount << ")";
    return out;
  }


  // Special ctor for "unlimited" CpuTimeLimit object;
  // called only by NoCpuTimeLimit (below).
  CpuTimeLimit::CpuTimeLimit(NO_LIMIT_t):
      myCountSinceLastCheck(0),
      myDoCheckAtThisCount(0),  // "IMPOSSIBLE VALUE" means no limit
      myLastCheckTime(0),
      myTargetCheckInterval(0),
      myTriggerTime(0),
      myTotalCount(0)
    {}
  

  const CpuTimeLimit& NoCpuTimeLimit()
  {
    static const CpuTimeLimit SingleCopy(CpuTimeLimit::NO_LIMIT);
    return SingleCopy;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/CpuTimeLimit.C,v 1.11 2018/06/27 10:20:16 abbott Exp $
// $Log: CpuTimeLimit.C,v $
// Revision 1.11  2018/06/27 10:20:16  abbott
// Summary: Updated
//
// Revision 1.10  2018/06/27 09:37:57  abbott
// Summary: More detailed printout (more helpful for debugging)
//
// Revision 1.9  2018/06/25 12:31:49  abbott
// Summary: Added overflow protection when increasing interval length
//
// Revision 1.8  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.7  2017/09/06 14:08:50  abbott
// Summary: Changed name to TimeoutException
//
// Revision 1.6  2017/07/23 15:32:32  abbott
// Summary: Fixed STUPID bug in myDeactivate
//
// Revision 1.5  2017/07/22 13:03:02  abbott
// Summary: Added new exception InterruptdByTimeout; changed rtn type of myOutputSelf
//
// Revision 1.4  2017/07/21 15:06:10  abbott
// Summary: Major revision -- no longer needs BOOST
//
// Revision 1.3  2017/07/21 13:21:22  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.2  2017/07/15 15:44:44  abbott
// Summary: Corrected error ID name (to ArgTooBig)
//
// Revision 1.1  2017/07/15 15:17:48  abbott
// Summary: Added CpuTimeLimit
//
//
