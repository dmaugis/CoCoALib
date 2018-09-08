#ifndef CoCoA_CpuTimeLimit_H
#define CoCoA_CpuTimeLimit_H

//   Copyright (c)  2017,2018  John Abbott, Anna Bigatti

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


#include "CoCoA/error.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class TimeoutException: public ErrorInfo
  {
  public:
    TimeoutException(const char* const context): ErrorInfo(ERR::TimedOut, context) {}
//    TimeoutException(const std::string& context): ErrorInfo(ERR::TimedOut, context) {}
    virtual ~TimeoutException();
    virtual std::ostream& myOutputSelf(std::ostream& out) const;
  private: // same data members as class ErrorInfo
  };

  //------------------------------------------------------------------

  class CpuTimeLimit
  {
  public:
    explicit CpuTimeLimit(double interval);
    void operator()(const  char* const context) const { if (/*!IamUnlimited() &&*/ ++myCountSinceLastCheck == myDoCheckAtThisCount && IamTimedOut()) throw TimeoutException(context); }
    void myPrepareForNewLoop() const { myCountSinceLastCheck = 0; if (!IamUnlimited()) myDoCheckAtThisCount = 1; }
    bool IamUnlimited() const { return (myDoCheckAtThisCount == 0); };
  private: // data members
    mutable long myCountSinceLastCheck;
    mutable long myDoCheckAtThisCount;  // 0 means unlimited, o/w always positive
    mutable double myLastCheckTime;
    double myTargetCheckInterval;
    double myTriggerTime;
    mutable long myTotalCount; // not used currently
    friend std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& );
  private: // impl details
    bool IamTimedOut() const;

    enum NO_LIMIT_t { NO_LIMIT };
    explicit CpuTimeLimit(NO_LIMIT_t);
    friend const CpuTimeLimit& NoCpuTimeLimit();
  };

  std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& TimeLimit);

  const CpuTimeLimit& NoCpuTimeLimit();

  inline bool IsUnlimited(const CpuTimeLimit& TimeLimit)
  { return TimeLimit.IamUnlimited(); }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/CpuTimeLimit.H,v 1.7 2018/06/27 10:20:16 abbott Exp $
// $Log: CpuTimeLimit.H,v $
// Revision 1.7  2018/06/27 10:20:16  abbott
// Summary: Updated
//
// Revision 1.6  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.5  2017/09/06 14:02:16  abbott
// Summary: Added new ERR::TimedOut; TimeoutException derives from ErrorInfo
//
// Revision 1.4  2017/07/22 12:56:34  abbott
// Summary: Updated exception type
//
// Revision 1.3  2017/07/21 15:06:09  abbott
// Summary: Major revision -- no longer needs BOOST
//
// Revision 1.2  2017/07/21 13:21:21  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.1  2017/07/15 15:18:43  abbott
// Summary: Added CpuTimeLimit
//
//

#endif
