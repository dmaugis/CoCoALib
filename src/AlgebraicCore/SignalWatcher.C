//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/SignalWatcher.H"
#include "CoCoA/error.H"

#include <iostream>

namespace CoCoA
{

  namespace // anonymous
  {

    volatile sig_atomic_t SignalReceived = 0; // GLOBAL VARIABLE!!!
    
  } // end of anonymous namespace


  SignalWatcher::SignalWatcher(int sig, void FnPtr(int)):
      mySig(sig)
  {
    struct sigaction sa;
    if (FnPtr == 0/*nullptr*/)
      sa.sa_handler = &SetSignalReceived; // default CoCoA signal handler
    else
      sa.sa_handler = FnPtr;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART; // Restart functions if interrupted by handler

    myPrevSigactionPtr = new struct sigaction;
    if (sigaction(sig, &sa, myPrevSigactionPtr) != 0)
    {
      delete myPrevSigactionPtr;
      CoCoA_ERROR("Unable to set signal handler", "SignalWatcher ctor");
    }
  }


  void SignalWatcher::myOutputSelf(std::ostream& out) const
  {
    out << "SignalWatcher(sig=" << mySig;
    if (IamActive())
      out << ')';
    out << ", DEACTIVATED)";
  }


  void SignalWatcher::myDeactivate()
  {
    if (!IamActive()) return;
    sigaction(mySig, myPrevSigactionPtr, 0/*nullptr*/);
    myPrevSigactionPtr = 0/*nullptr*/;
  }


  SignalWatcher::~SignalWatcher()
  {
    myDeactivate();
  }


  std::ostream& operator<<(std::ostream& out, const SignalWatcher& SW)
  {
    SW.myOutputSelf(out);
    return out;
  }


  int GetAndResetSignalReceived()
  {
    const int sig = SignalReceived;
    SignalReceived = 0;
    return sig;
  }

  void SetSignalReceived(int sig)
  {
    // SANITY CHECK on value of sig???

// Next 2 lines are for an old version of CpuTimeLimit which used SIGVTALRM
//    // SIGVTALRM has lowest priority -- ignore it if another signal is waiting.
//    if (sig == SIGVTALRM && SignalReceived != 0) return;
    SignalReceived = sig;
  }


  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  InterruptedBySignal::~InterruptedBySignal()
  {}

  std::ostream& InterruptedBySignal::myOutputSelf(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CoCoA::InterruptedBySignal(signal=" << mySignal
        << ", context=\"" << myContext << "\")";
    return out;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SignalWatcher.C,v 1.3 2017/07/23 15:28:47 abbott Exp $
// $Log: SignalWatcher.C,v $
// Revision 1.3  2017/07/23 15:28:47  abbott
// Summary: Renamed SetInterruptSignalReceived to SetSignalReceived
//
// Revision 1.2  2017/07/22 13:01:02  abbott
// Summary: Added new exception class InterruptedBySignal; some cleaning
//
// Revision 1.1  2017/07/21 13:21:23  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
//
