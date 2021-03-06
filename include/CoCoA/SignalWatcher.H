#ifndef CoCoA_SignalWatcher_H
#define CoCoA_SignalWatcher_H

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
//   -------------------------------------------------------------------------

#include "CoCoA/interrupt.H"

#include <iosfwd>
//using std::ostream;
#include <csignal>
//using sigaction; (perhaps not in namespace std?)

namespace CoCoA
{

  class InterruptedBySignal: public InterruptReceived
  {
  public:
    InterruptedBySignal(int sig, const char* const context): InterruptReceived(context), mySignal(sig) {}
    virtual ~InterruptedBySignal();
    virtual std::ostream& myOutputSelf(std::ostream& out) const;
  private: // data members (in addition to those of CoCoA::InterruptReceived)
    const int mySignal; // signal which caused the interrupt
    friend int TriggeringSignal(const InterruptedBySignal& intr);
  };

  inline int TriggeringSignal(const InterruptedBySignal& intr)
  { return intr.mySignal; }


  // RAII class: ctor sets new sigaction, dtor restores previous one.
  // NB myDeactivate restores previous sigaction.
  class SignalWatcher
  {
  public:
    explicit SignalWatcher(int sig, void FnPtr(int) = 0/*nullptr*/);
    ~SignalWatcher();
    void myOutputSelf(std::ostream& out) const;
    friend void deactivate(SignalWatcher& SW) { SW.myDeactivate(); }
  private:
    bool IamActive() const { return (myPrevSigactionPtr != 0/*nullptr*/); }
    void myDeactivate();
  private: // data members
    int mySig;
    struct sigaction* myPrevSigactionPtr;
  };

  void deactivate(SignalWatcher& SW);
  std::ostream& operator<<(std::ostream& out, const SignalWatcher& w);

  int GetAndResetSignalReceived();
//  void ResetInterruptSignalReceived();
  void SetSignalReceived(int sig); // used by interpreter & C5 (gui)

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SignalWatcher.H,v 1.3 2017/07/23 15:28:46 abbott Exp $
// $Log: SignalWatcher.H,v $
// Revision 1.3  2017/07/23 15:28:46  abbott
// Summary: Renamed SetInterruptSignalReceived to SetSignalReceived
//
// Revision 1.2  2017/07/22 13:00:59  abbott
// Summary: Added new exception class InterruptedBySignal; some cleaning
//
// Revision 1.1  2017/07/21 13:21:22  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
//

#endif
