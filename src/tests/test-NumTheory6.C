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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/error.H"

#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;

//----------------------------------------------------------------------
// This test checks that radical for integers gives the expected answer.
// It assumes that power(210,3) fits into a long.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  // Ad hoc impl -- no arg checks!
  long pwr(long base, int exp)
  {
    long ans = 1;
    for (int i=0; i < exp; ++i)
      ans *= base;
    return ans;
  }

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    const long p1 = 2;
    const long p2 = 3;
    const long p3 = 5;
    const long p4 = 7;
    for (int e1 = 0; e1 < 4; ++e1)    
    for (int e2 = 0; e2 < 4; ++e2)    
    for (int e3 = 0; e3 < 4; ++e3)    
    for (int e4 = 0; e4 < 4; ++e4)    
    {
      const long N = pwr(p1,e1) * pwr(p2,e2) * pwr(p3,e3) * pwr(p4,e4);
      const long radN = pwr(p1,min(e1,1)) * pwr(p2,min(e2,1)) * pwr(p3,min(e3,1)) * pwr(p4,min(e4,1));
      CoCoA_ASSERT_ALWAYS(radical(N) == radN);
      CoCoA_ASSERT_ALWAYS(radical(-N) == radN);
      CoCoA_ASSERT_ALWAYS(radical(radical(N)) == radN);
    }
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
