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

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// This test checks that IsSqFree gives the expected answer over the
// range 1 to 100.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    for (int i=1; i < 100; ++i)
    {
      const factorization<long> FacInfo = factor(i);
      const vector<long>& mult = FacInfo.myMultiplicities();
      int MaxExp = 0;
      for (int j=0; j < len(mult); ++j)
        if (mult[j] > MaxExp) MaxExp = mult[j];
      const bool ans1 = IsSqFree(i);
      const bool3 ans2 = IsSqFree(BigInt(i));
      const bool ans3 = IsSqFree(-i);
      const bool3 ans4 = IsSqFree(BigInt(-i));
      CoCoA_ASSERT_ALWAYS(ans1 == (MaxExp <= 1));
      CoCoA_ASSERT_ALWAYS(ans1 == IsTrue3(ans2) && ans1 == !IsFalse3(ans2));
      CoCoA_ASSERT_ALWAYS(ans1 == ans3);
      CoCoA_ASSERT_ALWAYS(ans1 == IsTrue3(ans4) && ans1 == !IsFalse3(ans4));
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
