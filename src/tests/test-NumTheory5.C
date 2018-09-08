//   Copyright (c)  2015  John Abbott

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
#include "CoCoA/error.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/NumTheory-prime.H"

#include <iostream>
using std::cerr;
using std::endl;

namespace CoCoA
{

  // Test for CRTMill (very simple test)

  void program()
  {
    GlobalManager CoCoAFoundations;

    // Copied from ex-NumTheory2.C
    const BigInt N = power(10,100);
    const BigInt UPB = 2*N+1;

    CRTMill crt;
    int p = 101;
    while (true)
    {
      p = NextPrime(p);
      crt.myAddInfo(N%p, p); // tell crt the new residue-modulus pair
      if (CombinedModulus(crt) > UPB) break;
    }

    // Check answer is correct.
    if (CombinedResidue(crt) != N)
      CoCoA_ERROR("Wrong answer", "CoCoA::Program");
  }

} // end of namespace CoCoA


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
