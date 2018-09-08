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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


////////////////////// See redmine issue 858 //////////////////////

namespace CoCoA
{

  void program()
  {
    // Check that floor for TwinFloats does not throw ERR::ShouldNeverGetHere
    GlobalManager CoCoAFoundations;

    const ring RR = NewRingTwinFloat(64);
    const RingElem x = one(RR);
    for (int k = 1; k < 200; ++k)
    {
      const RingElem smaller = x - power(BigRat(1,2), k); // k=145 used to trigger ERR::ShouldNeverGetHere in floor
      try
      {
        const BigInt N = floor(smaller);
        CoCoA_ASSERT_ALWAYS(IsZero(N) || IsOne(N));
      }
      catch (const RingTwinFloat::InsufficientPrecision& err)
      {}
    }
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
