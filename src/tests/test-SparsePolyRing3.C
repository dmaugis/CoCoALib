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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/time.H"

#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::clog;
using std::endl;

//----------------------------------------------------------------------
// This test checks that radical of a polynomial produces the expected result.
// Coeffs are QQ or ZZ/(2) or ZZ/(32003)
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.


  void test(const ring& k)
  {
    ring P = NewPolyRing(k, symbols("x,y,z"));
    const RingElem f1 = RingElem(P, "x");
    const RingElem f2 = RingElem(P, "3*y");
    const RingElem f3 = RingElem(P, "5*x+7");
    const RingElem f4 = RingElem(P, "11*x+13*z");
    const RingElem f5 = RingElem(P, "17*x+19*y+23*z+29");

    // First test just on monomials
    for (int e1 = 0; e1 < 4; ++e1)    
    for (int e2 = 0; e2 < 4; ++e2)    
    {
      const RingElem f = power(f1,e1) * power(f2,e2);
      const RingElem expected = power(f1,min(e1,1)) * power(f2,min(e2,1));
      const RingElem radF = radical(f);
      CoCoA_ASSERT_ALWAYS(monic(radF) == monic(expected));
      CoCoA_ASSERT_ALWAYS(monic(radical(-f)) == monic(radF));
      CoCoA_ASSERT_ALWAYS(monic(radical(radF)) == monic(radF));
    }

    // More general test on polynomials
    for (int e1 = 0; e1 < 1; ++e1)    // effectively ignore f1 to make the test faster.
    for (int e2 = 0; e2 < 4; ++e2)    
    for (int e3 = 0; e3 < 4; ++e3)    
    for (int e4 = 0; e4 < 4; ++e4)    
    for (int e5 = 0; e5 < 3; ++e5)    
    {
      const RingElem f = power(f1,e1) * power(f2,e2) * power(f3,e3) * power(f4,e4) * power(f5,e5);
      const RingElem expected = power(f1,min(e1,1)) * power(f2,min(e2,1)) * power(f3,min(e3,1)) * power(f4,min(e4,1)) * power(f5,min(e5,1));
      const RingElem radF = radical(f);
      CoCoA_ASSERT_ALWAYS(monic(radF) == monic(expected));
      CoCoA_ASSERT_ALWAYS(monic(radical(-f)) == monic(radF));
      CoCoA_ASSERT_ALWAYS(monic(radical(radF)) == monic(radF));
    }
  }

  void program()
  {
    // You may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

   test(RingQQ());
   test(NewZZmod(2)); 
   test(NewZZmod(32003));
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
