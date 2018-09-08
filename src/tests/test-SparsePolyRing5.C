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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-RealRadical.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test checks the function RealRadical in some simple cases.
// It also checks HasRealRoots3.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingElem f1(P, "(x^2+1)*(x^4+4)*(x^5-5)");
    CoCoA_ASSERT_ALWAYS(RealRadical(f1) == RingElem(P, "x^5-5"));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f1) == RingElem(P, "x^5-5"));

    RingElem f2(P, "x^2+y^2+1");
    CoCoA_ASSERT_ALWAYS(RealRadical(f2) == one(P));
    CoCoA_ASSERT_ALWAYS(RealRadical(f2*f2) == one(P));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f2) == RingElem(P, "x^5-5"));

    RingElem f3(P, "x^2+y^2");
    CoCoA_ASSERT_ALWAYS(RealRadical(f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f3*f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f2*f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f3) == f3*RingElem(P,"x^5-5"));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f2*f3) == f3*RingElem(P,"x^5-5"));


    RingElem g1(P, "x^2+y^4+z^6");
    CoCoA_ASSERT_ALWAYS(IsTrue3(HasRealRoot3(g1)));
    CoCoA_ASSERT_ALWAYS(IsTrue3(HasRealRoot3(g1-1)));
    CoCoA_ASSERT_ALWAYS(IsFalse3(HasRealRoot3(g1+1)));
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
