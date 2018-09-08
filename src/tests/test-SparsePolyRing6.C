//   Copyright (c)  2017  John Abbott

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
#include "CoCoA/OrthogonalPolys.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-SturmSeq.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test checks that NumRealRoots gives the correct answer;
// internally NumRealRoots calls SturmSeq, so hopefully this also
// also checks that SturmSeq is correct.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    const ring P = NewPolyRing(RingQQ(), symbols("x"));
    RingElem f(P, "x^4+x^3-x-1");
    CoCoA_ASSERT_ALWAYS(NumRealRoots(f) == 2);
    
    RingElem g = ChebyshevPoly(50, RingElem(P,"x"));
    CoCoA_ASSERT_ALWAYS(NumRealRoots(g) == 50);
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
