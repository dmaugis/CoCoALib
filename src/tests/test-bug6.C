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
#include "CoCoA/ideal.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


////////////////////// See redmine issue 870 //////////////////////

namespace CoCoA
{

  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring Qxy = NewPolyRing(RingQQ(), symbols("x,y"));
    const RingElem x = indet(Qxy,0);
    const RingElem y = indet(Qxy,1);
    vector<RingElem> v;
    v.push_back(x); v.push_back(y);
    const ideal I1 = ideal(v);
    v.clear();
    v.push_back(x-1); v.push_back(y-1);
    const ideal I2 = ideal(v);
    // These two ideal should be equal:
    const ideal I3 = I1*I2;
    const ideal I4 = I2*I1;
    // Check explicitly that GBasis elems reduce to zero
    const vector<RingElem>& GB3 = GBasis(I3);
    const vector<RingElem>& GB4 = GBasis(I4);
    CoCoA_ASSERT_ALWAYS(len(GB3) >= 2);
    CoCoA_ASSERT_ALWAYS(len(GB4) >= 2);
    for (int i=0; i < len(GB3); ++i)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I3)));
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I4)));
    }
    for (int i=0; i < len(GB4); ++i)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I3)));
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I4)));
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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

