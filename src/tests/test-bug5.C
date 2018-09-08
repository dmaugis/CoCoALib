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
#include "CoCoA/BigRat.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/matrix.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // This test checks that a matrix with non integer entries
    // is rejected by NewMatrixOrdering.
    GlobalManager CoCoAFoundations;

    matrix M = NewDenseMat(RingQQ(),2,2);
    SetEntry(M,0,0,BigRat(1,2));
    SetEntry(M,0,1,BigRat(1,1));
    SetEntry(M,1,0,BigRat(3,2));
    SetEntry(M,1,1,BigRat(2,1));

    try
    {
      PPOrdering ord = NewMatrixOrdering(M,1);
      CoCoA_ASSERT_ALWAYS(false && "SHOULD NEVER GET HERE!");
    }
    catch (const ErrorInfo& err)
    {
      // Dispose of expected error; o/w rethrow
      if (err != ERR::BadArg) throw;
    }

    // This test checks that a matrix over a ring of non-zero
    // characteristic is rejected by Newmatrixordering.
    matrix M2 = NewDenseMat(NewZZmod(5),2,2);
    SetEntry(M,0,0,1);
    SetEntry(M,0,1,1);
    SetEntry(M,1,0,1);
    SetEntry(M,1,1,2);

    try
    {
      PPOrdering ord = NewMatrixOrdering(M2,1);
      CoCoA_ASSERT_ALWAYS(false && "SHOULD NEVER GET HERE!");
    }
    catch (const ErrorInfo& err)
    {
      // Dispose of expected error; o/w rethrow
      if (err != ERR::BadRing) throw;
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
