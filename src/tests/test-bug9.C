//   Copyright (c)  2018  John Abbott

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


#include "CoCoA/RingQQ.H"
#include "CoCoA/factor.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


////////////////////// See redmine issue 1185 //////////////////////

namespace CoCoA
{

  void program()
  {
    // See redmine #1185: RemainingFactor used to have wrong sign
    GlobalManager CoCoAFoundations;

    const ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingElem f(P, "(y+z)*(y^2-x)");
    factorization<RingElem> facs = factor(f);
    const std::vector<RingElem>& fac = facs.myFactors();
    const std::vector<long>& mult = facs.myMultiplicities();
    CoCoA_ASSERT_ALWAYS(len(fac) == 2);
    CoCoA_ASSERT_ALWAYS(mult[0] == 1);
    CoCoA_ASSERT_ALWAYS(mult[1] == 1);

    RingElem g = facs.myRemainingFactor() * fac[0] * fac[1];
    CoCoA_ASSERT_ALWAYS(f == g);
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
