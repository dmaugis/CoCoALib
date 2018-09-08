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


#include "CoCoA/QuotientRing.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


////////////////////// See redmine issue 951 comment 7 //////////////////////

namespace CoCoA
{

  void program()
  {
    // gcd call in last line used to SEGV.
    GlobalManager CoCoAFoundations;

    const ring F101 = NewZZmod(101);
    const ring P = NewPolyRing(F101, symbols("z"));
    const RingElem z = indet(P,0);
    const RingElem f = RingElem(P, "-19*z^111 -26*z^110 +23*z^109 -28*z^108 -36*z^107 +10*z^106 -47*z^105 +16*z^104 -19*z^103 -20*z^102 +20*z^101 -39*z^100 +45*z^99 +42*z^98 +15*z^97 -36*z^96 +40*z^95 -11*z^94 -24*z^93 -21*z^92 -8*z^91 +38*z^90 -44*z^89 -50*z^88 +13*z^87 +12*z^86 +37*z^85 +46*z^84 +37*z^83 -45*z^82 +35*z^81 -3*z^80 +39*z^79 -12*z^78 +45*z^77 +10*z^76 -22*z^75 -25*z^74 +20*z^73 -2*z^72 +18*z^71 +34*z^70 +41*z^69 -38*z^68 +19*z^67 -49*z^66 -48*z^65 +47*z^64 -20*z^63 -31*z^62 -8*z^61 -7*z^60 -4*z^59 -19*z^58 -18*z^57 +28*z^56 +9*z^55 +49*z^54 -13*z^53 +21*z^52 -14*z^51 +8*z^50 +22*z^49 -9*z^48 -18*z^47 -43*z^46 +48*z^45 -11*z^44 +21*z^43 +45*z^42 +20*z^41 -10*z^40 +37*z^39 -49*z^38 +5*z^37 -23*z^36 +38*z^35 -13*z^34 +12*z^33 +32*z^32 -44*z^31 +19*z^30 -25*z^29 -21*z^28 -36*z^27 -47*z^26 +19*z^25 +17*z^24 -30*z^23 +34*z^22 +43*z^21 +5*z^20 +3*z^19 +15*z^18 -36*z^17 +28*z^16 -49*z^15 +23*z^14 -32*z^13 -10*z^12 -50*z^11 -29*z^10 +30*z^9 -30*z^8 +5*z^7 +34*z^6 +44*z^5 -50*z^4 -44*z^3 +44*z^2");

    const RingElem Df = deriv(f, z);
    CoCoA_ASSERT_ALWAYS(monic(gcd(f,Df)) == z);
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
