#ifndef CoCoA_OrthogonalPolys_H
#define CoCoA_OrthogonalPolys_H

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
//   -------------------------------------------------------------------------

#include "CoCoA/ring.H"

namespace CoCoA
{

  RingElem ChebyshevPoly(long n, ConstRefRingElem x);  // 1st type
  RingElem ChebyshevPoly2(long n, ConstRefRingElem x); // 2nd type

  RingElem HermitePoly(long n, ConstRefRingElem x);    // physics
  RingElem HermitePoly2(long n, ConstRefRingElem x);   // probability

  RingElem LaguerrePoly(long n, ConstRefRingElem x);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/OrthogonalPolys.H,v 1.1 2017/10/16 19:53:49 abbott Exp $
// $Log: OrthogonalPolys.H,v $
// Revision 1.1  2017/10/16 19:53:49  abbott
// Summary: Added new fns ChebyshevPoly, HermitePoly, LaguerrePoly
//
//

#endif
