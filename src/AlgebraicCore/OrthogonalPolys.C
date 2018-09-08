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

#include "CoCoA/OrthogonalPolys.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/PolyRing.H"

#include <utility>
// using std::pair
using namespace std;


namespace CoCoA
{

  //-------------------------------------------------------
  // Chebyshev polys of 1st type
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem ChebyshevPoly_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c = power(2, n-1);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k*(n-k));
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  RingElem ChebyshevPoly_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = one(owner(x));
    RingElem curr = x;
    const RingElem twox = 2*x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem tmp = twox*curr - prev;
      swap(prev, curr);
      swap(curr, tmp);
    }
    return curr;
  }


  std::pair<RingElem, RingElem> ChebyshevPoly_recursive(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return std::make_pair(one(owner(x)), x);
    if (n == 1) return std::make_pair(x, 2*x*x-1);
    const long halfn = n/2; // integer division
    const std::pair<RingElem,RingElem> RecAns = ChebyshevPoly_recursive(halfn, x);
    const RingElem& Ta = RecAns.first;
    const RingElem& Tb = RecAns.second;
    if (IsEven(n))
      return make_pair(2*power(Ta,2)-1, 2*Ta*Tb - x);
    return make_pair(2*Ta*Tb-x, 2*power(Tb,2)-1);
  }


  // Dispatch function
  RingElem ChebyshevPoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "ChebyshevPoly: index must be >= 0");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return ChebyshevPoly_explicit(n, x);
    return ChebyshevPoly_recursive(n, x).first;
  }


  //-------------------------------------------------------
  // Chebyshev polys of 2nd type
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem ChebyshevPoly2_explicit(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "ChebyshevPoly: index must be >= 0");
    BigInt c = power(2, n);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k*(n-k+1));
      ans += c*power(x,n-2*k);
    }
    return ans;
  }

  RingElem ChebyshevPoly2_iter(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "ChebyshevPoly2: index must be >= 0");
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = one(owner(x));
    RingElem curr = 2*x;
    const RingElem twox = 2*x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem tmp = twox*curr - prev;
      swap(prev, curr);
      swap(curr, tmp);
    }
    return curr;
  }

  // Dispatch function
  RingElem ChebyshevPoly2(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "ChebyshevPoly2: index must be >= 0");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return ChebyshevPoly2_explicit(n, x);
    return ChebyshevPoly2_iter(n, x);
  }


  //=======================================================

  //-------------------------------------------------------
  // Hermite polys  (physics)
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem HermitePoly_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c = power(2, n);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(4*k);
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  // Physicists' Hermite polys:
  RingElem HermitePoly_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    const RingElem twox = 2*x;
    if (n == 1) return twox;
    RingElem prev = twox;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = twox*prev - deriv(prev,x);
      swap(curr, prev); // really just prev = curr;
    }
    return prev;
  }


  // Dispatch function
  RingElem HermitePoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "HermitePoly: index must be >= 0");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return HermitePoly_explicit(n, x);
    return HermitePoly_iter(n, x);
  }


  //-------------------------------------------------------
  // Hermite polys  (probability)
  //-------------------------------------------------------

  // Explicit formula from Wikipedia
  RingElem HermitePoly2_explicit(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    BigInt c(1);
    RingElem ans = c*power(x,n);
    for (int k=1; k <= n/2; ++k) // integer division!
    {
      c = -(c*(n-2*k+2)*(n-2*k+1))/(2*k);
      ans += c*power(x,n-2*k);
    }
    return ans;
  }


  RingElem HermitePoly2_iter(long n, ConstRefRingElem x)
  {
    CoCoA_ASSERT(n >= 0);
    if (n == 0) return one(owner(x));
    if (n == 1) return x;
    RingElem prev = x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = x*prev - deriv(prev,x);
      swap(curr, prev); // really just assignment prev = curr
    }
    return prev;
  }
  
  // Dispatch function
  RingElem HermitePoly2(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "HermitePoly2: index must be >= 0");
    const ring& P = owner(x);
    if (IsPolyRing(P) && IsIndet(x))
      return HermitePoly2_explicit(n, x);
    return HermitePoly2_iter(n, x);
  }

  //=======================================================

  //-------------------------------------------------------
  // Laguerre polys  (mult by factorial so coeffs are integer)
  //-------------------------------------------------------

  // Laguerre poly ***MULT BY FACTORIAL(N)***
  // Inefficient (roughly cubic in n)
  RingElem LaguerrePoly(long n, ConstRefRingElem x)
  {
    if (n < 0)
      CoCoA_ERROR(ERR::BadArg, "LaguerrePoly: index must be >= 0");
    if (n == 0) return one(owner(x));
    if (n == 1) return 1-x;
    RingElem prev2 = one(owner(x));
    RingElem prev1 = 1-x;
    for (int k = 2; k <= n; ++k)
    {
      RingElem curr = (2*k-1-x)*prev1 - power(k-1,2)*prev2;
      swap(prev1, prev2); // really just prev2 = prev1;
      swap(curr, prev1);  // really just prev1 = curr;
    }
    return prev1;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/OrthogonalPolys.C,v 1.4 2018/05/22 14:16:40 abbott Exp $
// $Log: OrthogonalPolys.C,v $
// Revision 1.4  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.3  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.2  2018/04/18 14:28:20  abbott
// Summary: Removed unused arg
//
// Revision 1.1  2017/10/16 19:53:50  abbott
// Summary: Added new fns ChebyshevPoly, HermitePoly, LaguerrePoly
//
//
