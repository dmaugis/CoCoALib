//   Copyright (c)  2016,2018  John Abbott, Anna M. Bigatti

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


#include "CoCoA/GCDFreeBasis.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

// #include<vector>
using std::vector;

#include <iostream>

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {

    // Result is X/Y^k where k is chosen max poss.
    RingElem RemoveBiggestPower(RingElem& X, const RingElem& Y)
    {
      CoCoA_ASSERT(owner(X) == owner(Y));
      RingElem quot = zero(owner(X));
      while (true)
      {
        if (!IsDivisible(quot, X,Y)) return X;
        swap(X, quot); // really assignment X = quot;
      }
    }
  
  } // end of anonymous namespace


  GCDFreeBasis_RingElem::LCR GCDFreeBasis_RingElem::myLCR(RingElem A, RingElem B) const
  {
    CoCoA_ASSERT(owner(A) == owner(B));
    const RingElem g = gcd(A,B);
    if (IsInvertible(g)) return LCR(A, vector<RingElem>(), B);
    A /= g;
    B /= g;
    if (IsInvertible(A)) // special case A == 1
    {
      RemoveBiggestPower(B, g);
      struct LCR tmp = myLCR(g, B);
      if (!IsInvertible(tmp.myL)) tmp.myC.push_back(tmp.myL);
      return LCR(A, tmp.myC, tmp.myR);
    }
    if (IsInvertible(B)) // special case B == 1
    {
      RemoveBiggestPower(A, g);
      struct LCR tmp = myLCR(A, g);
      if (!IsInvertible(tmp.myR)) tmp.myC.push_back(tmp.myR);
      return LCR(tmp.myL, tmp.myC, B);
    }
    // General case  A != 1 && B != 1
    struct LCR LCR_Ag = myLCR(A, g);
    if (IsInvertible(LCR_Ag.myR)) { swap(LCR_Ag.myR, B); return LCR_Ag; }
    struct LCR LCR_RB = myLCR(LCR_Ag.myR, B);
    vector<RingElem>& C = LCR_Ag.myC;
    if (!IsInvertible(LCR_RB.myL)) C.push_back(LCR_RB.myL);
    C.insert(C.end(), LCR_RB.myC.begin(), LCR_RB.myC.end()); // concat
    swap(LCR_Ag.myR, LCR_RB.myR); // equiv LCR_Ag.myR = LCR_RB.myR  really assignment
    return LCR_Ag;
  }


  void GCDFreeBasis_RingElem::myRefineBasis(RingElem X)
  {
    CoCoA_ASSERT(myCoprimeBasis.empty() || owner(X) == owner(myCoprimeBasis[0]));
    if (IsZero(X) || IsInvertible(X)) return;
    const int sz = len(myCoprimeBasis);
    if (sz == 0) { myCoprimeBasis.push_back(X); return; } // BUG??? gcd(X,0)???
    vector<RingElem> NewBasis; NewBasis.reserve(sz+1);
    for (int i=0; i < sz; ++i)
    {
      if (IsInvertible(X)) { NewBasis.push_back(myCoprimeBasis[i]); continue; }
      struct LCR tmp = myLCR(myCoprimeBasis[i], X);
      swap(X, tmp.myR); // X = tmp.myR; really assignment
      if (!IsInvertible(tmp.myL)) NewBasis.push_back(tmp.myL);
      if (!tmp.myC.empty())
        NewBasis.insert(NewBasis.end(), tmp.myC.begin(), tmp.myC.end());
    }
    if (!IsInvertible(X)) NewBasis.push_back(X); // BUG????  gcd(X,0)
    swap(myCoprimeBasis, NewBasis);
  }

// #if 0
//   // Alternative impl, should copy fewer values, but is slower ?!?
//   void GCDFreeBasis_RingElem::myRefineBasis(RingElem N)
//   {
//     CoCoA_ASSERT(N >= 0);
//     if (IsZero(N) || IsOne(N)) return;
//     const int sz = len(myCoprimeBasis);
//     if (sz == 0) { myCoprimeBasis.push_back(N); return; }
//     vector<RingElem> NewElems;
//     for (int i=0; i < sz; ++i)
//     {
//       if (IsOne(N)) break;
//       struct LCR tmp = myLCR(myCoprimeBasis[i], N);
//       if (IsOne(tmp.myL) && IsOne(tmp.myR)) { swap(myCoprimeBasis[i],tmp.myC[0])/*really assignment*/; break; }
//       if (!IsOne(tmp.myL))
//         swap(myCoprimeBasis[i], tmp.myL); /*really assignment*/
//       else
//       {
//         swap(myCoprimeBasis[i], tmp.myC.back()); /*really assignment*/
//         tmp.myC.resize(tmp.myC.size()-1); // delete last elem
//       }
//       swap(N, tmp.myR); // N = tmp.myR; (really assignment)
//       if (!tmp.myC.empty())
//         NewElems.insert(NewElems.end(), tmp.myC.begin(), tmp.myC.end()); // append
//     }
//     if (!IsOne(N)) { myCoprimeBasis.push_back(N); /*std::clog << '+';*/ }
//     if (!NewElems.empty())
//     {
//       /*std::clog<<"New elems: " << len(NewElems) << endl;*/
//       myCoprimeBasis.insert(myCoprimeBasis.end(), NewElems.begin(), NewElems.end());
//     }
//   }
// #endif


  void GCDFreeBasis_RingElem::myAddInfo(const RingElem& X)
  {
    if (!IsTrueGCDDomain(owner(X))) CoCoA_ERROR(ERR::NotTrueGCDDomain, "GCDFreeBasis_RingElem::myAddInfo");
    if (!myCoprimeBasis.empty() && owner(X) != owner(myCoprimeBasis[0]))
      CoCoA_ERROR(ERR::MixedRings, "GCDFreeBasis_RingElem::myAddInfo");
    myRefineBasis(X);
  }
  

  void GCDFreeBasis_RingElem::myAddInfo(const std::vector<RingElem>& v)
  {
    const int sz = len(v);
    for (int i=0; i < sz; ++i)
      myRefineBasis(v[i]);
  }


  std::ostream& operator<<(std::ostream& out, const GCDFreeBasis_RingElem& gcdfb)
  {
    if (!out) return out;
    out << "GCDFreeBasis_RingElem object";
    return out;
  }



// OLD IMPL
  // // std::vector<RingElem> GCDFreeBasis(const std::vector<RingElem>& L)
  // // {
  // //   // assume L is list of elems of ring with GCD
  // //   if (L.empty()) CoCoA_ERROR(ERR::Empty, "GCDFreeBasis");
  // //   const int n = len(L);
  // //   const ring& R = owner(L[0]);
  // //   vector<RingElem> ans;
  // //   ans.push_back(gcd(L[0], zero(R)));
  // //   if (n == 1) return ans;
  // //   for (int i=1; i < n; ++i)
  // //   {
  // //     const RingElem RemainingFactor = RefineGCDFreeBasis(ans, L[i]);
  // //     if (IsInvertible(RemainingFactor)) continue;
  // //     ans.push_back(gcd(RemainingFactor, zero(R)));
  // //   }
  // //   return ans;
  // // }

  // // struct GFB2
  // // {
  // //   GFB2(const RingElem& A, const RingElem& B): Afactor(A), Bfactor(B), GFBCommonFactors(0) { CoCoA_ASSERT(IsCoprime(A,B));}
  // //   RingElem Afactor;
  // //   RingElem Bfactor;
  // //   vector<RingElem> GFBCommonFactors;
  // // };

  // // GFB2 GCDFreeBasis2(RingElem A, RingElem B)
  // // {
  // //   const RingElem g = gcd(A,B);
  // //   if (IsOne(g)) return GFB2(A,B);
  // //   A /= g;
  // //   B /= g;
  // //   if (IsInvertible(A))
  // //   {
  // //     while (true)
  // //     {
  // //       RingElem quot;
  // //       std::clog << "Test DIV  B=" << B << "  in ring " << owner(B) << std::endl
  // //                 << "by g=" << g << "  in ring " << owner(g) << std::endl;
  // //       if (!IsDivisible(quot, B, g)) break;
  // //       swap(B,quot); // cheap assignment (maybe std::move?)
  // //     }
  // //     GFB2 ans = GCDFreeBasis2(g, B);
  // //     if (!IsInvertible(ans.Afactor)) ans.GFBCommonFactors.push_back(ans.Afactor);
  // //     ans.Afactor = A;
  // //     return ans;
  // //   }
  // //   if (IsInvertible(B))
  // //   {
  // //     while (true)
  // //     {
  // //       RingElem quot;
  // //       if (!IsDivisible(quot, A,g)) break;
  // //       swap(A,quot); // cheap assignment (maybe std::move?)
  // //     }
  // //     GFB2 ans = GCDFreeBasis2(g, A);
  // //     if (!IsInvertible(ans.Bfactor)) ans.GFBCommonFactors.push_back(ans.Bfactor);
  // //     ans.Bfactor = B;
  // //     return ans;
  // //   }
  // //   // general case: A, B are non-trivial
  // //   GFB2 ans = GCDFreeBasis2(A, g);
  // //   if (IsInvertible(ans.Bfactor)) { ans.Bfactor = B; return ans; }
  // //   GFB2 tmp = GCDFreeBasis2(ans.Bfactor, B);
  // //   if (!IsInvertible(tmp.Afactor)) ans.GFBCommonFactors.push_back(tmp.Afactor);
  // //   /// concat -- SLUG SLUG SLUG, do this better!
  // //   for (int i=0; i < len(tmp.GFBCommonFactors); ++i)
  // //     ans.GFBCommonFactors.push_back(tmp.GFBCommonFactors[i]);
  // //   ans.Bfactor = tmp.Bfactor;
  // //   return ans;
  // // }

  // // RingElem RefineGCDFreeBasis(std::vector<RingElem>& basis, RingElem N)
  // // {
  // //   const ring& R = owner(N);
  // //   if (basis.empty()) { basis.push_back(gcd(N,zero(R))); return one(R); }
  // //   CoCoA_ASSERT(owner(basis[0]) == R);
  // //   RingElem RemainingFactor = one(R);
  // //   vector<RingElem> NewBasis;
  // //   const int n = len(basis);
  // //   for (int i=0; i < n; ++i)
  // //   {
  // //     if (IsInvertible(N)) { NewBasis.push_back(basis[i]); continue; }
  // //     GFB2 tmp = GCDFreeBasis2(basis[i], N);
  // //     N = tmp.Bfactor;
  // //     if (!IsInvertible(tmp.Afactor)) NewBasis.push_back(tmp.Afactor);
  // //     // concat  SLUG SLUG SLUG
  // //     for (int j=0; j < len(tmp.GFBCommonFactors); ++j)
  // //       NewBasis.push_back(tmp.GFBCommonFactors[j]);
  // //   }
  // //   swap(basis, NewBasis);
  // //   return N;
  // // }
  
} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/GCDFreeBasis.C,v 1.2 2018/06/25 12:31:13 abbott Exp $
// $Log: GCDFreeBasis.C,v $
// Revision 1.2  2018/06/25 12:31:13  abbott
// Summary: Restructured code
//
// Revision 1.1  2017/02/01 10:36:49  abbott
// Summary: IMpl of GCDFreeBasis (transl from CoCoA-5)
//
//
