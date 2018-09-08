//   Copyright (c)  2017-2018  John Abbott,  Anna Bigatti
//   Original authors: Marvin Brandenstein, Alice Moallemy and Carsten Dettmar

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


#include "CoCoA/SparsePolyOps-RadicalMembership.H"

#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/apply.H"
#include "CoCoA/factor.H"
#include "CoCoA/ideal.H"
#include "CoCoA/interrupt.H"

#include <vector>
using std::vector;


namespace CoCoA
{

  namespace // anonymous
  {

    vector<RingElem> RadicalHelpers(const vector<RingElem>& G)
    {
      const long n = len(G);
      vector<RingElem> ans;
      for (long i=0; i < n; ++i)
      {
        if (NumTerms(G[i]) > 100) continue; // might take too long? (heuristic)
        factorization<RingElem> FacInfo = SqFreeFactor(G[i]);
        int nfacs = len(FacInfo.myMultiplicities());
        int MaxMult = 1;
        for (int j=0; j < nfacs; ++j)
        {
          if (MaxMult < FacInfo.myMultiplicities()[j])
            MaxMult = FacInfo.myMultiplicities()[j];
        }
        if (MaxMult == 1) continue;
        RingElem NewPoly = FacInfo.myFactors()[0];
        for (int j=1; j < nfacs; ++j)
          NewPoly *= FacInfo.myFactors()[j];
        ans.push_back(NewPoly);
      }
      return ans;
    }


    bool RabinovichTrick(ConstRefRingElem f, ideal I)
    {
      // Make 2 attempts to spot polys in radical(I):
      {
        const vector<RingElem> ExtraGens = RadicalHelpers(gens(I));
        if (!ExtraGens.empty())
        {
          I += ideal(ExtraGens);
        }
      }
      {
        const vector<RingElem> ExtraGens = RadicalHelpers(ReducedGBasis(I));
        if (!ExtraGens.empty())
        {
          I += ideal(ExtraGens);
        }
      }
      
      if (IsElem(f,I)) return true; // simple special case
      PolyRing P = RingOf(I);
      ring R = CoeffRing(P);
      PolyRing RabinovichPolyRing = NewPolyRing(R, NumIndets(P)+1);
      const RingElem& t = indet(RabinovichPolyRing,NumIndets(P)); 
      vector<RingElem> image = indets(RabinovichPolyRing);
      image.pop_back();
      RingHom Phi = PolyAlgebraHom(P, RabinovichPolyRing, image);
      vector<RingElem> newGens = apply(Phi, gens(I));
      newGens.push_back(1-t*Phi(f));
      return IsOne(ideal(newGens)); 
    }

  } // end of namespace anonymous
  

  bool IsInRadical(ConstRefRingElem f, const ideal& I)
  {
    const char* const FnName = "IsInRadical";
    if (owner(f) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_ERROR(ERR::NotPolyRing, FnName); // BUG BUG err should say "polyring over field"

    if (!IsHomog(I))
      return RabinovichTrick(f, I);

    RingElem g = f;
    while (!IsZero(g))
    {
      if (!RabinovichTrick(CutLF(g), I))
        return false;
    }
    return true;
  }

  bool IsInRadical(const ideal& I, const ideal& J)
  {
    const vector<RingElem>& G = gens(I);
    const long n = len(G);
    for (int i=0; i < n; ++i)
    {
      if (!IsInRadical(G[i], J)) return false;
    }
    return true;
  }


  // // Same as IsInRadical, without check for homogeneity
  // bool IsInRadical2(ConstRefRingElem f, const ideal& I)
  // {
  //   return RabinovichTrick(f,I);
  // }

  // // Same as IsInRadical, without check for homogeneity
  // bool IsInRadical2(const ideal& I, const ideal& J)
  // {
  //   for (int i=0; i < len(gens(I)); i++)
  //   {
  //     if (!IsInRadical2(gens(I)[i], J)) return false;
  //   }
  //   return true;
  // }
  

// MinPowerInIdeal: Determines the minimal power m such that f^m is in the Ideal J
  long MinPowerInIdeal_naive(ConstRefRingElem f, const ideal& I)
  {
    const char* const FnName = "MinPowerInIdeal";
    if (owner(f) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_ERROR(ERR::NotPolyRing, FnName); // BUG BUG err should say "polyring over field"

    if (!IsInRadical(f, I)) return -1;
    long D = 1;
    RingElem FpowD = f;
    while (!IsElem(FpowD, I))
    {
      CheckForInterrupt(FnName);
      FpowD *= f;
      ++D;
    }
    return D;
  }
  

// Exercise 5b): Alternative for MinPowerInIdeal using normal form
  long MinPowerInIdeal(RingElem f, const ideal& I)
  {
    const char* const FnName = "MinPowerInIdeal";
    if (owner(f) != RingOf(I))
      CoCoA_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_ERROR(ERR::NotPolyRing, FnName); // BUG BUG err should say "polyring over field"

    f = NF(f,I);
    if (!IsInRadical(f, I)) return -1;
    long D = 1;
    RingElem FpowD = f;
    while (!IsZero(FpowD))
    {
      CheckForInterrupt(FnName);
      FpowD = NF(FpowD*f, I);
      ++D;
    }
    return D;
  }


// // HomogComp: Splits the polynomial f into its homogeneous components
//   vector<RingElem> HomogComp(RingElem f)
//   {
//     vector<RingElem> components;
//     while (!IsZero(f))
//       components.push_back(CutLF(f));

//     return components;
//   }
  

// // following three methods only used for conjectures 
//   long MaxElem(const vector<long>& v)
//   {
//     long max = v[0];
//     for (int i = 1; i < len(v); ++i)
//     {
//       if (v[i] >= max)
//       {
//         max = v[i];
//       }
//     }
//     return max;
//   }


//   long Prod(const vector<long>& v)
//   {
//     long prod = v[0];
//     for(int i = 1; i < len(v); ++i)
//     {
//       prod *= v[i];
//     }
//     return prod;
//   }

  
//   long Lcm(const vector<long>& v)
//   {
//     long tmp = v[0];
//     for (int i=1; i < len(v); ++i)
//     {
//       tmp = lcm(tmp, v[i]);
//     }
//     return tmp;
//   }

  
  // // MinPowerInIdealH (for homogeneous ideals): For every homogeneous
  // // component f_i, it determines the minimal exponent m_i, such that
  // // (f_i)^(m_i) is in the ideal J.
  // // Question: Knowing that, can we derive the minimal exponent m, s.t. f^m is in J?
  // vector<long> MinPowerInIdealH(ConstRefRingElem f, const ideal& J)
  // {
  //   vector<long> presumedExponent;
  //   vector<long> homogExponents; // alternative
  //   if (!IsInRadical(f, J))
  //   {
  //     presumedExponent.push_back(-1);
  //   }
  //   else if (!IsHomog(J))
  //   {
  //     presumedExponent.push_back(MinPowerInIdeal(f, J));
  //   }
  //   // f is in the radical of J, J is a radical ideal
  //   else {
  //     vector<RingElem> homogComp = HomogComp(f);	 

  //     for (int i = 0; i < len(homogComp); i++)
  //     {     
  //       homogExponents.push_back(MinPowerInIdeal(homogComp[i], J)); // alternative
  //     }

  //     // first conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the minimal power such that every homogenous component of f is in the Ideal J
  //     // second conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the product of the minimal powers for the homogenous components of f
  //     // such that every homogenous component of f is in the Ideal J
  //     // third conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the lcm of the minimal powers of the homogenous components 
  //     // such that every homogenous component of f is in the Ideal J
  //     // turns out: conjectures are wrong!
  //     // we could not find a simple connection between those powers
  //     presumedExponent.push_back(MaxElem(homogExponents));
  //     presumedExponent.push_back(Prod(homogExponents));
  //     presumedExponent.push_back(Lcm(homogExponents));
  //     cout << "Here the vector of powers for the homogenous components of the polynomial: " << homogExponents << endl; // alternative     
  //   }
  //   return presumedExponent;
  // }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-RadicalMembership.C,v 1.3 2018/07/24 13:46:11 abbott Exp $
// $Log: SparsePolyOps-RadicalMembership.C,v $
// Revision 1.3  2018/07/24 13:46:11  abbott
// Summary: Added CheckForInterrupt to potentially long loops
//
// Revision 1.2  2018/05/18 16:38:52  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.1  2018/04/06 15:14:10  bigatti
// -- renamed RadicalMembership.C
//
// Revision 1.2  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.1  2017/07/19 16:39:02  abbott
// Summary: Added RadicalMembership
//
//
