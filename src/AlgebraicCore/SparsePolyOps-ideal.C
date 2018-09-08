//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2018  John Abbott and Anna M. Bigatti
//             2017  Alice Moallemy (translation CoCoA5-->CoCoALib)

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


// Source code for ideals in SparsePolyRing (functions and member functions)

#include "CoCoA/SparsePolyOps-ideal.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/DenseMatrix.H" // for MultiplicationMat/myDiv
#include "CoCoA/MatrixOps.H" // for LinSolve
#include "CoCoA/MatrixView.H" // for ZeroMat
#include "CoCoA/NumTheory.H" // for CRT
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-MinPoly.H" // for MinPolyDef
#include "CoCoA/SparsePolyOps-MonomialIdeal.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGOperations.H"  // for myIntersect, my Elim..
#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H" // for ideal ctor
#include "CoCoA/apply.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"  // for myGcd
#include "CoCoA/random.H" // for RandomLongStream
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"
//#include "CoCoA/ideal.H"     // for myGcd
//#include "CoCoA/matrix.H" // for OrdMat, myDivMod

#include <algorithm>
//using std::max;     // for MaxExponent, StdDeg
using std::remove;  // for myColon
using std::sort;    // for AreGoodIndetNames, QuotientBasisSorted
//#include <functional>
//using std::not1;    // for AreLPPSqFree
//using std::ptr_fun; // for AreLPPSqFree
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
//#include <iterator>
//using std::back_inserter;
#include <list>
//#include <vector>
using std::vector;

namespace CoCoA
{


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGens() const
  { return myGensValue; }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myTidyGens(const CpuTimeLimit& CheckForTimeOut) const
  { return myGBasis(CheckForTimeOut); }


  const std::vector<RingElem>& GBasis(const ideal& I,
                                      const CpuTimeLimit& timeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasis(I)");
    return TidyGens(I, timeout);
  }


  const std::vector<RingElem>& GBasisByHomog(const ideal& I,
                                              const CpuTimeLimit& timeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasisByHomog(I, timeout)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisByHomog(timeout);
  }


  std::vector<RingElem> GBasisSelfSatCore(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasisSelfSatCore(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisSelfSatCore();
  }


  std::vector<RingElem> GBasisRealSolve(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasisRealSolve(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisRealSolve();
  }


  const std::vector<RingElem>& ReducedGBasis(const ideal& I)
  {
    if (IsCommutative(RingOf(I))) return GBasis(I); // the same (2017)
    CoCoA_ERROR(ERR::NYI, "ReducedGBasis non-commutative");
    return GBasis(I); // just to keep the compiler quiet
  }


  const std::vector<RingElem>& MinGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "MinGens(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myMinGens();
  }


  std::vector<ideal> SparsePolyRingBase::IdealImpl::myPrimaryDecomposition() const
  {
    //if (IhaveMonomialGens()) return myPrimaryDecomposition_MonId();
    if (IhaveSqFreeMonomialGens()) return myPrimaryDecomposition_MonId();
    if (IamZeroDim()) return myPrimaryDecomposition_0dim();
    CoCoA_ERROR(ERR::NYI, "myPrimaryDecomposition() -- general ideal");
    return vector<ideal>(); // just to keep compiler quiet
  }


  ideal LT(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "LT(I)");
    std::vector<RingElem> GB = TidyGens(I);
    std::vector<RingElem> v;
    const SparsePolyRing P = RingOf(I);
    for (long i=0; i<len(GB); ++i)
      v.push_back(monomial(P, LPP(GB[i])));
    return ideal(P, v);
  }


  ideal LF(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "LF(I)");
    const SparsePolyRing P = RingOf(I);
    if ( GradingDim(P)==0 ) CoCoA_ERROR("GradingDim must be non-0", "LF(I)");
    std::vector<RingElem> GB = TidyGens(I);
    std::vector<RingElem> v;
    for (long i=0; i<len(GB); ++i)  v.push_back(LF(GB[i]));
    return ideal(P, v);
  }


  ideal homog(const ideal& I, ConstRefRingElem x)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "homog(I, x)");
    if (AreGensMonomial(I)) return I;
    if (IsZero(I)) return I;
    std::vector<RingElem> HomogIdealGens;
    std::vector<RingElem> v(1,x);
    ComputeHomogenization(HomogIdealGens, gens(I), v);
    return ideal(HomogIdealGens);
  }


  ideal IdealOfGBasis(const ideal& I)
  {
    ideal J(RingOf(I), GBasis(I));
    SetGBasisAsGens(J);
//     if (!uncertain3(IamPrime3Flag)) ...
//     if (!uncertain3(IamMaximal3Flag
    return J;
  }


  ideal IdealOfMinGens(const ideal& I)
  {
    ideal J = I;
    MinGens(J); // so GBasis and such are stored in original ideal
    MakeUnique(J)->myMinimalize();
    return J;
  }


  // Anna: must be friend for ourGetPtr
  bool HasGBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "HasGBasis(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveGBasis();
  }
  

  // Anna: must be friend for ourGetPtr
  bool AreGensMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "AreGensMonomial(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveMonomialGens();
  }
  

  // Anna: must be friend for ourGetPtr
  bool AreGensSqFreeMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "AreGensSqFreeMonomial(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveSqFreeMonomialGens();
  }
  

  // Anna: must be friend for ourGetPtr
  void SetGBasisAsGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "GBasisAsGens(I)");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    ptrI->mySetGBasisAsGens();
  }
  

//   namespace  //anonymous  for QuotientBasis
//   {

//     bool IsDivisible(ConstRefPPMonoidElem pp, const std::list<PPMonoidElem>& ByL)
//     {
//       for (std::list<PPMonoidElem>::const_iterator i=ByL.begin(); i != ByL.end(); ++i)
//         if (IsDivisible(pp, *i)) return true;
//       return false;
//     }
    

// //     // ??? copied from ex-QuotientBasis.C
// //     bool IsDivisible(ConstRefPPMonoidElem pp, const std::vector<PPMonoidElem>& ByL)
// //     {
// //       const long n = len(ByL);
// //       for (long i=0; i < n; ++i)
// //         if ( IsDivisible(pp, ByL[i]) ) return true;
// //       return false;
// //     }


//     // Apparently STL bind2nd cannot have (const?) reference arguments.
//     // Apparently BOOST bind would work: what about C++-11?
//     // This is a hack for STL
//     class IsDivHack
//     {
//     public:
//       IsDivHack (const std::list<PPMonoidElem>& L): myL(L) {}
//       bool operator () (ConstRefPPMonoidElem pp) {return IsDivisible(pp, myL);}
      
//     private:
//       const std::list<PPMonoidElem> & myL;
//     };


//     void QuotientBasisRec(std::vector<PPMonoidElem>& ans, 
//                           const std::list<PPMonoidElem>& L, 
//                           ConstRefPPMonoidElem prefix, 
//                           long idx)
//     {
//       PPMonoid PPM = owner((L.front()));
//       const PPMonoidElem& X = indets(PPM)[idx];
//       PPMonoidElem prefixXd(prefix);  // prefix * x[idx]^d
//       int MaxDeg = 0;
//       if (idx == NumIndets(PPM)-1)
//       {
//         MaxDeg = exponent((L.front()), idx);
//         for (int d=0; d < MaxDeg;  ++d, prefixXd *= X)  ans.push_back(prefixXd);
//         return;
//       }
//       for (std::list<PPMonoidElem>::const_iterator it=L.begin(); it != L.end() ; ++it)
//         if (exponent((*it),idx) > MaxDeg)  MaxDeg = exponent((*it),idx);
//       std::list<PPMonoidElem> CutOff, tmp;
//       PPMonoidElem Xd(PPM);  // x[idx]^d
//       for (int d=0; d < MaxDeg; ++d, prefixXd *= X, Xd *= X)
//       {
//         for (std::list<PPMonoidElem>::const_iterator it=L.begin(); it != L.end(); ++it)
//           if (exponent(*it,idx) == d)  tmp.push_back(((*it)/Xd));
//         CutOff.remove_if(IsDivHack(tmp));
//         CutOff.splice(CutOff.end(), tmp);
//         QuotientBasisRec(ans, CutOff, prefixXd, idx+1);
//       }
//     }
//   } // anonymous namespace


//   std::vector<PPMonoidElem> QuotientBasis(const ideal& I)
//   {
//     const char* const fn = "QuotientBasis";
//     if (!IsSparsePolyRing(RingOf(I))) CoCoA_ERROR(ERR::NotSparsePolyRing, fn);
//     if (!IsZeroDim(I)) CoCoA_ERROR("ideal must be 0-dimensional", fn);
//     vector<RingElem> GB = GBasis(I);
//     std::list<PPMonoidElem> LeadingPPs;
//     vector<PPMonoidElem> ans;
//     for (long i=0; i < len(GB); ++i)  LeadingPPs.push_back(LPP(GB[i]));
//     QuotientBasisRec(ans, LeadingPPs, PPMonoidElem(PPM(RingOf(I))), 0);
//     return ans;
//   }


//   std::vector<PPMonoidElem> QuotientBasisSorted(const ideal& I)
//   {
//     const char* const fn = "QuotientBasis";
//     if (!IsSparsePolyRing(RingOf(I))) CoCoA_ERROR(ERR::NotSparsePolyRing, fn);
//     if (!IsZeroDim(I)) CoCoA_ERROR("ideal must be 0-dimensional", fn);
//     vector<PPMonoidElem> QB = QuotientBasis(I);
//     std::sort(QB.begin(), QB.end());
//     return QB;
//   }





  //-- IdealImpl ----------------------------------------

  ideal SparsePolyRingBase::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(SparsePolyRing(this), gens)); //??? ugly ???
  }


  SparsePolyRingBase::IdealImpl::IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens):
      myP(P),
      myGensValue(gens),
      myInvBasisContainerPtr(new Involutive::UniversalInvolutiveBasisContainer(gens)),
      IhaveGBasisValue(false)
  {
    // IhaveMonomialGens3Value = uncertain3; // default for bool3
    // IhaveSqFreeMonomial3Gens = uncertain3; // default for bool3
    if (!IsField(CoeffRing(P)))
      CoCoA_ERROR("ERR:NYI ideal of polynomials with coeffs not in a field", "ideal(SparsePolyRing, gens)");//???
  }


  IdealBase* SparsePolyRingBase::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


  const SparsePolyRing& SparsePolyRingBase::IdealImpl::myRing() const
  {
    return myP;
  }


  bool SparsePolyRingBase::IdealImpl::IamZero() const
  {
    for (long i=0; i<len(myGens()); ++i)
      if (!IsZero(myGens()[i])) return false;
    return true;
  }


  
  namespace // anonymous --------------------------------
  {
    bool AreLPPSqFree(const std::vector<RingElem>& v)
    {
      const long n = len(v);
      for (long i=0; i < n; ++i)
        if (!IsSqFree(LPP(v[i]))) return false;
      return true;
//   We *DO NOT USE* STL algorithm because std::ptr_fun does not work if the fn has formal params which are of reference type
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(CoCoA::IsSqFreeLPP)))
// 	//                     not1(ptr_fun(static_cast<bool(*)(ConstRefRingElem)>(CoCoA::IsRadLPP))))
//         == v.end();
    }

    } // anonymous end ----------------------------------


  bool SparsePolyRingBase::IdealImpl::IhaveGBasis() const
  { return IhaveGBasisValue; }


  bool SparsePolyRingBase::IdealImpl::IhaveMonomialGens() const
  {
    if (IsUncertain3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = AreMonomials(myGensValue);
    return IsTrue3(IhaveMonomialGens3Value);
  }


  bool SparsePolyRingBase::IdealImpl::IhaveSqFreeMonomialGens() const
  {
    if (IsUncertain3(IhaveSqFreeMonomialGens3Value))
    {
      if (!IhaveMonomialGens()) IhaveSqFreeMonomialGens3Value = false3;
      else IhaveSqFreeMonomialGens3Value = AreLPPSqFree(myGensValue);
    }
    return IsTrue3(IhaveSqFreeMonomialGens3Value);
  }


  void SparsePolyRingBase::IdealImpl::mySetGBasisAsGens() const
  {
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
  }


  void SparsePolyRingBase::IdealImpl::myReset() const
  {
    IamMaximal3Flag = uncertain3;
    IamPrimary3Flag = uncertain3;
    IamPrime3Flag = uncertain3;
    IamRadical3Flag = uncertain3;
    IhaveSqFreeMonomialGens3Value = uncertain3;
    IhaveMonomialGens3Value = uncertain3;
    IhaveGBasisValue = false;
    myGBasisValue.clear();
    myMinGensValue.clear();
  }


  bool SparsePolyRingBase::IdealImpl::myTestIsMaximal() const
  {
    if (IamZero()) { return myAssignMaximalFlag(false); }// ideal <0> in K[X]
    if (!IsField(CoeffRing(myP)))
      CoCoA_ERROR(ERR::ExpectedCoeffsInField, "myTestIsMaximal()");//???
    if (NumIndets(myP) == 1)  // Simple case: univariate poly ring 
      return myAssignMaximalFlag(IsIrred(myGBasis(NoCpuTimeLimit())[0]));
    if (!IsZeroDim(ideal(const_cast<IdealImpl*>(this))))  return false;
    // Harder case: poly ring has at least 2 indets.
    return myAssignMaximalFlag(myTestIsMaximal_0dim());
  }


  bool SparsePolyRingBase::IdealImpl::myTestIsPrimary() const
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_ERROR(ERR::ExpectedCoeffsInField, "myTestIsPrimary()");//???
    if (NumIndets(myP) == 1)  // Simple case: univariate poly ring 
    {
      const factorization<RingElem> F = factor(myGBasis(NoCpuTimeLimit())[0]); // PID
      if (len(F.myFactors()) != 1) return myAssignPrimeFlag(false);
      if (F.myMultiplicities()[0] == 1) myAssignMaximalFlag(true);
      return true;
    }
    if (!IsZeroDim(ideal(const_cast<IdealImpl*>(this))))
      CoCoA_ERROR("Not Yet Implemented for general ideals", "myTestIsPrimary");
    return myTestIsPrimary_0dim();
  }


  bool SparsePolyRingBase::IdealImpl::myTestIsPrime() const
  {
    if (NumIndets(myP) == 1 && IsField(CoeffRing(myP)))
      return myTestIsMaximal();
    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myTestIsPrime()");//???
    return false; // just to keep the compiler quiet
  }


  bool SparsePolyRingBase::IdealImpl::myTestIsRadical() const
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_ERROR(ERR::ExpectedCoeffsInField, "myTestIsRadical()");//???
    if (IhaveMonomialGens()) return myTestIsRadical_MonId();
    if (NumIndets(myP) == 1)  // Simple case: univariate poly ring 
    {
      const factorization<RingElem> F = factor(myGBasis(NoCpuTimeLimit())[0]); // PID
      if (len(F.myFactors()) != 1) myAssignPrimeFlagPID(false);
      if (F.myMultiplicities()[0] == 1) myAssignMaximalFlagPID(true);
      return true;
    }
    if (!IsZeroDim(ideal(const_cast<IdealImpl*>(this))))
      CoCoA_ERROR("Not Yet Implemented for general ideals", "myTestIsRadical");
    return myTestIsRadical_0dim();
  }


  void SparsePolyRingBase::IdealImpl::myReduceMod(RingElemRawPtr rawf) const
  {
    //??? very basic default implementation
    RingElem tmp = NR(RingElemAlias(myP, rawf), myGBasis(NoCpuTimeLimit()));
    myP->mySwap(rawf, raw(tmp));
  }


  bool SparsePolyRingBase::IdealImpl::IhaveElem(RingElemConstRawPtr rawf) const
  {
    RingElem g = RingElemAlias(myP, rawf);
    myReduceMod(raw(g));
    return IsZero(g);
  }


  bool SparsePolyRingBase::IdealImpl::IamZeroDim() const
  {
    const vector<RingElem>& GB = myTidyGens(NoCpuTimeLimit());
    const long GBlen = len(GB); // MUST BE A REDUCED GBASIS !!!
    long NumIndetPowers = 0;
    for (long i=0; i < GBlen; ++i)
      if (IsIndetPosPower(LPP(GB[i])))  ++NumIndetPowers;
    return (NumIndetPowers == NumIndets(myRing()));
  }
  

  const SparsePolyRingBase::IdealImpl* SparsePolyRingBase::IdealImpl::ourGetPtr(const ideal& I)
  {
    return dynamic_cast<const SparsePolyRingBase::IdealImpl*>(I.myIdealPtr());
  }


  void SparsePolyRingBase::IdealImpl::myAdd(const ideal& Jin)
  {
    const IdealImpl* const J = ourGetPtr(Jin);
    if (IsZero(Jin)) return;
    myGensValue.insert(myGensValue.end(), gens(Jin).begin(), gens(Jin).end());
    bool3 IhaveMG_old = IhaveMonomialGens3Value;
    bool3 IhaveSFMG_old = IhaveSqFreeMonomialGens3Value;
    myReset();
    // we can recover some info about monomial gens:
    IhaveMonomialGens3Value = IhaveMG_old;
    IhaveSqFreeMonomialGens3Value = IhaveSFMG_old;
    if (IsTrue3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = J->IhaveMonomialGens3Value;
    if (IsTrue3(IhaveSqFreeMonomialGens3Value))
      IhaveSqFreeMonomialGens3Value = J->IhaveSqFreeMonomialGens3Value;
  }


  void SparsePolyRingBase::IdealImpl::myMul(const ideal& Jin)
  {
    if (IhaveMonomialGens() && AreGensMonomial(Jin))
    {
      myMul_MonId(Jin);
      return;
    }
    vector<RingElem> tmpV;
    const SparsePolyRingBase::IdealImpl* const J = ourGetPtr(Jin);
    for (vector<RingElem>::const_iterator itI=myGensValue.begin(); itI!=myGensValue.end(); ++itI)
      for (vector<RingElem>::const_iterator itJ=J->myGensValue.begin(); itJ!=J->myGensValue.end(); ++itJ)
        tmpV.push_back((*itI)*(*itJ));
    swap(tmpV, myGensValue);  // ANNA does this make copies???  2010-02-03
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::myIntersect(const ideal& J)
  {
    if (IamZero()) return;
    CoCoA_ASSERT(!IsZero(J));
    if (IhaveMonomialGens() && AreGensMonomial(J))
    {
      myIntersect_MonId(J);
      return;
    }
    ComputeIntersection(myGensValue, myGensValue, gens(J));
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::myColon(const ideal& J)
  {
    if (IsZero(J))
      myGensValue = vector<RingElem>(1, one(myRing()));
    else
    {
      if (IhaveMonomialGens() && AreGensMonomial(J))
      {
        myColon_MonId(J);
        return;
      }
      const RingElem Z(zero(myRing()));
      myGensValue.erase(remove(myGensValue.begin(), myGensValue.end(),Z),
                        myGensValue.end());
      ComputeColon(myGensValue, myGensValue, gens(J));
    }
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::mySaturate(const ideal& J)
  {
    ComputeSaturation(myGensValue,myGensValue, gens(J));
    myReset();
  }


  void SparsePolyRingBase::IdealImpl::myMinimalize()
  {
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_ERROR("Input is not homogeneous", "myMinimalize");
    myGBasis(NoCpuTimeLimit()); // this sets GBasis and MinGens
    myGensValue = myMinGens(); // if monomial ideal min gens are in GBasis
    if (IsFalse3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = uncertain3;
    if (IsFalse3(IhaveSqFreeMonomialGens3Value))
      IhaveSqFreeMonomialGens3Value = uncertain3;
  }


  void SparsePolyRingBase::IdealImpl::myElim(const std::vector<RingElem>& ElimIndets)
  {
    //    if (IhaveGBasisValue) and elim ordering ...
    if (IhaveMonomialGens())
    {
      myElim_MonId(ElimIndets);
      return;
    }
    const RingElem Z(zero(myRing()));
    myGensValue.erase(remove(myGensValue.begin(), myGensValue.end(),Z),
                      myGensValue.end());
    PPMonoidElem ElimIndetsProd(myRing()->myPPM());
    const long n = len(ElimIndets);
    for (long i=0 ; i<n ; ++i)
    {
      if (!IsIndet(ElimIndets[i]))
        CoCoA_ERROR(ERR::NotIndet, "myElim");
      ElimIndetsProd *= LPP(ElimIndets[i]);
    }
    ComputeElim(myGensValue, myGensValue, ElimIndetsProd);
    myReset(); // can we recover some info?
  }


  namespace {  // anonymous ------------------------------
    RingElem CoeffOfTermSparse(ConstRefRingElem f, ConstRefPPMonoidElem pp)
    {
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        if (PP(itf) == pp) return coeff(itf);
      return zero(CoeffRing(owner(f)));
    }

    matrix MultiplicationMat(ConstRefRingElem f, const ideal& I)
    {
      std::vector<PPMonoidElem> QB = QuotientBasis(I);
      SparsePolyRing P = owner(f);
      matrix Mf = NewDenseMat(CoeffRing(P), len(QB), len(QB));
      for (long j=0; j<len(QB); ++j)
      {
        RingElem tmpf = f;
        P->myMulByPP(raw(tmpf), raw(QB[j]));
        tmpf = NF(tmpf, I);
        for (long i=0; i<len(QB); ++i)
          SetEntry(Mf,i,j, CoeffOfTermSparse(tmpf,QB[i]));
      }
      return Mf;
    }

  } // anonymous end -------------------------------------------


  bool SparsePolyRingBase::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
// check for IsZeroDivisor(den) is done by operator/
    const SparsePolyRing P = myRing();
    //    if (IsField(CoeffRing(P)) && P->myIsConstant(rawden))
    if (P->myIsConstant(rawden))
    {
      RingElem GDivF = RingElemAlias(P, rawnum);
      P->myDivByCoeff(raw(GDivF), raw(P->myLC(rawden))); // exc safe?
      P->mySwap(rawlhs, raw(GDivF));
      return true;
    }
//  if (IsZeroDim(ideal(this))) // !!! DOES NOT COMPILE BECAUSE const!!!
    ideal I = ideal(const_cast<IdealImpl*>(this)); // Tappullus Horribilis
    if (IsZeroDim(I))
    {
      std::vector<PPMonoidElem> QB = QuotientBasis(I);
      long LenQB = len(QB);
      RingElem G = RingElemAlias(P,rawnum);
      RingElem F = RingElemAlias(P,rawden);
      matrix coeffsG = NewDenseMat(ZeroMat(CoeffRing(P), len(QB), 1));
      for (long i=0; i<LenQB; ++i)
        SetEntry(coeffsG,i,0, CoeffOfTermSparse(G, QB[i]));
      matrix coeffsGDivF = LinSolve(MultiplicationMat(F,I), coeffsG);
      RingElem GDivF(P);
      for (long i=0; i<LenQB; ++i)
        GDivF += monomial(P, coeffsGDivF(i,0), QB[i]);
      P->mySwap(rawlhs, raw(GDivF));
      return true;
    }

// The following should be a general solution...
///    auto tmp = GenRepr(rawnum, ideal(rawden)+this);
///    return tmp[0];

    CoCoA_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myDivMod");
    return false; // just to keep the compiler quiet!!
  }




  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasis(const CpuTimeLimit& CheckForTimeOut) const
  {
    if (IhaveMonomialGens()) return myGBasis_MonId();
    if (IhaveGBasis()) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    vector<RingElem> MinGens;
    ComputeGBasis(myGBasisValue, MinGens, myGensValue, CheckForTimeOut);
    if (!MinGens.empty())  myMinGensValue = MinGens;
    IhaveGBasisValue = true;
    return myGBasisValue;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasisByHomog(const CpuTimeLimit& CheckForTimeOut) const
  {
    if (IhaveMonomialGens()) return myGBasis_MonId();
    if (IhaveGBasis()) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    // check GradingDim
    // make correct ordering
    SparsePolyRing  P = myRing();
    if (!IsStdDegRevLex(ordering(PPM(P))))
      CoCoA_ERROR(ERR::NYI, "GBasisByHomog not drl");
    long n = NumIndets(P);
    SparsePolyRing Ph = NewPolyRing(CoeffRing(P), NewSymbols(n+1));
    const std::vector<RingElem>& X = indets(Ph);
    RingHom phi =PolyAlgebraHom(P,Ph,vector<RingElem>(X.begin(),X.begin()+n));
    std::vector<RingElem> v = indets(P);  v.push_back(one(P));
    RingHom dehom = PolyAlgebraHom(Ph,P,v);
    const RingElem& h = indet(Ph, n);
    const std::vector<RingElem>& g = myGens();
    std::vector<RingElem> gh;
    for (long i=0; i<len(g); ++i)  gh.push_back(homog(phi(g[i]), h));
    std::vector<RingElem> GB = apply(dehom, GBasis(ideal(gh)));
    // trick: recompute to remove redundant elements in GB (interreduction)
    gh.clear();
    for (long i=0; i<len(GB); ++i)  gh.push_back(homog(phi(GB[i]), h));
    GB = apply(dehom, GBasis(ideal(gh), CheckForTimeOut));
    //
    swap(myGBasisValue, GB);
    IhaveGBasisValue = true;
    return myGBasisValue;
  }


  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisSelfSatCore() const
  {
    if (IhaveMonomialGens()) return myGBasis_MonId();
    if (IhaveGBasis()) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    vector<RingElem> GBSelfSat;
    ComputeGBasisSelfSatCore(GBSelfSat, myGensValue);
    return GBSelfSat;
  }


  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisRealSolve() const
  {
    if (IamZero()) return myGBasisValue;
    vector<RingElem> GBRealSolve;
    ComputeGBasisRealSolve(GBRealSolve, myGensValue);
    return GBRealSolve;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myMinGens() const
  {
    if (IhaveMonomialGens()) return myGBasis_MonId(); // interreduced
    if (!myMinGensValue.empty()) return myMinGensValue;
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_ERROR("Input is not homogeneous", "myMinGens");
    myGBasis(NoCpuTimeLimit());
    return myMinGensValue;
  }


  bool IsZeroDim(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsZeroDim(I)");
    if (IsZero(I) || IsOne(I)) return false;
    // Now we know I is non-trivial.
//     const SparsePolyRing P = RingOf(I);
//     const vector<RingElem>& GB = TidyGens(I);
//     const long GBlen = len(GB); // MUST BE A REDUCED GBASIS !!!
//     long NumIndetPowers = 0;
//     for (long i=0; i < GBlen; ++i)
//       if (IsIndetPosPower(LPP(GB[i])))
//         ++NumIndetPowers;
//     return (NumIndetPowers == NumIndets(P));
    return SparsePolyRingBase::IdealImpl::ourGetPtr(I)->IamZeroDim();
  }


  bool IsHomog(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsHomog(ideal)");
    if (GradingDim(RingOf(I))==0)
      CoCoA_ERROR(ERR::ZeroGradingDim, "IsHomog(ideal)");
    if (IsZero(I) || IsOne(I)) return true;
    // Now we know I is non-trivial.
    const SparsePolyRing P = RingOf(I);
    const vector<RingElem>& GB = TidyGens(I);
    const long GBlen = len(GB); // MUST BE A REDUCED GBASIS !!!
    for (long i=0; i < GBlen; ++i)
      if (!IsHomog(GB[i]))  return false;
    return true;
  }


  RingElem DenSigma(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "DenSigma(I)");
    if (!IsFractionField(CoeffRing(RingOf(I))))
      CoCoA_ERROR("Coeff ring must be a fraction field of a GCDDomain", "DenSigma(I)");
    return CommonDenom(ReducedGBasis(I));
  }

  // template???
  bool IsSigmaGoodPrime(const BigInt& p, const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsSigmaGoodPrime(p,I)");
    if (!IsQQ(CoeffRing(RingOf(I))))
      CoCoA_ERROR("Coeff ring must be RingQQ", "IsSigmaGoodPrime(p,I)");
    if (!IsPrime(p))
      CoCoA_ERROR("First argument must be a prime number", "IsSigmaGoodPrime(p,I)");
    return !IsDivisible(CommonDenom(ReducedGBasis(I)), p);
  }

  // template???
  bool IsSigmaGoodPrime(const long p, const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_ERROR(ERR::NotSparsePolyRing, "IsSigmaGoodPrime(p,I)");
    if (!IsQQ(CoeffRing(RingOf(I))))
      CoCoA_ERROR("Coeff ring must be RingQQ", "IsSigmaGoodPrime(p,I)");
    if (!IsPrime(p))
      CoCoA_ERROR("First argument must be a prime number", "IsSigmaGoodPrime(p,I)");
    return !IsDivisible(CommonDenom(ReducedGBasis(I)), p);
  }


  //  namespace { // anonymous
    ideal radical_0dimDRL(const ideal& I)
    {
      const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
      return ptrI->myRadical_0dimDRL(); // behaves differently from other memb fns
    }
  //  }
  
  
// will be radical(I)
  ideal radical_0dim(const ideal& I)
  {
    //    const SparsePolyRingBase::IdealImpl* const ptrI = 
    //      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    if (IsTrue3(IsRadical3(I))) return I;
    PolyRing P = RingOf(I);
    //    if (HasStdDegRevLex(P))  return radical_0dimDRL(I);
    if (HasGBasis(I)) return radical_0dimDRL(I);
    PolyRing P_drl = NewPolyRing(CoeffRing(P), NumIndets(P));
    RingHom phi = PolyAlgebraHom(P, P_drl, indets(P_drl));
    RingHom psi = PolyAlgebraHom(P_drl, P, indets(P));
    ideal RadI = radical_0dimDRL(ideal(apply(phi, gens(I))));
    return ideal(apply(psi, gens(RadI)));
  } // radical_0dim(I)
  





} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-ideal.C,v 1.21 2018/08/06 09:38:29 bigatti Exp $
// $Log: SparsePolyOps-ideal.C,v $
// Revision 1.21  2018/08/06 09:38:29  bigatti
// -- renamed GBasisViaHomog --> GBasisByHomog
//
// Revision 1.20  2018/08/06 08:57:48  bigatti
// -- added timeout for GBasisByHomog
// -- now using GBasisByHomog in IsPrimary and radical
//
// Revision 1.19  2018/08/05 16:31:25  bigatti
// -- added GBasisByHomog
//
// Revision 1.18  2018/06/27 12:15:19  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.17  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.16  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.15  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.14  2018/05/17 15:49:00  bigatti
// -- sorted #includes
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.13  2018/04/18 15:39:34  abbott
// Summary: Added missing return in NYI fn.
//
// Revision 1.12  2018/04/16 21:49:26  bigatti
// -- added IamZeroDim
// -- added myPrimaryDecomposition_0dim
//
// Revision 1.11  2018/04/10 14:59:32  bigatti
// -- now member function for PrimaryDecomposition
//
// Revision 1.10  2018/04/10 14:51:43  bigatti
// -- added virtual myPrimaryDecomposition (with default implementation)
//
// Revision 1.9  2018/04/10 14:20:47  bigatti
// -- fixed includes
//
// Revision 1.8  2018/04/09 16:25:53  bigatti
// -- functions for zer-dim ideals moved into SparsePolyOps-IdealZeroDim
//
// Revision 1.7  2018/04/06 17:28:22  bigatti
// -- fixed includes
//
// Revision 1.6  2018/04/04 12:36:51  bigatti
// -- radical: returning I itself, in IsRadical3 = true
//
// Revision 1.5  2018/03/29 13:09:04  bigatti
// -- improved isprimary (with gbasis timeout)
//
// Revision 1.4  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.3  2018/03/20 11:48:12  bigatti
// -- new code for 0dim ideals (radical, primary...)
//
// Revision 1.2  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.1  2018/03/12 14:38:41  bigatti
// -- renamed SparsePoly-ideal/involutive into SparsePolyOps-ideal/involutive
//
// Revision 1.3  2018/03/09 14:54:01  bigatti
// -- removed useless includes
//
// Revision 1.2  2018/02/27 17:26:07  bigatti
// -- split off involutive part
//
// Revision 1.1  2018/02/27 17:12:37  bigatti
// -- Renamed SparsePolyRing_ideal.C --> SparsePoly-ideal.C
//
