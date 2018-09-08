//   Copyright (c)  1999,2009-2011  John Abbott

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


#include "CoCoA/NumTheory.H"
#include "CoCoA/NumTheory-prime.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
#include "CoCoA/config.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ToString.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H"

#include <algorithm>
using std::min;
using std::max;
using std::swap;
#include <cmath>
// myThreshold uses floor,pow,sqrt
#include <cstdlib>
using std::ldiv;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {

    // Function to compute base^exp mod modulus; assumes 0 <= base < modulus.
    // Unclever for exp >= modulus.  modulus need not be prime, but (modulus-1)^2
    // must fit inside an unsigned long.
    // This is an iterative implementation of binary powering; seems to be
    // usefully faster than the more obvious recursive implementation.
    // If modulus=1, always produces 0.  Produces 1 when asked to compute 0^0.
    unsigned long PowerModSmallModulus(unsigned long base, unsigned long exp, unsigned long modulus)
    {
      CoCoA_ASSERT(modulus != 0);
      CoCoA_ASSERT(modulus > 1);
      CoCoA_ASSERT(modulus-1 <= MaxSquarableInteger<unsigned long>());

      if (exp == 0 || base == 1) return 1;
      if (exp == 1 || base == 0) return base;
    
      unsigned long ans = 1;
      while (exp > 1)
      {
        if (exp&1) ans = (ans*base)%modulus;
        exp /= 2;
        base = (base*base)%modulus;
      }
      return (ans*base)%modulus;
    }


    // Common sanity check for the modulus
    void CheckModulus(const MachineInt& modulus, const char* const FnName)
    {
      if (IsNegative(modulus) || !IsSignedLong(modulus) || AsSignedLong(modulus) < 2)
        CoCoA_ERROR(ERR::BadModulus, FnName);
    }

    // Common sanity check for the modulus
    void CheckModulus(const BigInt& modulus, const char* const FnName)
    {
      if (modulus < 2)
        CoCoA_ERROR(ERR::BadModulus, FnName);
    }

  } // end of anonymous namespace


  // Compute (base^exp)%modulus  with exp >= 0  and modulus >= 2
  // The case 0^0 yields 1.
  long PowerMod(const MachineInt& base, const MachineInt& exponent, const MachineInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    if (IsNegative(exponent))
    {
      const long inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, negate(exponent), modulus);
    }
    const unsigned long e = AsUnsignedLong(exponent);
    if (e == 0) return 1; // Note that 0^0 gives 1
    const unsigned long m = AsUnsignedLong(modulus);
    unsigned long b = uabs(base)%m;
    if (b == 0) return 0;
    if (IsNegative(base) && (e%2 == 1))
      b = m-b;
    if (m <= MaxSquarableInteger<unsigned long>())
      return PowerModSmallModulus(b,e,m); // Call the fn in the anon namespace above.
    // Square of modulus doesn't fit into unsigned long, so compute with BigInts and convert answer back to unsigned long.
    return ConvertTo<long>(PowerMod(b,e,BigInt(m)));
  }

  long PowerMod(const MachineInt& base, const BigInt& exponent, const MachineInt& modulus)
  {
    // BUG/SLUG: horribly inefficient!!
    // What is the best way to implement this fn?????  Any volunteers?
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return ConvertTo<long>(PowerMod(BigInt(base), exponent, BigInt(modulus)));
  }

  long PowerMod(const BigInt& base, const MachineInt& exponent, const MachineInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return PowerMod(base%modulus, exponent, modulus);
  }

  long PowerMod(const BigInt& base, const BigInt& exponent, const MachineInt& modulus)
  {
    // BUG/SLUG: horribly inefficient!!
    // What is the best way to implement this fn?????  Any volunteers?
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return ConvertTo<long>(PowerMod(base, exponent, BigInt(modulus)));
  }

  BigInt PowerMod(const MachineInt& base, const MachineInt& exponent, const BigInt& modulus)
  {
    return PowerMod(BigInt(base), BigInt(exponent), modulus);
  }

  BigInt PowerMod(const BigInt& base, const MachineInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(BigInt, MachineInt, BigInt)");
    if (IsNegative(exponent))
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, negate(exponent), modulus);
    }
    BigInt ans;
    // ASSUME: line below gives 1 for 0^0
    mpz_powm_ui(mpzref(ans), mpzref(base), AsUnsignedLong(exponent), mpzref(modulus));
    return ans;
  }


  BigInt PowerMod(const MachineInt& base, const BigInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(MachineInt, BigInt, BigInt)");
    if (exponent < 0)
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, -exponent, modulus);
    }
    return PowerMod(BigInt(base), exponent, modulus);
  }


  BigInt PowerMod(const BigInt& base, const BigInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(BigInt, BigInt, BigInt)");
    if (IsZero(exponent)) return BigInt(1); // x^0 --> 1 even for x==0
    if (exponent < 0)
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, -exponent, modulus);
    }
    BigInt ans;
    mpz_powm(mpzref(ans), mpzref(base), mpzref(exponent), mpzref(modulus));
    return ans;
  }



  // Non-negative gcd of machine integers.
  // Emphasis is on simplicity rather than utmost speed.
  long gcd(const MachineInt& a, const MachineInt& b)
  {
    unsigned long A = uabs(a);
    unsigned long B = uabs(b);
    // Dispose of the trivial cases first.
    if (A == 0) return B;
    if (B == 0) return A;
    if (A == 1 || B == 1) return 1;
    if (A == B) return A;

    // General case.
    if (A < B) std::swap(A, B);
    while (B != 0)
    {
      A %= B;
      std::swap(A, B);
    }
    return NumericCast<long>(A);
  }

  const BigInt gcd(const BigInt& A, const MachineInt& b)
  { if (IsZero(b)) return abs(A); else return BigInt(gcd(A%uabs(b), b)); }

  const BigInt gcd(const MachineInt& a, const BigInt& B)
  { if (IsZero(a)) return abs(B); else return BigInt(gcd(a, B%uabs(a))); }

  const BigInt gcd(const BigInt& A, const BigInt& B)
  {
    BigInt ans;
    mpz_gcd(mpzref(ans), mpzref(A), mpzref(B));
    return ans;
  }


  // -------------------------------------------------------
  // IsCoprime (trivial defns)
  
  bool IsCoprime(const MachineInt& a, const MachineInt& b)
  { return (gcd(a,b) == 1); }
  
  bool IsCoprime(const BigInt& A,     const MachineInt& b)
  { return IsOne(gcd(A,b)); }

  bool IsCoprime(const MachineInt& a, const BigInt& B)
  { return IsOne(gcd(a,B)); }
  
  bool IsCoprime(const BigInt& A,     const BigInt& B)
  { return IsOne(gcd(A,B)); }


  namespace // anonymous
  {
    
    // The fn below is OK -- it is not obvious, but the longs cannot overflow: for
    // a proof see Theorem 4.3 & the comment after it in Shoup's book
    // "A Computational Introduction to Number Theory and Algebra".

    // InvModNoArgCheck is faster than this; perhaps mimic InvModNoArgCheck, and then compute directly
    // the other cofactor at the end??? see redmine!
    long ExtendedEuclideanAlg(long& CofacA, long& CofacB, unsigned long a, unsigned long b)
    {
      if (a < b) return ExtendedEuclideanAlg(CofacB, CofacA, b, a);

      // Now we are sure that a >= b.
      long m00 = 1;
      long m01 = 0;
      long m10 = 0;
      long m11 = 1;

      while (true)
      {
        if (b==0) { CofacA = m00; CofacB = m01; return NumericCast<long>(a); }
        long q = a/b; // can overflow (harmlessly) only if b == 1 upon 1st iteration
        a -= q*b;
        m00 -= q*m10;
        m01 -= q*m11;

        if (a == 0) { CofacA = m10; CofacB = m11; return NumericCast<long>(b); }
        q = b/a;
        b -= q*a;
        m10 -= q*m00;
        m11 -= q*m01;
      }
    }

  } // end of namespace anonymous


  long ExtGcd(long& CofacA, long& CofacB, const MachineInt& a, const MachineInt& b)
  {
    if (IsZero(a) && IsZero(b)) CoCoA_ERROR(ERR::NotNonZero, "ExtGcd (machine int)");
    const long g = ExtendedEuclideanAlg(CofacA, CofacB, uabs(a), uabs(b));
    if (IsNegative(a)) CofacA = -CofacA;
    if (IsNegative(b)) CofacB = -CofacB;
    return g;
  }

  BigInt ExtGcd(BigInt& CofacA, BigInt& CofacB, const BigInt& A, const BigInt& B)
  {
    if (IsZero(A) && IsZero(B)) CoCoA_ERROR(ERR::NotNonZero, "ExtGcd (BigInt)");
    BigInt ans;
    mpz_gcdext(mpzref(ans), mpzref(CofacA), mpzref(CofacB), mpzref(A), mpzref(B));
    // GMP guarantees only that abs(CofacA) < abs(B), but we want 2*abs(CofacA) <= abs(B)/ans
    // so the next few lines check, and if necessary tweak, the cofactors.
    // The code below is needed only for GMP-4.3.1; GMP-4.2.4 gave correctly reduced answers.
    if (A == 0 || B == 0 || abs(A) == abs(B)) return ans;
    const BigInt CofacAModulus = abs(B)/ans;
    const BigInt CofacBModulus = abs(A)/ans;
    const BigInt q = RoundDiv(abs(CofacA), CofacAModulus);
    if (q != 0)
    {
      if (CofacA < 0)
      { CofacA += q*CofacAModulus; if (sign(A) == sign(B)) CofacB -= q*CofacBModulus; else CofacB += q*CofacBModulus; }
      else
      { CofacA -= q*CofacAModulus; if (sign(A) == sign(B)) CofacB += q*CofacBModulus; else CofacB -= q*CofacBModulus; }
    }
    if (2*abs(CofacB) <= CofacBModulus) return ans; // NB already have 2*abs(CofacA) <= CofacAModulus
    if (CofacB < 0)
    { CofacB += CofacBModulus; if (sign(A) == sign(B)) CofacA -= CofacAModulus; else CofacA += CofacAModulus; }
    else
    { CofacB -= CofacBModulus; if (sign(A) == sign(B)) CofacA += CofacAModulus; else CofacA -= CofacAModulus; }
    return ans;
  }


  namespace // anonymous for file local fn
  {
    // On some platforms this "unsigned int" version is significantly faster
    // than the standard "unsigned long" version.  This version is "hidden
    // inside" the unsigned long version -- it calls this version if it can.

    // JAA believes that this routine always returns a reduced value,
    // but does not have a proof of this currently.
    unsigned int InvModNoArgCheck_int(unsigned int r, unsigned int m, const ErrorActionIndicator ErrorAction = ThrowOnError)
    {
//      CoCoA_ASSERT(r >= 0);
      CoCoA_ASSERT(r < m);
      CoCoA_ASSERT(m >= 2);
      const unsigned int p = m;
      unsigned int cr = 1;
      unsigned int cm = 0; // this is minus what you think it is!
      while (r != 0)
      {
        {
          const unsigned int q = m/r;
          m -= q*r;
          cm += q*cr;
        }
        if (m == 0) break;
        {
          const unsigned int q = r/m;
          r -= q*m;
          cr += q*cm;
        }
      }
      if (r+m != 1)
      {
        if (ErrorAction == ThrowOnError) CoCoA_ERROR(ERR::DivByZero, "InvMod");
        return 0;
      }
      if (r == 0) return p-cm;
      return cr;
    }

  } // end of anonymous namespace

  // JAA believes that this routine always returns a reduced value,
  // but does not have a proof of this currently.
  // If r has no inverse mod m then return 0 (if ErrorAction == RtnZeroOnError),
  // o/w throw ERR::DivByZero (if ErrorAction == ThrowOnError).
  unsigned long InvModNoArgCheck(unsigned long r, unsigned long m, const ErrorActionIndicator ErrorAction /*= ThrowOnError*/)
  {
//    CoCoA_ASSERT(r >= 0);
    CoCoA_ASSERT(r < m);
    CoCoA_ASSERT(m >= 2);
    CoCoA_ASSERT(m <= static_cast<unsigned long>(numeric_limits<long>::max())); // check fits into SIGNED long
    // special case if modulus fits into "unsigned int"
    if (m <= numeric_limits<unsigned int>::max())
      return InvModNoArgCheck_int(r, m);
    const unsigned long p = m;
    unsigned long cr = 1;
    unsigned long cm = 0; // this is minus what you think it is!
    while (r != 0)
    {
      {
        const unsigned long q = m/r;
        m -= q*r;
        cm += q*cr;
      }
      if (m == 0) break;
      {
        const unsigned long q = r/m;
        r -= q*m;
        cr += q*cm;
      }
    }
    if (r+m != 1)
    {
      if (ErrorAction == ThrowOnError) CoCoA_ERROR(ERR::DivByZero, "InvMod");
      return 0;
    }
    if (r == 0) return p-cm;
    return cr;
  }

  
  long InvMod(const MachineInt& r, const MachineInt& m, const ErrorActionIndicator ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)");
    const unsigned long modulus = uabs(m);
    unsigned long residue = LeastNNegRemainder(r,m);
    return InvModNoArgCheck(residue, modulus, ErrorAction);
  }

  long InvMod(const BigInt& r, const MachineInt& m, const ErrorActionIndicator ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    const long residue = LeastNNegRemainder(r,m);
    return InvModNoArgCheck(residue, AsUnsignedLong(m), ErrorAction);
  }

  BigInt InvMod(const MachineInt& r, const BigInt& m, const ErrorActionIndicator ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    const BigInt residue(r);
    BigInt ans;
    const int InverseExists = mpz_invert(mpzref(ans), mpzref(residue), mpzref(m));
    if (InverseExists) return ans;
    // Not invertible, so handle error.
    if (ErrorAction == ThrowOnError) CoCoA_ERROR(ERR::DivByZero, "InvMod");
    return BigInt(0);
  }

  BigInt InvMod(const BigInt& r, const BigInt& m, const ErrorActionIndicator ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    BigInt ans;
    const int InverseExists = mpz_invert(mpzref(ans), mpzref(r), mpzref(m));
    if (InverseExists) return ans;
    // Not invertible, so handle error.
    if (ErrorAction == ThrowOnError) CoCoA_ERROR(ERR::DivByZero, "InvMod");
    return BigInt(0);
  }


  long lcm(const MachineInt& a, const MachineInt& b)
  {
    if (IsZero(a) || IsZero(b)) return 0;
    const unsigned long AoverG = uabs(a)/gcd(a,b);
    if (uabs(b) > numeric_limits<long>::max()/AoverG)
      CoCoA_ERROR(ERR::ArgTooBig, "lcm(MachineInt,MachineInt)");
    return AoverG*uabs(b); // implicit case to long is safe (see overflow check above)
  }

  const BigInt lcm(const BigInt& A, const MachineInt& b)
  {
    return lcm(A, BigInt(b));
  }

  const BigInt lcm(const MachineInt& a, const BigInt& B)
  {
    return lcm(BigInt(a), B);
  }

  const BigInt lcm(const BigInt& a, const BigInt& b)
  {
    BigInt ans;
    mpz_lcm(mpzref(ans), mpzref(a), mpzref(b));
    return ans;
  }


// Do we want a procedural form???
//   void lcm(BigInt& lhs, const BigInt& a, const BigInt& b)
//   {
//     mpz_lcm(mpzref(lhs), mpzref(a), mpzref(b));
//   }


  //------------------------------------------------------------------

  long radical(const MachineInt& n)
  {
    if (IsZero(n)) return 0;
    const factorization<long> FacInfo = factor(n);
    const vector<long>& facs = FacInfo.myFactors();
    long ans = 1;
    const int nfacs = len(facs);
    for (int i=0; i < nfacs; ++i)
      ans *= facs[i];
    return ans;
  }

  BigInt radical(const BigInt& N)
  {
    if (IsZero(N)) return N;
    if (IsProbPrime(abs(N))) return abs(N);
    const factorization<BigInt> FacInfo = factor(N);
    const vector<BigInt>& facs = FacInfo.myFactors();
    BigInt ans(1);
    const int nfacs = len(facs);
    for (int i=0; i < nfacs; ++i)
      ans *= facs[i];
    return ans;
  }
  

  factorization<long> SmoothFactor(const MachineInt& n, const MachineInt& TrialLimit)
  {
    if (IsZero(n))
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(n,TrialLimit):  n must be non-zero");
    if (!IsSignedLong(TrialLimit) || AsSignedLong(TrialLimit) < 1)
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(n,TrialLimit):  TrialLimit must be at least 1 and fit into a machine long");
    if (!IsSignedLong(n))
      CoCoA_ERROR(ERR::ArgTooBig, "SmoothFactor(n,TrialLimit):  number to be factorized must fit into a machine long");
    // Below Pmax is unsigned long so that the code will work even if input TrialLimit is numeric_limits<long>::max()
    const unsigned long Pmax = AsUnsignedLong(TrialLimit);
    unsigned long RemainingFactor = uabs(n);


    // Main loop: we simply do trial divisions by primes up to 31, & by numbers not divisible by 2,3,5,...,31.
    vector<long> factors;
    vector<long> exponents;

    // Use "mostly primes"; faster than using primes
    FastMostlyPrimeSeq TrialFactorSeq;
    unsigned long p = *TrialFactorSeq;

    // NB in line below cannot test p*p <= RemainingFactor directly as the product may overflow.
    while (p <= Pmax && p <= RemainingFactor/p)
    {
      int exp = 0;
      unsigned long rem = RemainingFactor % p;
      while (rem == 0)
      {
        ++exp;
        RemainingFactor /= p;
        rem = RemainingFactor % p;
      }
      if (exp > 0)
      {
        factors.push_back(p);
        exponents.push_back(exp);
      }
      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }
    // if RemainingFactor is non-triv & below limit, add it to the list of factors found.
    if (RemainingFactor > 1 && RemainingFactor <= Pmax)
    {
      factors.push_back(RemainingFactor);
      exponents.push_back(1);
      RemainingFactor = 1;
    }
    if (IsNegative(n))
      return factorization<long>(factors, exponents, -static_cast<long>(RemainingFactor));
    return factorization<long>(factors, exponents, RemainingFactor);
  }


  // This is very similar to the function above -- but I don't see how to share code.
  factorization<BigInt> SmoothFactor(const BigInt& N, const MachineInt& TrialLimit)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(N,TrialLimit):  N must be non-zero");
    if (!IsSignedLong(TrialLimit) || AsSignedLong(TrialLimit) < 2)
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(N,TrialLimit):  TrialLimit must be at least 2 and fit into a machine long");
    // Below Pmax is unsigned long so that the code will work even if input TrialLimit is numeric_limits<long>::max()
    const unsigned long Pmax = AsUnsignedLong(TrialLimit);
    BigInt RemainingFactor = abs(N);

    // Main loop: we simply do trial divisions by primes up to 31 & then by numbers not divisible by 2,3,5,...,31.
    vector<BigInt> factors;
    vector<long> exponents;

    // Use "mostly primes"; faster than using primes
    FastMostlyPrimeSeq TrialFactorSeq;
    unsigned long p = *TrialFactorSeq;

    long LogRemainingFactor = FloorLog2(RemainingFactor);
    long CountDivTests = 0;
    // NB in line below cannot test p*p <= RemainingFactor directly as the product p*p may overflow.
    while (p <= Pmax && p <= RemainingFactor/p)
    {
      int exp = 0;
      BigInt quo,rem;
      quorem(quo,rem, RemainingFactor, p);
      while (rem == 0)
      {
        ++exp;
        RemainingFactor = quo;
        quorem(quo,rem, RemainingFactor, p);
      }
      if (exp > 0)
      {
        factors.push_back(BigInt(p));
        exponents.push_back(exp);
        if (quo < p) break; // quo was set by last "failed" call to quorem
        LogRemainingFactor = FloorLog2(RemainingFactor);
        // If several div tests have found no factors, check whether RemainingFactor is prime...
        if (p > 256 && CountDivTests == LogRemainingFactor && LogRemainingFactor < 64)
        {
          if (IsPrime(RemainingFactor)) break;
        }
        CountDivTests = 0;
      }
      ++CountDivTests;
      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }
    // if RemainingFactor is non-triv & below limit, add it to the list of factors found.
    if (RemainingFactor > 1 && RemainingFactor <= Pmax)
    {
      factors.push_back(RemainingFactor);
      exponents.push_back(1);
      RemainingFactor = 1;
    }
    if (N < 0)
      return factorization<BigInt>(factors, exponents, -RemainingFactor);
    return factorization<BigInt>(factors, exponents, RemainingFactor);
  }

  factorization<BigInt> SmoothFactor(const BigInt& N, const BigInt& TrialLimit)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(N,TrialLimit):  N must be non-zero");
    if (TrialLimit < 2)
      CoCoA_ERROR(ERR::BadArg, "SmoothFactor(N,TrialLimit):  TrialLimit must be at least 2");
    
    // Not implemented for large TrialLimit because it would be hideously slow...
    // A naive implementation could simply copy code from SmoothFactor(N,pmax) above.

    long pmax;
    if (!IsConvertible(pmax, TrialLimit))
      CoCoA_ERROR(ERR::NYI, "SmoothFactor(N,TrialLimit) with TrialLimit greater than largest signed long");
    return SmoothFactor(N, pmax);
  }


  factorization<long> factor(const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_ERROR(ERR::BadArg, "factor(n):  n must be non-zero");
    if (!IsSignedLong(n))
      CoCoA_ERROR(ERR::ArgTooBig, "factor(n):  n must fit into a signed long");
    // Simple rather than efficient.
    if (uabs(n) < 2) return SmoothFactor(n,2);
    return SmoothFactor(n,uabs(n));
  }

  factorization<BigInt> factor(const BigInt& N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::BadArg, "factor(N):  N must be non-zero");
    const long PrimeLimit = FactorBigIntTrialLimit; // defined in config.H
    factorization<BigInt> ans = SmoothFactor(N, PrimeLimit);
    const BigInt& R = ans.myRemainingFactor();
    if (abs(R) == 1) return ans;
    if (abs(R) < power(PrimeLimit,2) || IsPrime(abs(R))) // ??? IsPrime or IsProbPrime ???
    {
      ans.myAppend(R,1);
      ans.myNewRemainingFactor(BigInt(sign(R)));
      return ans;
    }

    // Could check for abs(R) being a perfect power...
    CoCoA_ERROR(ERR::NYI, "factor(N) unimplemented in this case -- too many large factors");
    return ans;
  }


  BigInt SumOfFactors(const MachineInt& n, long k)
  {
    if (IsNegative(n) || IsZero(n)) CoCoA_ERROR(ERR::NotPositive, "SumOfFactors");
    if (IsNegative(k)) CoCoA_ERROR(ERR::NotNonNegative, "SumOfFactors");
    const factorization<long> FacInfo = factor(n);
    const vector<long>& mult = FacInfo.myMultiplicities();
    const int s = len(mult);
    if (k == 0)
    {
      long ans = 1;
      for (int i=0; i < s; ++i)
        ans *= 1+mult[i];
      return BigInt(ans);
    }
    // Here k > 0
    BigInt ans(1);
    const vector<long>& fac = FacInfo.myFactors();
    for (int i=0; i < s; ++i)
      ans *= (power(fac[i], k*(mult[i]+1)) - 1)/(power(fac[i],k) - 1);
    return ans;
  }


  long SmallestNonDivisor(const MachineInt& N)
  {
    if (IsZero(N)) CoCoA_ERROR(ERR::NotNonZero, "SmallestNonDivisor");
    unsigned long n = uabs(N);
    if (IsOdd(n)) return 2;
    FastFinitePrimeSeq TrialDivisorList;
    while (n%(*TrialDivisorList) == 0)
      ++TrialDivisorList;
    return *TrialDivisorList;
  }
  
  long SmallestNonDivisor(const BigInt& N)
  {
    if (IsZero(N)) CoCoA_ERROR(ERR::NotNonZero, "SmallestNonDivisor");
    if (IsOdd(N)) return 2;
    // SLUG! simple rather than quick
    FastMostlyPrimeSeq TrialDivisorList;
    while (N%(*TrialDivisorList) == 0)
      ++TrialDivisorList;
    return *TrialDivisorList;
  }


  bool IsSqFree(const MachineInt& N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::BadArg, "IsSqFree(N):  N must be non-zero");
    if (!IsSignedLong(N))
      CoCoA_ERROR(ERR::BadArg, "IsSqFree(N):  N must be non-zero");
    long n = uabs(N); // implicit cast to long is safe
    if (n < 4) return true;

    // Main loop: we simply do trial divisions by the first few primes, then divisors from NoSmallFactorSeq
    const long Pmax = min(FactorBigIntTrialLimit, MaxSquarableInteger<long>());
    long counter = 0;
    const long LogN = FloorLog2(N);
    const long TestIsPrime = (LogN>1000)?1000000:LogN*LogN;

    // Use "mostly primes"; much faster than using primes.
    FastMostlyPrimeSeq TrialFactorSeq;
    long p = *TrialFactorSeq;

    while (p <= Pmax)
    {
      ldiv_t qr = ldiv(n, p);
      if (p > qr.quot) return true;  // test  equiv to p^2 > N
      if (qr.rem == 0)
      {
        n = qr.quot;  // equiv. to N /= p;
        if (n%p == 0) return false;
      }
      else
      {
        ++counter;
        if (counter == TestIsPrime && IsProbPrime(n))
          return true;
      }

      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }

    if (counter < TestIsPrime && IsProbPrime(n)) return true;
    if (IsPower(n)) return false;
    return true;
  }


  bool3 IsSqFree(BigInt N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::BadArg, "IsSqFree(N):  N must be non-zero");
    N = abs(N);
    if (N < 4) return true3;

    // Main loop: we simply do trial divisions by the first few primes, then divisors from NoSmallFactorSeq
    const long Pmax = FactorBigIntTrialLimit;
    long counter = 0;
    const long LogN = FloorLog2(N);
    const long TestIsPrime = (LogN>1000)?1000000:LogN*LogN;

    // Use "mostly primes"; much faster than using primes.
    FastMostlyPrimeSeq TrialFactorSeq;
    long p = *TrialFactorSeq;

    BigInt quot, rem;
    while (p <= Pmax)
    {
      quorem(quot, rem, N, p);
      if (p > quot) return true3; // test  equiv to p^2 > N
      if (IsZero(rem))
      {
        swap(N, quot);  // equiv. to N /= p;
        if (N%p == 0) return false3;
      }
      else
      {
        ++counter;
        if (counter == TestIsPrime && IsProbPrime(N))
          return true3;
      }

      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }

    if (IsPower(N)) return false3;
    if (counter < TestIsPrime && IsProbPrime(N)) return true3;
    return uncertain3;
  }


  long FactorMultiplicity(const MachineInt& b, const MachineInt& n)
  {
    if (IsZero(n)) CoCoA_ERROR(ERR::NotNonZero, "FactorMultiplicity");
    if (IsNegative(b)) CoCoA_ERROR(ERR::BadArg, "FactorMultiplicity");
    const unsigned long base = AsUnsignedLong(b);
    if (base == 0 || base == 1) CoCoA_ERROR(ERR::BadArg, "FactorMultiplicity");
    unsigned long m = uabs(n);
    long mult = 0;
    while (m%base == 0)
    {
      m /= base;
      ++mult;
    }
    return mult;
  }


  long FactorMultiplicity(const MachineInt& b, BigInt N)
  {
    if (IsZero(N)) CoCoA_ERROR(ERR::NotNonZero, "FactorMultiplicity");
    if (IsNegative(b)) CoCoA_ERROR(ERR::BadArg, "FactorMultiplicity");
    const unsigned long base = AsUnsignedLong(b);
    if (base == 0 || base == 1) CoCoA_ERROR(ERR::BadArg, "FactorMultiplicity");

    long mult = 0;
    if (FloorLog2(N) > 2*numeric_limits<unsigned long>::digits)
    {
      // Case: N is fairly large
      // Repeatedly divide by largest power of prime which fits in ulong:
      unsigned long pwr = base;
      int exp = 1;
      while (numeric_limits<unsigned long>::max()/base >= pwr)
      {
        pwr *= base;
        ++exp;
      }
      BigInt r;
      while (true)
      {
        quorem(N, r, N, pwr);
        if (!IsZero(r)) break;
        mult += exp;
      }
      N = r; // any remaining powers of base are now in r
    }

    // Case: N is fairly small
    while (N%base == 0)
    {
      N /= base;
      ++mult; // BUG??? could conceivably overflow if N is a high power of 2???
    }
    return mult;
  }

  long FactorMultiplicity(const BigInt& B, BigInt N)
  {
    if (IsZero(N)) CoCoA_ERROR(ERR::NotNonZero, "FactorMultiplicity");
    long b;
    if (IsConvertible(b, B)) return FactorMultiplicity(b,N);
    if (B < 2) CoCoA_ERROR(ERR::BadArg, "FactorMultiplicity");
    long mult = 0;
    BigInt r;
    while (true)
    {
      quorem(N, r, N, B);
      if (!IsZero(r)) break;
      ++mult;
    }
    return mult;
  }


  long EulerPhi(const MachineInt& n)
  {
    if (IsZero(n) || IsNegative(n))
      CoCoA_ERROR(ERR::BadArg, "EulerPhi(n):  n must be strictly positive");
    if (!IsSignedLong(n)) CoCoA_ERROR(ERR::ArgTooBig, "EulerPhi(n)");
    const factorization<long> facpows = factor(n);
    const vector<long>& primes = facpows.myFactors();
    long ans = AsSignedLong(n);
    const int NumPrimes = len(primes);
    for (int i=0; i < NumPrimes; ++i)
    {
      const long p = primes[i];
      ans = (ans/p)*(p-1);
    }
    return ans;
  }

  BigInt EulerPhi(const BigInt& N)
  {
    if (N <= 0)
      CoCoA_ERROR(ERR::BadArg, "EulerPhi(N):  N must be strictly positive");
    const factorization<BigInt> facpows = factor(N);
    const vector<BigInt>& primes = facpows.myFactors();

    BigInt ans = abs(N);
    const int NumPrimes = len(primes);
    for (int i=0; i < NumPrimes; ++i)
    {
      const BigInt p = primes[i];
      ans = (ans/p)*(p-1);
    }
    return ans;
  }



  // Taken from Cohen's book (page 24).
  long MultiplicativeOrderModPrime(unsigned long residue, unsigned long p)
  {
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(p-1 <= MaxSquarableInteger<unsigned long>());
    CoCoA_ASSERT(IsCoprime(residue, p));

    if (p ==  2) return 1;
    const factorization<long> facpows = factor(p-1);
    const vector<long>& q = facpows.myFactors();
    const vector<long>& pwr = facpows.myMultiplicities();
    const int n = len(q);
    long e = p-1;
    for (int i=0; i < n; ++i)
    {
      e /= SmallPower(q[i], pwr[i]);
      unsigned long rpower = PowerMod(residue, e, p);
      while (rpower != 1)
      {
        rpower = PowerMod(rpower, q[i], p);
        e *= q[i];
      }
    }
    return e;
  }

  long MultiplicativeOrderModPrimePower(unsigned long residue, unsigned long p, int e)
  {
    CoCoA_ASSERT(e > 0);
    CoCoA_ASSERT(residue >= 1 && residue < p);
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(p-1 <= MaxSquarableInteger<unsigned long>());
    long ord =  MultiplicativeOrderModPrime(residue, p);
    if (e == 1) return ord;
    const unsigned long q = SmallPower(p,e);
    unsigned long rpower = PowerMod(residue, ord, q);
    while (rpower != 1)
    {
      rpower = PowerMod(rpower, p, q);
      ord *= p;
    }
    return ord;
  }


  // Taken from Cohen's book (page 24).
  BigInt MultiplicativeOrderModPrime(const BigInt& residue, const BigInt& p)
  {
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(IsCoprime(residue, p));

// ???BUG NYI: use machine int code if possible???
    if (p == 2) return BigInt(1);
    const factorization<BigInt> facpows = factor(p-1);
    const vector<BigInt>& q = facpows.myFactors();
    const vector<long>& pwr = facpows.myMultiplicities();
    const int n = len(q);
    BigInt e = p-1;
    for (int i=0; i < n; ++i)
    {
      e /= power(q[i], pwr[i]);
      BigInt rpower = PowerMod(residue, e, p);
      while (rpower != 1)
      {
        rpower = PowerMod(rpower, q[i], p);
        e *= q[i];
      }
    }
    return e;
  }

  BigInt MultiplicativeOrderModPrimePower(const BigInt& residue, const BigInt& p, long e)
  {
    CoCoA_ASSERT(e > 0);
    CoCoA_ASSERT(residue >= 1 && residue < p);
    CoCoA_ASSERT(IsPrime(p));
    BigInt ord =  MultiplicativeOrderModPrime(residue, p);
    if (e == 1) return ord;
    const BigInt q = power(p,e);
    BigInt rpower = PowerMod(residue, ord, q);
    while (rpower != 1)
    {
      rpower = PowerMod(rpower, p, q);
      ord *= p;
    }
    return ord;
  }


  long MultiplicativeOrderMod(const MachineInt& residue, const MachineInt& modulus)
  {
    if (IsNegative(modulus) || AsUnsignedLong(modulus) < 2)
      CoCoA_ERROR(ERR::BadArg, "MultiplicativeOrderMod: must have modulus >= 2");
    if (gcd(residue, modulus) != 1)
      CoCoA_ERROR(ERR::BadArg, "MultiplicativeOrderMod: residue must be coprime to modulus");
    const unsigned long m = AsUnsignedLong(modulus);
    if (m-1 > MaxSquarableInteger<unsigned long>())
    {
      // return ConvertTo<long>(MultiplicativeOrderMod(residue, BigInt(modulus)));
      return ConvertTo<long>(MultiplicativeOrderMod(residue, BigInt(modulus)));
    }

    unsigned long r = uabs(residue)%m;
    if (r != 0 && IsNegative(residue)) r = m-r;

    const factorization<long> mfacs = factor(m);
    const vector<long>& factor = mfacs.myFactors();
    const vector<long>& exponent = mfacs.myMultiplicities();
    const int n = len(factor);
    long ord = 1;
    for (int i=0; i < n; ++i)
    {
      unsigned long facpow = SmallPower(factor[i], exponent[i]);
      ord = lcm(ord, MultiplicativeOrderModPrimePower(r%facpow, factor[i], exponent[i]));
    }
    return ord;
  }

  long MultiplicativeOrderMod(const BigInt& residue, const MachineInt& modulus)
  {
    return MultiplicativeOrderMod(residue%uabs(modulus), modulus);
  }

  BigInt MultiplicativeOrderMod(const MachineInt& residue, const BigInt& modulus)
  {
    return MultiplicativeOrderMod(BigInt(residue), modulus);
  }

  BigInt MultiplicativeOrderMod(const BigInt& residue, const BigInt& modulus)
  {
    if (modulus < 2)
      CoCoA_ERROR(ERR::BadArg, "MultiplicativeOrderMod: must have modulus >= 2");
    if (gcd(residue, modulus) != 1)
      CoCoA_ERROR(ERR::BadArg, "MultiplicativeOrderMod: residue must be coprime to modulus");
    unsigned long m;
    if (IsConvertible(m, modulus) && m-1 <= MaxSquarableInteger<unsigned long>())
    {
      return BigInt(MultiplicativeOrderMod(residue%m, m));
    }

    const factorization<BigInt> mfacs = factor(modulus);
    const vector<BigInt>& factor = mfacs.myFactors();
    const vector<long>& exponent = mfacs.myMultiplicities();
    const int n = len(factor);
    BigInt ord(1);
    for (int i=0; i < n; ++i)
    {
      const BigInt facpow = power(factor[i], exponent[i]);
      ord = lcm(ord, MultiplicativeOrderModPrimePower(residue%facpow, factor[i], exponent[i]));
    }
    return ord;
  }

  // BigInt version of MultiplicativeOrder???  Could be VERY SLOW!!!

  long PrimitiveRoot(const MachineInt& pp)
  {
    if (IsNegative(pp) || !IsSignedLong(pp) || !IsPrime(pp))
      CoCoA_ERROR(ERR::BadArg, "PrimitiveRoot(p):  p must be a (positive) prime");

    unsigned long p = AsUnsignedLong(pp);
    if (p == 2) return 1;
    const factorization<long> facpows = factor(p-1);
    const vector<long>& primes = facpows.myFactors();
    const int NumPrimes = len(primes);
    for (unsigned long root = 2; /*empty*/; ++root)
    {
      CoCoA_ASSERT(root < p);
      if (root == 4 || root == 8 || root == 9 || root == 16) continue; // skip the first few prime powers
      bool failed = false;
      for (int i=0; i < NumPrimes; ++i)
      {
	if (PowerMod(root, (p-1)/primes[i], p) == 1)
        { failed = true; break; } // effectively a "continue;" for the outer loop
      }
      if (!failed) return root;
    }
  }

  // Essentially identical to the function above.
  // !!!WARNING: calls factor(P-1) so could be VERY SLOW!!!
  long PrimitiveRoot(const BigInt& P)
  {
    if (P <= 0 || !IsPrime(P))
      CoCoA_ERROR(ERR::BadArg, "PrimitiveRoot(P):  P must be a (positive) prime");

    if (P == 2) return 1;
    const factorization<BigInt> facpows = factor(P-1);
    const vector<BigInt>& primes = facpows.myFactors();
    const int NumPrimes = len(primes);
    for (unsigned long root = 2; /*empty*/; ++root)
    {
      CoCoA_ASSERT(root < P);
      if (root == 4 || root == 8 || root == 9 || root ==16) continue; // skip the first few prime powers
      bool failed = false;
      for (int i=0; i < NumPrimes; ++i)
      {
	if (PowerMod(root, (P-1)/primes[i], P) == 1)
        { failed = true; break; } // effectively a "continue;" for the outer loop
      }
      if (!failed) return root;
    }
  }


  //////////////////////////////////////////////////////////////////////
  BigRat SimplestBigRatBetween(const BigRat& lo, const BigRat& hi)
  {
    if (IsZero(lo) || IsZero(hi)) return BigRat(0,1);
    if (sign(lo) != sign(hi)) return BigRat(0,1);
    if (lo == hi) return lo;
    if (sign(lo) < 0) return -SimplestBigRatBetween(-hi, -lo);
    if (lo > hi) return SimplestBigRatBetween(hi,lo);
    // Now we have 0 < lo < hi.
    if (IsOneDen(lo) || IsOneDen(hi)) return BigRat(ceil(lo), 1);
    // Now lo & hi have the same integer part.
    ContFracIter L(lo);
    ContFracIter R(hi);
    ContFracApproximant ans;
    while (!IsEnded(L) && !IsEnded(R) && quot(L) == quot(R))
    {
      ans.myAppendQuot(quot(L));
      ++L;
      ++R;
    }
    if (IsEnded(L)) return lo;
    if (IsEnded(R)) return hi;
    if (quot(L) < quot(R))
    {
      if (IsFinal(L))
        ans.myAppendQuot(quot(L));
      else
        ans.myAppendQuot(quot(L)+1);
    }
    else
    {
      if (IsFinal(R))
        ans.myAppendQuot(quot(R));
      else
        ans.myAppendQuot(quot(R)+1);
    }
    return ans.myRational();
  }


  // This is really a bit messy/involved; isn't there a neater way?
  BigRat SimplestBinaryRatBetween(const BigRat& lo, const BigRat& hi)
  {
    if (IsZero(lo) || IsZero(hi)) return BigRat(0,1);
    if (sign(lo) != sign(hi)) return BigRat(0,1);
    if (lo == hi) return lo;
    if (sign(lo) < 0) return -SimplestBinaryRatBetween(-hi, -lo);
    if (lo > hi) return SimplestBinaryRatBetween(hi,lo);

    // Now we have 0 < lo < hi.
    // First check whether there is an integer between lo and hi.
    const BigInt L = ceil(lo);
    const BigInt H = floor(hi);
    if (H >= L)
    {
      // There is an integer in the interval.
      if (H == L) return BigRat(L,1);
      const long exp = FloorLog2(H-L);
      const BigInt pwr2 = power(2,exp);
      const BigInt q = H/pwr2; // integer division!
      // BigInt q,r;
      // quorem(q,r, H,pwr2);
      // if (IsZero(r)) // both L and H are divisible by pwr2
      // {
      //   if (IsEven(L/pwr2)) return BigRat(L,1); else return BigRat(H,1);
      // }
      if (IsEven(q) || (q-1)*pwr2 < L)
        return BigRat(q*pwr2,1);
      return BigRat((q-1)*pwr2,1);
    }

    // lo & hi have the same integer part.
    const BigRat recip = 1/(hi-lo);
    long exp = FloorLog2(recip);
    if (!IsOneDen(recip) || !IsPowerOf2(num(recip))) ++exp;
    const BigInt lwb = ceil(power(2,exp)*lo);
    const BigInt upb = floor(power(2,exp)*hi);
    if (lwb == upb) return BigRat(lwb,power(2,exp));
    CoCoA_ASSERT(lwb+1 == upb);
    if (IsEven(lwb)) return BigRat(lwb,power(2,exp));
    return BigRat(upb,power(2,exp));
  }

  //////////////////////////////////////////////////////////////////////

  void CRTMill::myAddInfo(const MachineInt& res, const MachineInt& mod, CoprimeFlag check)
  {
    if (IsNegative(mod) ) CoCoA_ERROR(ERR::NotNonNegative, "CRTMill::myAddInfo");
    if (!IsSignedLong(mod) || !IsSignedLong(res)) CoCoA_ERROR(ERR::ArgTooBig, "CRTMill::myAddInfo");
    const long m = AsSignedLong(mod);
    if (check == CheckCoprimality && gcd(m,myM) != 1) CoCoA_ERROR("new modulus not coprime", "CRTMill::myAddInfo");
    CoCoA_ASSERT(IsCoprime(m,myM));
    long r = uabs(res)%m;
    if (r != 0 && IsNegative(res)) r = m-r;
    const long a = myR%m;
    long k;
    if (m <= MaxSquarableInteger<long>())
      k = SymmRemainder((r-a)*InvMod(myM,m),m); // if both (r-a) & InvMod(..) are SymmRem then can allow m <= 2*MaxSquarableInteger
    else
      k = SymmRemainder(BigInt(r-a)*InvMod(myM,m),m); // use BigInt to avoid possible overflow (don't care about speed)
    myR += k*myM;
    myM *= m;
  }

  void CRTMill::myAddInfo(const BigInt& r, const BigInt& m, CoprimeFlag check)
  {
//???JAA    CoCoA_ASSERT(m > 1);
    if (IsOne(m)) return; //???JAA
    if (IsOne(m)) { if (!IsOne(myM)) CoCoA_ERROR(ERR::BadModulus, "CRTMill::myAddInfo"); return; }
    CoCoA_ASSERT((r < m) && (r > -m));
    if (check == CheckCoprimality && gcd(m,myM) != 1) CoCoA_ERROR("new modulus not coprime", "CRTMill::myAddInfo");
    CoCoA_ASSERT(IsCoprime(m,myM));
    const BigInt a = myR%m;
    const BigInt k = SymmRemainder((r-a)*InvMod(myM,m), m); ///???BUG/SLUG??? SymmRemainder(r-a,m) ???
    myR += k*myM;
    myM *= m;
  }


  //////////////////////////////////////////////////////////////////////
  std::ostream& operator<<(std::ostream& out, const CRTMill& CRT)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CRT(residue=" << CRT.myR << ", modulus=" << CRT.myM << ")";
    return out;
  }

  //////////////////////////////////////////////////////////////////////

  ContFracIter::ContFracIter(const BigRat& Q)
  {
    myQuot = floor(Q);
    if (!IsOneDen(Q))
      myFrac = 1/(Q - myQuot); // non-negative!
  }

  const BigInt& ContFracIter::operator*() const
  {
    if (IsEnded(*this)) CoCoA_ERROR(ERR::IterEnded, "ContFracIter::operator*");
    return myQuot;
  }

  ContFracIter& ContFracIter::operator++()
  {
    if (IsEnded(*this)) CoCoA_ERROR(ERR::IterEnded, "ContFracIter::operator++");
    myQuot = floor(myFrac);
    if (IsOneDen(myFrac))
      myFrac = 0;
    else
      myFrac = 1/(myFrac - myQuot); // strictly positive!
    return *this;
  }

  ContFracIter ContFracIter::operator++(int)
  {
    ContFracIter prev = *this;
    operator++();
    return prev;
  }


  bool IsEnded(const ContFracIter& CFIter)
  {
    return IsZero(CFIter.myQuot) && IsZero(CFIter.myFrac);
  }


  bool IsFinal(const ContFracIter& CFIter)
  {
    return IsZero(CFIter.myFrac);
  }


  std::ostream& operator<<(std::ostream& out, const ContFracIter& CFIter)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ContFracIter(myFrac = " << CFIter.myFrac
        << ", myQuot = " << CFIter.myQuot << ")";
    return out;
  }


  //////////////////////////////////////////////////////////////////

  ContFracApproximant::ContFracApproximant():
      myCurr(BigRat::OneOverZero), // WARNING: anomalous value, 1/0
      myPrev(0,1)
  {}


  void ContFracApproximant::myAppendQuot(const MachineInt& q)
  {
    // Simple rather than fast.
    myAppendQuot(BigInt(q));
  }

  void ContFracApproximant::myAppendQuot(const BigInt& q)
  {
    // These 9 lines should avoid (all explicit) temporaries:
    // NB I have to use pointers to mpq_t because GMP's design won't let me use references.
    mpq_t* prev = &mpqref(myPrev);
    mpq_t* curr = &mpqref(myCurr);
    mpq_t* next = &mpqref(myNext);
    mpz_mul(mpq_numref(*next), mpq_numref(*curr), mpzref(q));
    mpz_add(mpq_numref(*next), mpq_numref(*next), mpq_numref(*prev));
    mpz_mul(mpq_denref(*next), mpq_denref(*curr), mpzref(q));
    mpz_add(mpq_denref(*next), mpq_denref(*next), mpq_denref(*prev));
    swap(myCurr, myPrev);
    swap(myNext, myCurr);
  }


  std::ostream& operator<<(std::ostream& out, const ContFracApproximant& CFConv)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ContFracApproximant(myCurr = " << CFConv.myCurr
        << ",  myPrev = " << CFConv.myPrev << ")";
    return out;
  }

  //////////////////////////////////////////////////////////////////


  // NB if (Q == 0) then myCFIter starts off "ended"
  CFApproximantsIter::CFApproximantsIter(const BigRat& Q):
      myCFIter(Q),
      myApproximant()
  {
    if (!IsEnded(myCFIter))
      myApproximant.myAppendQuot(quot(myCFIter));
  }

  CFApproximantsIter::CFApproximantsIter(const ContFracIter& CFIter):
      myCFIter(CFIter),
      myApproximant()
  {
    if (!IsEnded(myCFIter))
      myApproximant.myAppendQuot(quot(myCFIter));
  }


  CFApproximantsIter& CFApproximantsIter::operator++()
  {
    if (IsEnded(*this)) CoCoA_ERROR(ERR::IterEnded, "CFApproximantsIter::operator++");
    ++myCFIter;
    if (IsEnded(myCFIter)) return *this;
    myApproximant.myAppendQuot(quot(myCFIter));

    return *this;
  }

  CFApproximantsIter CFApproximantsIter::operator++(int)
  {
    CFApproximantsIter prev = *this;
    operator++();
    return prev;
  }


  std::ostream& operator<<(std::ostream& out, const CFApproximantsIter& CFAIter)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CFApproximantsIter(myApproximant = " << CFAIter.myApproximant
        << ",  myCFIter = " << CFAIter.myCFIter << ")";
    return out;
  }


  // Return first cont frac convergent having rel error at most MaxRelErr
  BigRat CFApprox(const BigRat& q, const BigRat& MaxRelErr)
  {
    // Simple rather than superfast.
    if (MaxRelErr < 0 || MaxRelErr > 1) CoCoA_ERROR(ERR::BadArg, "CFApprox: relative error must be between 0 and 1");
    if (IsZero(q) || IsZero(MaxRelErr)) return q;
    const BigRat MaxAbsErr = abs(q*MaxRelErr);
    if (MaxAbsErr >= 1) return BigRat(floor(q),1);
    CFApproximantsIter CFAIter(q);
    while (abs(q - *CFAIter) > MaxAbsErr)
    {
      ++CFAIter;
    }
    return *CFAIter;
  }


  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction (continued fraction method)

  RatReconstructByContFrac::RatReconstructByContFrac(MachineInt LogEps):
      myCRT(),
      myLogEps(AsSignedLong(LogEps)),
      myResultIsUpToDate(true),
      myResultIsConvincing(false),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {
    if (myLogEps < 3) CoCoA_ERROR(ERR::BadArg, "RatReconstructByContFrac ctor");
    if (myLogEps > 10000) CoCoA_ERROR(ERR::ArgTooBig, "RatReconstructByContFrac ctor");
  }



  void RatReconstructByContFrac::myAddInfo(const MachineInt& r, const MachineInt& m)
  {
    myCRT.myAddInfo(r,m);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }

  void RatReconstructByContFrac::myAddInfo(const BigInt& R, const BigInt& M)
  {
    myCRT.myAddInfo(R,M);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }


  void RatReconstructByContFrac::myUpdateResult() const
  {
    if (myResultIsUpToDate) return;
    CoCoA_ASSERT(!myResultIsConvincing);
//    VerboseLog VERBOSE("RatRecon");
    myResultIsUpToDate = true;
    const BigInt& X = CombinedResidue(myCRT);
    const BigInt& M = CombinedModulus(myCRT);
//    VERBOSE(50) << "X=" << FloatStr(X) << "  M=" << FloatStr(M) << std::endl;

    const long LogM = FloorLog2(M);
    const long LogExtraFactor = myLogEps + FloorLog2(LogM*LogM);
    // Check for zero (allowing some faulty residues)
    if (2*FloorLog2(gcd(X,M)) > LogExtraFactor+LogM)
    {
//      VERBOSE(50) << "ZERO" << std::endl;
      myResult = 0;
      myResultIsConvincing = true;
      return;
    }

    // Main algm
    const BigRat q(X, M);
    BigInt MaxQuot(1);
    for (ContFracIter it(q); !IsEnded(it); ++it)
    {
      MaxQuot = max(MaxQuot, quot(it));
//      if (quot(it) > MaxQuot) MaxQuot = quot(it);
    }

    // Result is "up-to-date".
    // Do first check whether it is "convincing"; if not, return.
    if (FloorLog2(MaxQuot) < myLogEps)
    {
//      VERBOSE(50) << "FAIL-1:  MaxQuot=" << MaxQuot << "  LogEps=" << myLogEps << std::endl;
      return;
    }

    ContFracApproximant CFA;
    for (ContFracIter it(q); quot(it) != MaxQuot; ++it)
      CFA.myAppendQuot(quot(it));

    const BigRat& RS = CFA.myRational();
    myResult = X - M*RS;
    myBadFactor = gcd(M, den(RS));
    // Second convincingness check:
    if (FloorLog2(num(myResult))+FloorLog2(den(myResult))+2*FloorLog2(myBadFactor)+LogExtraFactor > LogM)
    {
//      VERBOSE(50) << "FAIL-2: LogNum=" << FloorLog2(num(myResult)) << "  LogDen=" << FloorLog2(den(myResult)) << "   2*LogBad=" << 2*FloorLog2(myBadFactor) << "  LogExtra=" << LogExtraFactor << "   LogM=" << LogM << std:: endl;
      return;
    }
//    VERBOSE(50) << "SUCCESS  myResult=" << myResult << std::endl;
    myResultIsConvincing = true;
  }


  const BigRat& ReconstructedRat(const RatReconstructByContFrac& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_ERROR("Result is not convincing","ReconstructedRat(RatReconstructByContFrac)");
    return reconstructor.myResult;
  }


  bool IsConvincing(const RatReconstructByContFrac& reconstructor)
  {
    reconstructor.myUpdateResult();
    return reconstructor.myResultIsConvincing;
  }


  const BigInt& BadMFactor(const RatReconstructByContFrac& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_ERROR("Result is not convincing","BadMFactor(RatReconstructByContFrac)");
    return reconstructor.myBadFactor;
  }


  std::ostream& operator<<(std::ostream& out, const RatReconstructByContFrac& reconstructor)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "RatReconstructByContFrac(CRT=" << reconstructor.myCRT
//        << ", threshold=" << reconstructor.myThresholdValue
        << ", LogEps=" << reconstructor.myLogEps
        << ", IsConvincing=" << reconstructor.myResultIsConvincing;
    if (reconstructor.myResultIsConvincing)
      out << ", result=" << reconstructor.myResult;
    out << ")";
    return out;
  }


  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction (2x2 lattice method)


  const long RatReconstructByLattice::ourDefaultSafetyFactor = 4096;

  BigInt RatReconstructByLattice::myCheckSafetyFactor(const BigInt& SafetyFactor)
  {
    // SafetyFactor == 0 --> use default value (see static data mem ourDefaultSafetyFactor)
    if (SafetyFactor < 0) CoCoA_ERROR(ERR::NotNonNegative, "RatReconstructByLattice ctor");
    if (SafetyFactor > 0) return SafetyFactor;
    return BigInt(ourDefaultSafetyFactor);
  }

  RatReconstructByLattice::RatReconstructByLattice(const MachineInt& SafetyFactor):
      myCRT(),
      mySafetyFactor(myCheckSafetyFactor(BigInt(SafetyFactor))),
      myResultIsUpToDate(true),
      myResultIsConvincing(true),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {}

  RatReconstructByLattice::RatReconstructByLattice(const BigInt& SafetyFactor):
      myCRT(),
      mySafetyFactor(myCheckSafetyFactor(SafetyFactor)),
      myResultIsUpToDate(true),
      myResultIsConvincing(true),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {}


  void RatReconstructByLattice::myAddInfo(const MachineInt& r, const MachineInt& m)
  {
    myCRT.myAddInfo(r,m);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }

  void RatReconstructByLattice::myAddInfo(const BigInt& R, const BigInt& M)
  {
    myCRT.myAddInfo(R,M);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }


  void RatReconstructByLattice::myUpdateResult() const
  {
    if (myResultIsUpToDate) return;
    myResultIsUpToDate = true;
    const BigInt R = CombinedResidue(myCRT);
    const BigInt M = CombinedModulus(myCRT);

    BigInt a0 = M;
    BigInt b0;
    BigInt a1 = R;
    BigInt b1(1);

    do
    {
      const BigInt q = round(BigRat(a0*a1+b0*b1, a1*a1+b1*b1));
      a0 -= q*a1;
      b0 -= q*b1;
      swap(a0,a1);
      swap(b0,b1);
    } while (a1*a1+b1*b1 < a0*a0+b0*b0);
    if (mySafetyFactor*(a0*a0+b0*b0) >= M) return; // result is "up-to-date" but NOT "convincing"
    myResult = BigRat(a0,b0);
    myBadFactor = gcd(gcd(a0,b0),M);
    myResultIsConvincing = true;
  }


  const BigRat& ReconstructedRat(const RatReconstructByLattice& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_ERROR("Result is not convincing","ReconstructedRat(RatReconstructByLattice)");
    return reconstructor.myResult;
  }

  bool IsConvincing(const RatReconstructByLattice& reconstructor)
  {
    reconstructor.myUpdateResult();
    return reconstructor.myResultIsConvincing;
  }


  const BigInt& BadMFactor(const RatReconstructByLattice& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_ERROR("Result is not convincing","BadMFactor(RatReconstructByLattice)");
    return reconstructor.myBadFactor;
  }

  std::ostream& operator<<(std::ostream& out, const RatReconstructByLattice& reconstructor)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "RatReconstructByLattice(CRT=" << reconstructor.myCRT
        << ", SafetyFactor=" << reconstructor.mySafetyFactor
        << ", IsConvincing=" << reconstructor.myResultIsConvincing;
    if (reconstructor.myResultIsConvincing)
      out << ", result=" << reconstructor.myResult;
    out << ")";
    return out;
  }

  //////////////////////////////////////////////////////////////////
  // FTRR

  BigInt ComputeMmax(long e, vector<long> mod)
  {
    sort(mod.begin(), mod.end());
    const long s = len(mod);
    BigInt ans(1);
    for (long i=s-e; i < s; ++i)
      ans *= mod[i];
    return ans;
  }

  // BUG BUG BUG  STOPGAP impl; <---> generic impl in tmp.H
  BigInt product(const vector<long>& mod)
  {
    BigInt ans(1);
    const long s = len(mod);
    for (long i=0; i < s; ++i)
      ans *= mod[i];
    return ans;
  }

  BigRat RatReconstructWithBounds(long e, const BigInt& P, const BigInt& Q, const std::vector<long>& res, const std::vector<long>& mod)
  {
    const BigRat FAILURE(BigRat::OneOverZero); // "impossible rational" used to indicate failure
    const long s = len(res);
    if (len(mod) != s) CoCoA_ERROR(ERR::BadArg, "FTRR");
    if (e < 0 || 2*e >= s || P < 1 || Q < 1) CoCoA_ERROR(ERR::BadArg, "FTRR");
    const BigInt Mmax = ComputeMmax(e, mod);
    if (2*P*Q*power(Mmax,2) >= product(mod)) CoCoA_ERROR(ERR::BadArg, "FTRR");

    {
      long CountZeroes = 0;
      for (long i=0; i < s; ++i)
        if (res[i]%mod[i] == 0) ++CountZeroes;
      if (CountZeroes >= s-e) return BigRat(0,1);
    }

    CRTMill CRT;
    for (long i=0; i < s; ++i)
      CRT.myAddInfo(res[i], mod[i]);
    const BigInt X = CombinedResidue(CRT);
    const BigInt M = CombinedModulus(CRT);

    if (gcd(X,M) > P*Mmax) return FAILURE;

    BigInt u1(1); BigInt u2;    BigInt u3(M);
    BigInt v1;    BigInt v2(1); BigInt v3(X);
    while (abs(v2) <= Q*Mmax)
    {
      const BigInt q = u3/v3; // floor division!
      u1 = u1 - q*v1;  swap(u1, v1);
      u2 = u2 - q*v2;  swap(u2, v2);
      u3 = u3 - q*v3;  swap(u3, v3);
    }
    const BigRat r = X + M*BigRat(u1, u2);
    if (abs(num(r)) > P || den(r) > Q) return FAILURE;
    long CountBadModuli=0;
    for (long i=0; i < s; ++i)
      if (gcd(u2, mod[i]) > 1) ++CountBadModuli;
    if (CountBadModuli > e) return FAILURE;
    return r;
  }

  //////////////////////////////////////////////////////////////////
  // Binomial repr of an integer.
  // A fair compromise between simplicity and speed.

  namespace  // file local fn
  {

    BigInt SearchUpwards(const BigInt& N, BigInt n, long r)
    {
      BigInt step(1);
      while (binomial(n+step,r) <= N)
      {
        n += step;
        step *= 2;
      }
      step /= 2; // step is power of 2 (may also be 1)
      while (step >= 1)
      {
        if (binomial(n+step, r) <= N)
          n += step;
        step /= 2;
      }
      return n;
    }

  } // end of anonymous namespace

  std::vector<BigInt> BinomialRepr(BigInt N, long r)
  {
    if (N < 0) CoCoA_ERROR(ERR::NotNonNegative, "BinomialRepr(N,r)");
    if (r < 1) CoCoA_ERROR(ERR::NotPositive, "BinomialRepr(N,r)");
    vector<BigInt> ans(r+1);
    while (r > 0 && N > 0)
    {
      BigInt lwb = (r-1 + iroot(power(2,r)*factorial(r)*N, r))/2;
      if (N <= r || lwb < r) lwb = r;
      if (lwb == r && N > r) lwb = r+1;
      ans[r] = SearchUpwards(N, lwb, r);
      N -= binomial(ans[r], r);
      --r;
    }
    return ans;
  }

  BigInt BinomialReprShift(BigInt N, long r, long shift1, long shift2)
  {
    const vector<BigInt> n = BinomialRepr(N, r);
    BigInt ans;
    for (long i=1; i <= r; ++i)
    {
      if (n[i]==0) continue;
      if (i+shift2 < 0) continue;
      ans += binomial(n[i]+shift1, i+shift2);
    }
    return ans;
  }


  // Procedure fills slots 0 to N (included!)
  // Source: http://www.mathpages.com/home/kmath383.htm
  // IDEA: could easily write a version which EXTENDS an existing table.
  void NumPartitionsTbl(vector<BigInt>& tbl, long N)
  {
    tbl.resize(N+1);
    tbl[0] = 1;
    for (long j=1; j <= N; ++j)
    {
      int sign = 1;
      long k = 1;
      long i = 1;
      BigInt sum;
      while (true)
      {
        if (i > j) break;
        if (sign == 1) sum += tbl[j-i]; else sum -= tbl[j-i];
        i += k; // should not overflow if tbl fits into RAM
        if (i > j) break; 
        if (sign == 1) sum += tbl[j-i]; else sum -= tbl[j-i];
        sign = -sign;
        ++k;
        i += 2*k-1; // should not overflow if tbl fits into RAM
      }
      tbl[j] = sum;
    }
  }

  BigInt NumPartitions(const MachineInt& n)
  {
    if (IsNegative(n)) return BigInt(0); // or give error???
    if (IsZero(n)) return BigInt(1);
    const long N = AsSignedLong(n);
    vector<BigInt> tbl(N+1);
    NumPartitionsTbl(tbl, N);
    return tbl[N];
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory.C,v 1.86 2018/05/22 14:16:39 abbott Exp $
// $Log: NumTheory.C,v $
// Revision 1.86  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.85  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.84  2018/05/05 15:23:46  abbott
// Summary: Added primality test in SmoothFactor if there were several tests without finding a factor, and the RemainingFactor is not too big
//
// Revision 1.83  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.82  2018/04/18 14:27:57  abbott
// Summary: Removed dirty trick from SmoothFactor
//
// Revision 1.81  2018/03/13 16:15:06  abbott
// Summary: Removed ProbPrimeIters (it is already in NumTheory-prime.C)
//
// Revision 1.80  2018/03/13 14:52:33  abbott
// Summary: Cleaned impl of RatReconstructByContFrac; added (commented-out) code for verbose
//
// Revision 1.79  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.78  2018/02/27 10:50:56  abbott
// Summary: Moved stuff about primes into NumTheory_prime
//
// Revision 1.77  2018/02/19 10:17:33  abbott
// Summary: Added comment
//
// Revision 1.76  2018/02/15 17:26:11  abbott
// Summary: Added EratosthenesRange, and PrimeSeq
//
// Revision 1.75  2018/01/17 10:29:52  abbott
// Summary: Changed InvModNoCheck into InvModNoArgCheck
//
// Revision 1.74  2018/01/16 11:42:51  abbott
// Summary: Changed NoThrow into RtnZeroOnError
//
// Revision 1.73  2017/11/08 14:03:56  abbott
// Summary: Added new fn SumOfFactors
//
// Revision 1.72  2017/10/17 15:51:00  abbott
// Summary: Replaced gcd(...)==1 by IsCoprime(...)
//
// Revision 1.71  2017/10/17 15:44:26  abbott
// Summary: Added new fn IsCoprime
//
// Revision 1.70  2017/04/18 09:59:30  bigatti
// -- separated errors in myAddInfo
//
// Revision 1.69  2016/11/23 15:37:25  abbott
// Summary: Now SmallestNonDivisor returns 0 if NextPrime reaches itslimit
//
// Revision 1.68  2016/11/22 14:31:01  abbott
// Summary: Added SmallestNonDivisor
//
// Revision 1.67  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.66  2016/10/25 20:54:09  abbott
// Summary: Added new fn IsSqFree (for BigInt and ringelem of PolyRing over field)
//
// Revision 1.65  2016/10/20 18:04:17  abbott
// Summary: Added radical; improved slightly factor and SmoothFactor
//
// Revision 1.64  2016/10/19 13:42:41  abbott
// Summary: Added radical for integers; fixed a bug in factor for MachineInt
//
// Revision 1.63  2016/07/04 11:05:39  abbott
// Summary: Removed two unnecessary assertions
//
// Revision 1.62  2016/06/29 13:11:37  abbott
// Summary: Added "NoThrow" option to InvMod; much cleaning inside NumTheory.C
//
// Revision 1.61  2016/04/19 08:53:26  bigatti
// -- fixed BinomialReprShift for 0 entries (case (45_2)^(+1) --> 55)
//
// Revision 1.60  2016/03/25 20:04:39  abbott
// Summary: Moved ExtendedEuclideanAlg into anon namespace
//
// Revision 1.59  2015/11/24 12:48:32  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"; also "ILogBase" --> "FloorLog"
//
// Revision 1.58  2015/11/21 19:17:53  abbott
// Summary: Added SimplestBinaryRatBetween; put InvModNoArgCheck_int into anon namespace
//
// Revision 1.57  2015/11/05 19:20:11  abbott
// Summary: Corrected two embarrassing bugs
//
// Revision 1.56  2015/11/05 14:18:19  abbott
// Summary: Changed InvMod to give error rather than return 0 if inv does not exist; added InvModNoArgCheck
//
// Revision 1.55  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.54  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.53  2015/07/30 15:57:43  abbott
// Summary: FactorMultiplicity now allows nonprime base
//
// Revision 1.52  2015/06/29 13:26:04  abbott
// Summary: Changed name "valuation" --> "FactorMultiplicity"; added new fns "BadMFactor"
// Author: JAA
//
// Revision 1.51  2014/10/28 15:12:08  abbott
// Summary: Renamed modulus --> CombinedModulus, residue --> CombinedResidue (for CRTMill)
// Author: JAA
//
// Revision 1.50  2014/09/16 10:41:41  abbott
// Summary: Added new fn eratosthenes (with doc, example, test)
// Author: JAA
//
// Revision 1.49  2014/09/01 16:25:52  abbott
// Summary: New condition for NextProbPrime/PrevProbPrime to use NextPrime/PrevPrime (previously allowed too large values); use new constant FactorBigIntTrialLimit; ExtGcd now gives error with args (0,0)
// Author: JAA
//
// Revision 1.48  2014/08/29 16:04:55  abbott
// Summary: Added optional 3rd arg to myAddInfo (so coprimality check is skipped)
// Author: JAA
//
// Revision 1.47  2014/05/06 13:13:36  abbott
// Summary: Removed useless fn PowerModLargeModulus
// Author: JAA
//
// Revision 1.46  2014/05/02 13:54:06  abbott
// Summary: Simplified ctor interface for RatReconstruct* (need explicit arg 0 for default behaviour)
// Author: JAA
//
// Revision 1.45  2014/04/24 09:54:50  abbott
// Summary: Corrected arg check in RatReconstructWithBounds
// Author: JAA
//
// Revision 1.44  2014/04/15 13:27:19  abbott
// Summary: Changed rtn type of PrimitiveRoot to long (for CoCoA-5/BuiltinOneLiners)
// Author: JAA
//
// Revision 1.43  2014/04/11 13:33:59  abbott
// Summary: Updated several assertions following revision of MaxSquarableInteger
// Author: JAA
//
// Revision 1.42  2014/04/04 10:15:17  abbott
// Summary: Updated to new interface for MaxSquarableInteger
// Author: JAA
//
// Revision 1.41  2014/03/24 12:09:21  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.40  2014/01/16 16:09:54  abbott
// Added NumPartitions.
//
// Revision 1.39  2013/10/15 16:20:03  abbott
// Added valuation.
//
// Revision 1.38  2013/05/21 14:31:45  abbott
// Added BinomialRepr and BinomialReprShift to CoCoALib and CoCoA-5.
//
// Revision 1.37  2013/05/20 15:47:09  abbott
// Added new fn BinomialRepr (placed in NumTheory).
//
// Revision 1.36  2013/03/26 14:58:59  abbott
// Replaced calls to obsolete proc "convert" by calls to "ConvertTo<...>".
//
// Revision 1.35  2013/02/26 11:29:17  abbott
// Added impl of RatReconstructWithBounds
//
// Revision 1.34  2013/02/22 22:43:56  abbott
// Consequential change: new syntax for getting result from CRTMill.
//
// Revision 1.33  2013/02/22 18:56:50  abbott
// Added feature that RatReconstructByContFrac & RatReconstructByLattice
// ctors accept arg 0 to mean "use default value".
//
// Revision 1.32  2013/02/19 18:48:15  abbott
// Added printing for CRTMill and RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.31  2013/02/15 17:46:00  abbott
// Added RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.30  2012/12/12 18:25:06  abbott
// Added new fn IsFinal for ContFracIter.
// Corrected (subtle) bug in SimplestBigRatBetween.
//
// Revision 1.29  2012/12/12 10:39:06  abbott
// Corrected (embarrassing) typo in an assertion.
//
// Revision 1.28  2012/12/11 17:30:30  abbott
// Changed name from SimplestRationalInInterval to SimplestBigRatBetween.
// Also fixed a bug in the impl.
//
// Revision 1.27  2012/12/05 15:09:24  abbott
// Added new class ContFracApproximant.
// Added new fn SimplestBigRatBetween (NB name changed on 2012-12-11).
// Some minor cleaning.
//
// Revision 1.26  2012/12/04 20:14:11  abbott
// Added new class CRTMill.
// Improved impl of ContFracIter class.
// Fixed a bug in CFApproxIter class.
//
// Revision 1.25  2012/10/22 10:30:56  abbott
// SmoothFactor now allows limit to be 1 (previously it wanted limit >= 2).
//
// Revision 1.24  2012/10/05 09:30:35  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.23  2012/09/26 12:51:26  abbott
// Clarified a comment about IsProbPrime.
//
// Revision 1.22  2012/06/29 15:17:27  abbott
// Fixed issue #199.
//
// Revision 1.21  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.20  2012/03/16 15:40:12  abbott
// Merged contents of NumTheoryQQ (continued fraction functions) into NumTheory.
// Merged the doc too.
//
// Revision 1.19  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.18  2011/09/06 13:39:08  abbott
// Minor cleaning: now uses "int" for increments and decrements in NextPrime & PrevPrime.
// Fixed (embarrassing) overflow bug.
//
// Revision 1.17  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.16  2011/03/23 21:00:46  abbott
// Removed FindPrimRoot from NumTheory.H because it was already
// available as PrimitiveRoot (a better name).
// Updated documentation for NumTheory.
//
// Revision 1.15  2011/03/22 20:17:18  abbott
// Added fn FindPrimRoot.
// Merged impls from obsolescent SmallPrime.C.
//
// Revision 1.14  2011/03/16 13:26:36  abbott
// Removed all "unsigned" from fn interfaces, and many unsigned from inside fn impls.
//
// Revision 1.13  2011/01/19 16:12:18  bigatti
// -- added ERR::ModulusLT2
//
// Revision 1.12  2010/03/05 21:37:58  abbott
// Completed implementation of MultiplicativeOrderMod.
//
// Revision 1.11  2010/03/03 15:02:09  abbott
// Corrected bug in PowerMod: forgot to give error for negative power of non-invertible base.
// Added tests of PowerMod to test-NumTheory1.C (and consequent change to expected output).
//
// Revision 1.10  2010/03/03 10:43:34  abbott
// Added PrimitiveRoot for big primes (might be very slow).
// Added MultiplicativeOrderMod (currently very SLUGGY implementation).
//
// Revision 1.9  2009/12/29 22:44:32  abbott
// Removed buggy proxy class ZZ::rtn.
// Consequent changes for function prototypes also in NumTheory.
// Corrected some minor buglets in NumTheory.
//
// Revision 1.8  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.7  2009/12/03 17:41:47  abbott
// Minor correction to a comment.
//
// Revision 1.6  2009/10/08 13:39:47  abbott
// Renamed "round" into "RoundDiv".
// Added some new versions of "RoundDiv".
// Added a test for "RoundDiv".
//
// Revision 1.5  2009/09/29 12:26:13  abbott
// Added patch/workaround for a problem in gmp-4.3.1, so that exgcd gives
// the right answers (gmp-4.3.1 no longer produces reduced cofactors).
//
// Revision 1.4  2009/07/06 12:30:48  abbott
// Fixed incorrect assertion in IsBigPrime.
// Corrected a bug in SmoothFactor (had > instead of >=), and tidied up a little.
//
// Revision 1.3  2009/07/02 16:28:10  abbott
// Fairly extensive change to NumTheory (at least internally and philosophically).
// Improved and cleaned NumTheory.  Added documentation.
// Clarified the exact behaviour of most functions.
//
// Revision 1.2  2009/06/11 14:10:58  abbott
// Added commented out procedural forms for gcd/lcm, in case we should
// later want to activate them.
//
// Revision 1.1  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
//
