//   Copyright (c)  2009-2010,2013,2018  John Abbott and Anna Bigatti

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


#include "CoCoA/BigRatOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"
#include "CoCoA/utils-gmp.H"

#include <cmath>
using std::abs;
using std::floor;
#include <iostream>
using std::ostream;
using std::istream;
#include <limits>
using std::numeric_limits;
#include <sstream>
using std::istringstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

namespace CoCoA
{
  

  // STANDARD ARITHMETIC OPERATIONS

  const BigRat abs(const BigRat& Q)
  {
    BigRat ans;
    mpq_abs(mpqref(ans), mpqref(Q));
    return ans;
  }

  const BigRat operator-(const BigRat& Q)
  {
    BigRat ans;
    mpq_neg(mpqref(ans), mpqref(Q));
    return ans;
  }

  const BigRat operator+(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_add(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator-(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_sub(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator*(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_mul(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  const BigRat operator/(const BigRat& Q1, const BigRat& Q2)
  {
    if (IsZero(Q2))
      CoCoA_ERROR(ERR::DivByZero, "Q1/Q2");
    BigRat ans;
    mpq_div(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }


  const BigRat operator+(const BigRat& Q, const BigInt& N)
  {
    BigRat ans = Q;
    return ans += N;
//   THE LINES BELOW SHOULD BE MORE EFFICIENT (but they're not very readable).
//     BigRat ans;
//     mpz_mul(mpq_denref(ans), mpzref(N), mpq_denref(mpqref(Q)));
//     mpz_add(mpq_numref(mpqref(ans)), mpq_denref(mpqref(ans)), mpq_numref(mpqref(Q)));
//     mpz_set(mpq_denref(mpqref(ans)), mpq_denref(mpqref(Q)));
//     return ans;
  }

  const BigRat operator-(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans -= N;
  }

  const BigRat operator*(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans *= N;
  }

  const BigRat operator/(const BigRat& Q, const BigInt& N)
  {
    if (IsZero(N))
      CoCoA_ERROR(ERR::DivByZero, "Q/N");
    BigRat ans(Q);
    return ans /= N;
  }

  const BigRat operator+(const BigInt& N, const BigRat& Q)
  {
    return Q+N;
  }

  const BigRat operator-(const BigInt& N, const BigRat& Q)
  {
    BigRat ans = Q-N;
    return -ans;
  }

  const BigRat operator*(const BigInt& N, const BigRat& Q)
  {
    return Q*N;
  }

  const BigRat operator/(const BigInt& N, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::DivByZero, "N/Q");
    return BigRat(N,1)/Q;
  }


  const BigRat operator+(const BigRat& Q, const MachineInt& n)
  {
    return Q + BigRat(n,1);
  }

  const BigRat operator-(const BigRat& Q, const MachineInt& n)
  {
    return Q - BigRat(n,1);
  }

  const BigRat operator*(const BigRat& Q, const MachineInt& n)
  {
    return Q * BigRat(n,1);
  }

  const BigRat operator/(const BigRat& Q, const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_ERROR(ERR::DivByZero, "Q/n");
    return Q / BigRat(n,1);
  }


  const BigRat operator+(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) + Q;
  }

  const BigRat operator-(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) - Q;
  }

  const BigRat operator*(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) * Q;
  }

  const BigRat operator/(const MachineInt& n, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::DivByZero, "n/Q");
    return BigRat(n,1) / Q;
  }

  const BigRat power(const BigRat& base, const BigInt& exponent)
  {
    if (exponent >= 0)
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_ERROR(ERR::BadPwrZero, "power(BigRat,BigInt)");
    return BigRat(power(den(base), -exponent), power(num(base), -exponent), BigRat::AlreadyReduced);
  }

  const BigRat power(const BigRat& base, const MachineInt& exponent)
  {
    if (!IsNegative(exponent))
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_ERROR(ERR::BadPwrZero, "power(BigRat,MachineInt)");
    return BigRat(power(den(base), negate(exponent)), power(num(base), negate(exponent)), BigRat::AlreadyReduced);
  }



  bool IsZero(const BigRat& Q)
  {
    return IsZero(num(Q));
  }


  bool IsOne(const BigRat& Q)
  {
    return mpq_cmp_ui(mpqref(Q), 1,1) == 0;
  }


  bool IsMinusOne(const BigRat& Q)
  {
    return mpq_cmp_si(mpqref(Q), -1,1) == 0;
  }


  bool IsOneNum(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_numref(mpqref(Q)), 1) == 0;
  }


  bool IsOneDen(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_denref(mpqref(Q)), 1) == 0;
  }


  int sign(const BigRat& Q)
  {
    return mpq_sgn(mpqref(Q));
  }


  // COMPARISON FUNCTIONS

  int cmp(const BigRat& Q1, const BigRat& Q2)
  {
    return sign(mpq_cmp(mpqref(Q1), mpqref(Q2)));
  }


  int cmp(const BigRat& Q, const BigInt& N)
  {
    return cmp(num(Q), N*den(Q));
  }


  int cmp(const BigInt& N, const BigRat& Q)
  {
    return cmp(N*den(Q), num(Q));
  }


  int cmp(const BigRat& Q, const MachineInt& n)
  {
    if (IsNegative(n))
      return sign(mpq_cmp_si(mpqref(Q), AsSignedLong(n),1));
    return sign(mpq_cmp_ui(mpqref(Q), AsUnsignedLong(n),1));
  }


  int cmp(const MachineInt& n, const BigRat& Q)
  {
    return -cmp(Q, n);
  }


  int CmpAbs(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_cmpabs(mpqref(Q1), mpqref(Q2));
  }

  // The next 4 prefer simplicity over speed.
  int CmpAbs(const BigRat& Q, const BigInt& N)
  {
    return CmpAbs(Q, BigRat(N,1));
  }

  int CmpAbs(const BigInt& N, const BigRat& Q)
  {
    return CmpAbs(BigRat(N,1), Q);
  }

  int CmpAbs(const BigRat& Q, const MachineInt& n)
  {
    return CmpAbs(Q, BigRat(n,1));
  }

  int CmpAbs(const MachineInt& n, const BigRat& Q)
  {
    return CmpAbs(BigRat(n,1), Q);
  }



  bool operator==(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_equal(mpqref(Q1), mpqref(Q2));
  }

  bool operator!=(const BigRat& Q1, const BigRat& Q2)
  {
    return !(Q1 == Q2);
  }

  bool operator> (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) > 0;
  }

  bool operator>=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) >= 0;
  }

  bool operator< (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) < 0;
  }

  bool operator<=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const BigInt& N)
  {
    return IsOneDen(Q) && num(Q) == N;
  }

  bool operator!=(const BigRat& Q, const BigInt& N)
  {
    return !(Q == N);
  }

  bool operator> (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) > 0;
  }

  bool operator>=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) >= 0;
  }

  bool operator< (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) < 0;
  }

  bool operator<=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) <= 0;
  }

                        
  bool operator==(const BigInt& N, const BigRat& Q)
  {
    return Q == N;
  }

  bool operator!=(const BigInt& N, const BigRat& Q)
  {
    return !(Q == N);
  }

  bool operator> (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) > 0;
  }

  bool operator>=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) >= 0;
  }

  bool operator< (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) < 0;
  }

  bool operator<=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const MachineInt& n)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const BigRat& Q, const MachineInt& n)
  {
    return !(Q == n);
  }

  bool operator> (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) > 0;
  }

  bool operator>=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) >= 0;
  }

  bool operator< (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) < 0;
  }

  bool operator<=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) <= 0;
  }

                
  bool operator==(const MachineInt& n, const BigRat& Q)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const MachineInt& n, const BigRat& Q)
  {
    return !(n == Q);
  }

  bool operator> (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) > 0;
  }

  bool operator>=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) >= 0;
  }

  bool operator< (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) < 0;
  }

  bool operator<=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) <= 0;
  }

                        

  // MISCELLANEOUS FUNCTIONS

  double mantissa(const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "mantissa(BigRat)");
    return 0.0;
  }

  long exponent(const BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "exponent(BigRat)");
    return 0;
  }

  // BUG BUG BUG Is this log(Q)  or log(abs(Q))???
  double log(const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_ERROR(ERR::BadArg, "log(Q)");
    return log(num(Q)) - log(den(Q));
  }

  long FloorLog2(const BigRat& Q) { return FloorLogBase(Q,2); }
  long FloorLog10(const BigRat& Q)  { return FloorLogBase(Q,10); }

  long FloorLogBase(const BigRat& Q, const MachineInt& base)
  {
    return FloorLogBase(Q, BigInt(base));
  }

  // BUGLY: Can we merge the two almost identical sections?
  long FloorLogBase(const BigRat& Q, BigInt base)
  {
    base = abs(base);
    if (base < 2) CoCoA_ERROR(ERR::BadArg, "FloorLogBase: base must be at least 2");
    if (IsZero(Q)) CoCoA_ERROR(ERR::BadArg, "FloorLogBase: cannot compute log(0)");
    const double ApproxLog = log(Q)/log(base);
    const double delta = 5 * std::abs(ApproxLog) * numeric_limits<double>::epsilon(); // ***ASSUME*** numerical error in ApproxLog is less than delta
///???BUG not yet fully implemented    if (ApproxLog > numeric_limits<long>::max()) CoCoA_ERROR(ERR::ArgTooBig, "FloorLogBase");
    const long NearestInt = static_cast<long>(std::floor(ApproxLog+0.5)); // probably right but could be too big by 1

    // Easy case if ApproxLog is far from an integer...
    if (std::abs(ApproxLog - NearestInt) > delta)
      return static_cast<long>(std::floor(ApproxLog));

    // Hard case: ApproxLog is very close to integer NearestInt
    if (NearestInt >= 0)
    {
      const BigInt pwr = power(base, NearestInt);
      const int test = cmp(abs(Q), pwr);
      if (test >= 0) return NearestInt;
      return NearestInt-1;
    }
    // else NearestInt < 0
    const BigInt pwr = power(base, -NearestInt);
    const BigRat shifted = abs(Q)*pwr;
    const int test = cmp(shifted, 1);
    if (test >= 0) return NearestInt;
    return NearestInt-1;
  }


  BigInt floor(const BigRat& Q)
  {
    BigInt ans;
    mpz_fdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }

  BigInt ceil(const BigRat& Q)
  {
    BigInt ans;
    mpz_cdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }


  // This impl clearly guarantees that rounding is compatible with RoundDiv!
  BigInt round(const BigRat& Q)
  {
    if (IsOneDen(Q)) return num(Q);
    return RoundDiv(num(Q), den(Q));
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigRatOps.C,v 1.1 2018/05/22 14:15:14 abbott Exp $
// $Log: BigRatOps.C,v $
// Revision 1.1  2018/05/22 14:15:14  abbott
// Summary: Split BigRat into BigRat (class defn) and BigRatOps
//
//
