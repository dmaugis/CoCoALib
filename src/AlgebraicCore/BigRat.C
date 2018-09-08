//   Copyright (c)  2009-2010,2013  John Abbott

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


#include "CoCoA/BigRat.H"

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
  
  BigRat::BigRat()
  {
    mpq_init(myRep);
  }


  BigRat::BigRat(const mpq_t q, CopyFromMPQ_t /*NotUsed*/)
  {
    if (q == 0/*nullptr*/)
      CoCoA_ERROR(ERR::NullPtr, "ctor BigRat(mpq_t)");

    mpq_init(myRep);
    mpq_set(myRep, q);
  }



  BigRat::BigRat(const MachineInt& n1, const MachineInt& n2, ReduceFlag status)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(n1,n2)");
    mpq_init(myRep);
    myAssign(BigInt(n1), BigInt(n2), status);
//    BELOW IS ORIGINAL CODE -- slightly better 'cos does not create temporaries.
//     const bool IsNegativeFraction = IsNegative(n1) ^ IsNegative(n2);
//     mpq_set_ui(myRep, uabs(n1), uabs(n2));
//     if (status == NotReduced)
//       mpq_canonicalize(myRep);
//     else
//       CoCoA_ASSERT(IsCoprime(n1,n2));
//     if (IsNegativeFraction)
//       mpq_neg(myRep, myRep);
  }

  BigRat::BigRat(const MachineInt& n1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(n1,N2)");
    mpq_init(myRep);
    myAssign(BigInt(n1), N2, status);
  }

  BigRat::BigRat(const BigInt& N1, const MachineInt& n2, ReduceFlag status)
  {
    if (IsZero(n2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(N1,n2)");
    mpq_init(myRep);
    myAssign(N1, BigInt(n2), status);
  }

  BigRat::BigRat(const BigInt& N1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_ERROR(ERR::DivByZero, "BigRat(N1,N2)");
    mpq_init(myRep);
    myAssign(N1, N2, status);
  }


  BigRat::BigRat(const std::string& str, ReadFromString_t /*NotUsed*/, ReduceFlag status)
  {
    mpq_init(myRep);
//     if (base != 0 && (base < 2 || base > 36))
//       CoCoA_ERROR(ERR::BadNumBase, "BigRat(string,int)");
    if (mpq_set_str(myRep, str.c_str(), 10) != 0)
      CoCoA_ERROR(ERR::BadArg, "BigRat(string)");
    if (status == NotReduced)
      mpq_canonicalize(myRep);
  }


  BigRat::BigRat(const MantExp2& ME)
  {
    mpq_init(myRep);
    if (IsZero(ME.myMantissa)) return;
    const long exp = ME.myExponent-ME.myNumDigits+1;
    mpz_set(mpq_numref(myRep), mpzref(ME.myMantissa));
    if (ME.mySign == -1) mpq_neg(myRep, myRep);
    if (exp >= 0)
      mpq_mul_2exp(myRep, myRep, exp);
    else
      mpq_div_2exp(myRep, myRep, -exp);
  }
  
  BigRat::BigRat(const MantExp10& ME)
  {
    mpq_init(myRep);
    if (IsZero(ME.myMantissa)) return;
    const long exp = ME.myExponent-ME.myNumDigits+1;
    mpz_set(mpq_numref(myRep), mpzref(ME.myMantissa));
    if (ME.mySign == -1) mpq_neg(myRep, myRep);
    if (exp >= 0)
    {
      const BigInt scale = power(10, exp);
      mpz_mul(mpq_numref(myRep), mpq_numref(myRep), mpzref(scale));
    }
    else
    {
      const BigInt scale = power(10,exp);
      mpz_set(mpq_denref(myRep), mpzref(scale));
      mpq_canonicalize(myRep);
    }
  }


  BigRat::BigRat(OneOverZero_t /*NotUsed*/)
  {
    mpq_init(myRep);
    myAssign(BigInt(1), BigInt(0), AlreadyReduced);  // AlreadyReduced diables check that denom is non-zero!!
  }


  BigRat::BigRat(const BigRat& from)
  {
    mpq_init(myRep);
    mpq_set(myRep, from.myRep);
  }


  BigRat::~BigRat()
  {
    mpq_clear(myRep);
  }


  // NOTE: This is NOT EXCEPTION CLEAN if the GMP fns can throw.
  void BigRat::myAssign(const BigInt& N1, const BigInt& N2, ReduceFlag status/*=NotReduced*/)
  {
    CoCoA_ASSERT(!IsZero(N2) || status == AlreadyReduced);
    const bool IsNegativeFraction = (N1 < 0) ^ (N2 < 0);
    mpz_abs(mpq_numref(myRep), mpzref(N1));
    mpz_abs(mpq_denref(myRep), mpzref(N2));
    if (status == NotReduced)
      mpq_canonicalize(myRep);
    else
      CoCoA_ASSERT(IsCoprime(N1,N2));
    if (IsNegativeFraction)
      mpq_neg(myRep, myRep);
  }


  BigRat& BigRat::operator=(const BigRat& rhs)
  {
    mpq_set(myRep, rhs.myRep);
    return *this;
  }

    // -------- functions that modify at least one argument or `*this' ----------

  BigRat& BigRat::operator+=(const BigRat& rhs)
  {
    mpq_add(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator-=(const BigRat& rhs)
  {
    mpq_sub(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator*=(const BigRat& rhs)
  {
    mpq_mul(myRep, myRep, rhs.myRep);
    return *this;
  }

  BigRat& BigRat::operator/=(const BigRat& rhs)
  {
    if (mpz_sgn(mpq_numref(rhs.myRep)) == 0)
      CoCoA_ERROR(ERR::DivByZero, "q1 /= q2");
    mpq_div(myRep, myRep, rhs.myRep);
    return *this;
  }
                        
    // Same but with RHS a BigInt...
  BigRat& BigRat::operator=(const BigInt& rhs)
  {
    mpq_set_z(myRep, mpzref(rhs));
    return *this;
  }


  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator+=(const BigInt& rhs)
  {
    const BigInt D = BigIntFromMPZ(mpq_denref(myRep));
    const BigInt tmp = rhs*D;
    mpz_add(mpq_numref(myRep), mpq_numref(myRep), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator-=(const BigInt& rhs)
  {
    const BigInt D = BigIntFromMPZ(mpq_denref(myRep));
    const BigInt tmp = rhs*D;
    mpz_sub(mpq_numref(myRep), mpq_numref(myRep), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  BigRat& BigRat::operator*=(const BigInt& rhs)
  {
    return operator*=(BigRat(rhs,1));
  }

  BigRat& BigRat::operator/=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "Q /= N");
    // Could be more efficient if "*this" is 0.
    return operator/=(BigRat(rhs,1));
  }
                        
    // Same but with RHS a MachinInteger...
  BigRat& BigRat::operator= (const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpq_set_si(myRep, AsSignedLong(rhs), 1);
    else
      mpq_set_ui(myRep, AsUnsignedLong(rhs), 1);
    return *this;
  }

  BigRat& BigRat::operator+=(const MachineInt& rhs)
  {
    return operator+=(BigInt(rhs));
  }

  BigRat& BigRat::operator-=(const MachineInt& rhs)
  {
    return operator-=(BigInt(rhs));
  }

  BigRat& BigRat::operator*=(const MachineInt& rhs)
  {
    return operator*=(BigInt(rhs));
  }

  BigRat& BigRat::operator/=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_ERROR(ERR::DivByZero, "Q /= n");
    return operator/=(BigInt(rhs));
  }


  const BigRat& BigRat::operator++()
  {
    mpz_add(mpq_numref(myRep), mpq_numref(myRep), mpq_denref(myRep)); // no need to reduce
    return *this;
  }

  const BigRat& BigRat::operator--()
  {
    mpz_sub(mpq_numref(myRep), mpq_numref(myRep), mpq_denref(myRep));
    return *this;
  }

  const BigRat BigRat::operator++(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator++();
    return ans;
  }

  const BigRat BigRat::operator--(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator--();
    return ans;
  }



  // I/O FUNCTIONS

  string ConvertToString(const BigRat& src, int base/*=10*/)
  {
    if (base < 2 || base > 36)
      CoCoA_ERROR(ERR::BadNumBase, "IsConvertible(string,BigRat,int)");
//    const long digits = SizeInBase(num(src),base) + SizeInBase(den(src),base);
    const long digits = mpz_sizeinbase(mpq_numref(mpqref(src)),base) +
                        mpz_sizeinbase(mpq_denref(mpqref(src)),base);
    vector<char> buffer(digits+3); // +3 to allow for minus sign, "/" character and terminating NUL
    mpq_get_str(&buffer[0], base, mpqref(src));
    return &buffer[0];
  }


  std::ostream& operator<<(std::ostream& out, const BigRat& Q)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << num(Q);
//    if (!IsOneDen(Q))
      if (!IsOne(den(Q)))
      out << "/" << den(Q);
    return out;
  }

  std::istream& operator>>(std::istream& in, BigRat& Q)
  {
    if (!in.good()) return in;
    BigInt N;
    in >> N;
    if (!in.good()) return in;
    Q = N;
    const char SlashOrDot = in.peek();  // might trigger EOF
    if (!in) { in.clear(); return in; }
    if (SlashOrDot != '/' && SlashOrDot != '.')
    {
      return in;
    }
    in.ignore(); // cannot trigger EOF
    const char AfterSlashOrDot = in.peek(); // might trigger EOF
    if (SlashOrDot == '.')
    {
      if (!in.good()) { in.clear(); return in; }
      if (!isdigit(AfterSlashOrDot))
      {
        return in;
      }

      const string AfterDot = ScanUnsignedIntegerLiteral(in);
      const long NumPlaces = len(AfterDot);
      if (NumPlaces == 0) return in;

      const int base = (in.flags() & std::ios::oct) ? 8 : (in.flags() & std::ios::hex) ? 16 : 10;
      istringstream FracDigits(AfterDot);
      if (base == 8) FracDigits >> std::oct;
      if (base == 16) FracDigits >> std::hex;
      BigRat FracPart;
      FracDigits >> FracPart;
      FracPart /= power(base, NumPlaces);
      Q += FracPart;
      return in;
    }

    // Found a slash
    if (!in.good()) { in.clear(); in.putback(SlashOrDot); return in; }
    if (!isdigit(AfterSlashOrDot))
    {
      in.putback(SlashOrDot);
      return in;
    }
    BigInt D;
    in >> D;
///???    if (!in.good()) return in;
    Q /= D; // Might throw DivByZero
    return in;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigRat& Q)
  {
    OMOut->mySendApplyStart();
    OMOut->mySendSymbol("nums1","rational");
    OMOut->mySend(num(Q));
    OMOut->mySend(den(Q));
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  OpenMathInput& operator>>(OpenMathInput& OMIn, BigRat& /*Q*/)
  {
    CoCoA_ERROR(ERR::NYI, "OpenMathInput fn for BigRat");
    return OMIn;
  }


  void swap(BigRat& a, BigRat& b)
  {
    mpq_swap(mpqref(a), mpqref(b));
  }

  const BigInt num(const BigRat& Q)
  {
    return BigIntFromMPZ(mpq_numref(mpqref(Q)));
  }

  const BigInt den(const BigRat& Q)
  {
    return BigIntFromMPZ(mpq_denref(mpqref(Q)));
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigRat.C,v 1.27 2018/05/22 14:16:39 abbott Exp $
// $Log: BigRat.C,v $
// Revision 1.27  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.26  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.25  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.24  2018/04/18 14:16:18  abbott
// Summary: Ctor from mpq_t now checks if arg is nullptr, and if so, throws helpful error.
//
// Revision 1.23  2017/10/17 15:51:00  abbott
// Summary: Replaced gcd(...)==1 by IsCoprime(...)
//
// Revision 1.22  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.21  2016/10/08 19:46:35  abbott
// Summary: Correct handling of oct/hex base when reading a "decimal"
//
// Revision 1.20  2016/08/02 12:49:34  abbott
// Summary: Renamed NumDigits to SizeInBase; updated doc.
//
// Revision 1.19  2016/07/21 14:14:22  abbott
// Summary: Changed op>> so it can read decimals
//
// Revision 1.18  2016/03/25 20:42:37  abbott
// Summary: Renamed utils_gmp to utils-gmp
//
// Revision 1.17  2016/03/25 19:58:40  abbott
// Summary: Added BigRat ctors from MantExp2 and MantExp10 structures
//
// Revision 1.16  2015/11/23 18:21:10  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.15  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.14  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.13  2015/09/17 18:10:53  abbott
// Summary: Cleaned impl of ILogBase; fixed redmine 776
//
// Revision 1.12  2014/07/09 11:41:35  abbott
// Summary: Corrected (embarassing) bug in ILogBase
// Author: JAA
//
// Revision 1.11  2014/07/07 12:11:12  abbott
// Summary: Corrected operator>> (forgot to ignore the '/')
// Author: JAA
//
// Revision 1.10  2014/06/14 19:25:11  abbott
// Summary: Added new fn CmpAbs (for BigRat)
// Author: JAA
//
// Revision 1.9  2014/05/16 12:02:28  abbott
// Summary: Changed comment about fn "round"
// Author: JAA
//
// Revision 1.8  2014/01/28 09:58:30  abbott
// Revised impl of ctor from std::string so that it accepts and respects 2nd arg saying whether the fraction should be canonicalized.
//
// Revision 1.7  2013/05/20 15:50:20  abbott
// Added new ctor for BigRat from mpq_t.
//
// Revision 1.6  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.5  2012/12/12 10:38:35  abbott
// Changed assertion to allow creation of 1/0 if marked as AlreadyReduced.
//
// Revision 1.4  2012/12/04 20:14:49  abbott
// Modified BigRat ctor to allow one to create 1/0 (if specified as AleadyReduced).
//
// Revision 1.3  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.2  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.1  2011/09/23 13:20:35  bigatti
// -- QQ.C renamed into BigRat.C
//
// Revision 1.18  2011/09/06 15:21:53  abbott
// Changed "cmp" functions so that the return value is in {-1,0,+1}.
//
// Revision 1.17  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/23 16:18:38  abbott
// Simplified defn of round; added comment about rounding halves towards zero.
//
// Revision 1.15  2011/08/17 11:57:39  abbott
// Added static_cast to keep compiler quiet.
//
// Revision 1.14  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.13  2011/06/23 16:01:07  abbott
// Removed single arg ctor QQ(MachineInteger), & consequential changes.
//
// Revision 1.12  2011/03/01 15:26:10  abbott
// Improved impl of ILogBase -- faster in most cases.
//
// Revision 1.11  2011/02/25 12:06:51  abbott
// Added new fn IsOneNum; also some minor code cleaning in QQ.C
//
// Revision 1.10  2011/01/14 17:23:19  abbott
// Fixed a minor bug in power.
//
// Revision 1.9  2010/12/26 13:03:16  abbott
// Added ILogBase function (to ZZ & QQ).
//
// Revision 1.8  2010/05/07 14:57:52  abbott
// Two main changes:
//   power(QQ,ZZ) now allows negative exponent
//   renamed QQ::AlreadyNormalized to QQ::AlreadyReduced
//           (and allowed denoms to be negative; the ctor then makes them positive).
//
// Revision 1.7  2010/03/22 11:49:28  abbott
// Added ctor from a string.
//
// Revision 1.6  2010/03/18 16:40:42  abbott
// Added missing include directive.
//
// Revision 1.5  2010/03/18 16:34:10  abbott
// Added new pseudo-ctors for QQ with optional flag to indicate that value is already normalized.
// Added OpenMath I/O operators.
//
// Revision 1.4  2009/10/26 15:39:24  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.3  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for QQs.
// Added example program for QQs.
// Minor correction to convert.C; minor cleaning to ex-ZZ1.C
//
// Revision 1.2  2009/07/06 12:31:26  abbott
// Commented out two unused function arguments (to keep compiler quiet
// when compiling a debugging version).
//
// Revision 1.1  2009/07/02 16:29:42  abbott
// Added new class QQ to represent rational numbers.
// Consequent change to the Makefile.
//
//
