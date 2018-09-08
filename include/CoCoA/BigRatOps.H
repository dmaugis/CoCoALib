#ifndef CoCoA_BigRatOps_H
#define CoCoA_BigRatOps_H

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
//   -------------------------------------------------------------------------


#include "CoCoA/BigRat.H"


namespace CoCoA
{

//  class MantExp10;       // fwd decl -- defined in FloatApprox.H
//  class MantExp2;        // fwd decl -- defined in FloatApprox.H

  // STANDARD ARITHMETIC OPERATIONS

  void swap(BigRat& a, BigRat& b);

  const BigRat abs(const BigRat& Q);
  const BigRat operator-(const BigRat& Q);

  const BigRat operator+(const BigRat& Q1, const BigRat& Q2);
  const BigRat operator-(const BigRat& Q1, const BigRat& Q2);
  const BigRat operator*(const BigRat& Q1, const BigRat& Q2);
  const BigRat operator/(const BigRat& Q1, const BigRat& Q2);

  const BigRat operator+(const BigRat& Q, const BigInt& N);
  const BigRat operator-(const BigRat& Q, const BigInt& N);
  const BigRat operator*(const BigRat& Q, const BigInt& N);
  const BigRat operator/(const BigRat& Q, const BigInt& N);

  const BigRat operator+(const BigInt& N, const BigRat& Q);
  const BigRat operator-(const BigInt& N, const BigRat& Q);
  const BigRat operator*(const BigInt& N, const BigRat& Q);
  const BigRat operator/(const BigInt& N, const BigRat& Q);

  const BigRat operator+(const BigRat& Q, const MachineInt& n);
  const BigRat operator-(const BigRat& Q, const MachineInt& n);
  const BigRat operator*(const BigRat& Q, const MachineInt& n);
  const BigRat operator/(const BigRat& Q, const MachineInt& n);

  const BigRat operator+(const MachineInt& n, const BigRat& Q);
  const BigRat operator-(const MachineInt& n, const BigRat& Q);
  const BigRat operator*(const MachineInt& n, const BigRat& Q);
  const BigRat operator/(const MachineInt& n, const BigRat& Q);

  const BigRat power(const BigRat& base, const BigInt& exponent);
  const BigRat power(const BigRat& base, const MachineInt& exponent);


  bool IsZero(const BigRat& Q);
  bool IsOne(const BigRat& Q);
  bool IsMinusOne(const BigRat& Q);
  bool IsOneNum(const BigRat& Q);
  bool IsOneDen(const BigRat& Q);
  int sign(const BigRat& Q);

  // COMPARISON FUNCTIONS

  int cmp(const BigRat& Q1, const BigRat& Q2);
  int cmp(const BigRat& Q, const BigInt& N);
  int cmp(const BigInt& N, const BigRat& Q);
  int cmp(const BigRat& Q, const MachineInt& n);
  int cmp(const MachineInt& n, const BigRat& Q);

  int CmpAbs(const BigRat& Q1, const BigRat& Q2);
  int CmpAbs(const BigRat& Q, const BigInt& N);
  int CmpAbs(const BigInt& N, const BigRat& Q);
  int CmpAbs(const BigRat& Q, const MachineInt& n);
  int CmpAbs(const MachineInt& n, const BigRat& Q);

  bool operator==(const BigRat& Q1, const BigRat& Q2);
  bool operator!=(const BigRat& Q1, const BigRat& Q2);
  bool operator> (const BigRat& Q1, const BigRat& Q2);
  bool operator>=(const BigRat& Q1, const BigRat& Q2);
  bool operator< (const BigRat& Q1, const BigRat& Q2);
  bool operator<=(const BigRat& Q1, const BigRat& Q2);
			
  bool operator==(const BigRat& Q, const BigInt& N);
  bool operator!=(const BigRat& Q, const BigInt& N);
  bool operator> (const BigRat& Q, const BigInt& N);
  bool operator>=(const BigRat& Q, const BigInt& N);
  bool operator< (const BigRat& Q, const BigInt& N);
  bool operator<=(const BigRat& Q, const BigInt& N);
			
  bool operator==(const BigInt& N, const BigRat& Q);
  bool operator!=(const BigInt& N, const BigRat& Q);
  bool operator> (const BigInt& N, const BigRat& Q);
  bool operator>=(const BigInt& N, const BigRat& Q);
  bool operator< (const BigInt& N, const BigRat& Q);
  bool operator<=(const BigInt& N, const BigRat& Q);
			
  bool operator==(const BigRat& Q, const MachineInt& n);
  bool operator!=(const BigRat& Q, const MachineInt& n);
  bool operator> (const BigRat& Q, const MachineInt& n);
  bool operator>=(const BigRat& Q, const MachineInt& n);
  bool operator< (const BigRat& Q, const MachineInt& n);
  bool operator<=(const BigRat& Q, const MachineInt& n);
      		
  bool operator==(const MachineInt& n, const BigRat& Q);
  bool operator!=(const MachineInt& n, const BigRat& Q);
  bool operator> (const MachineInt& n, const BigRat& Q);
  bool operator>=(const MachineInt& n, const BigRat& Q);
  bool operator< (const MachineInt& n, const BigRat& Q);
  bool operator<=(const MachineInt& n, const BigRat& Q);
			

  // MISCELLANEOUS FUNCTIONS

  double mantissa(const BigRat& Q);
  long exponent(const BigRat& Q);
  double log(const BigRat& Q);
  long FloorLog2(const BigRat& Q);
  long FloorLog10(const BigRat& Q);
  long FloorLogBase(const BigRat& Q, const MachineInt& base);
  long FloorLogBase(const BigRat& Q, BigInt base);

  BigInt floor(const BigRat& Q);   ///< greatest integer <= Q
  BigInt ceil(const BigRat& Q);    ///< least integer >= Q
  BigInt round(const BigRat& Q);   ///< round rational Q to nearest integer (halves round toward +infinity)


  // Inline fns


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigRatOps.H,v 1.1 2018/05/22 14:15:14 abbott Exp $
// $Log: BigRatOps.H,v $
// Revision 1.1  2018/05/22 14:15:14  abbott
// Summary: Split BigRat into BigRat (class defn) and BigRatOps
//
//

#endif
