//   Copyright (c)  2005,2008,2016  John Abbott

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

#include "CoCoA/MatrixOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MatrixFp.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/ToString.H"
#include "CoCoA/apply.H"
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/geobucket.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/time.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/verbose.H"

#include <algorithm>
using std::min;
#include <cmath>
using std::ldexp;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    int signature(const vector<int>& perm)
    {
      int ans = 1;
      int n = len(perm);
      for (int i=0; i < n; ++i)
        for (int j=i+1; j < n; ++j)
          if (perm[i] > perm[j]) ans = -ans;
      return ans;
    }

  } // end of namespace anonymous

  // Naive dense matrix multiplication
  // Currently just creates a DenseMat to contain the answer.
  // BUG: must make this behave "intelligently" when multiplying two sparse matrices
  // (what does "sparse" mean?  e.g. consider  transpose(SparseMat)*SparseMat
  //  or even DiagMat(...)*AnyMat, or even ZeroMat*AnyMat,... lots of cases!!!??? BUG
  matrix operator*(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    if (NumCols(Mleft) != NumRows(Mright))
      CoCoA_ERROR(ERR::BadMatrixSize, "Mat1*Mat2");
    const ring R = RingOf(Mleft);
    if (RingOf(Mright) != R)
      CoCoA_ERROR(ERR::MixedRings, "Mat1*Mat2");
    const long N = NumCols(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mright);
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
      {
        RingElem tmp(R);
        for(long k=0; k < N; ++k)
          tmp += Mleft(i,k)*Mright(k,j);
        SetEntry(ans, i, j, tmp);
      }
    return ans;
  }


  matrix operator+(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const ring R = RingOf(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1+Mat2");
    if (NumCols(Mright) != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1+Mat2");
    if (RingOf(Mright) != R) CoCoA_ERROR(ERR::MixedRings, "Mat1+Mat2");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)+Mright(i,j));
    return ans;
  }


  matrix operator-(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const ring R = RingOf(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1-Mat2");
    if (NumCols(Mright) != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1-Mat2");
    if (RingOf(Mright) != R) CoCoA_ERROR(ERR::MixedRings, "Mat1-Mat2");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)-Mright(i,j));
    return ans;
  }


  // Innermost loop is ugly but rather faster than a "clean" implementation.
  void mul(matrix& lhs, ConstMatrixView M1, ConstMatrixView M2)
  {
    const ring& R = RingOf(lhs);
    if (NumCols(M1) != NumRows(M2))
      CoCoA_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
    if ( NumRows(lhs) != NumRows(M1) || NumCols(lhs) != NumCols(M2) )
      CoCoA_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
    if (RingOf(M1) != R || RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "mul(Mat1, Mat2)");
    // Use of the temporary ans avoids aliasing problems and makes the code exception safe.
    matrix ans(lhs->myZeroClone(RingOf(lhs), NumRows(M1), NumCols(M2)));
    RingElem sum(R), prod(R);
    for (long i=0; i < NumRows(M1); ++i)
      for (long j=0; j < NumCols(M2); ++j)
      {
        sum = 0;
        for (long k=0; k < NumCols(M1); ++k)
        {
          // Next 2 lines just do:  sum += M1(i,k) * M2(k,j);
          R->myMul(raw(prod), raw(M1(i,k)), raw(M2(k,j)));
          R->myAdd(raw(sum), raw(sum), raw(prod));
        }
        SetEntry(ans, i, j, sum);
      }

    // The answer is in ans; now swap it with the entries of lhs.
    swap(lhs, ans);
  }

  matrix operator*(ConstRefRingElem x, ConstMatrixView M)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "x*M");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, x*M(i,j));
    return ans;
  }
  
  matrix operator*(const BigRat& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const BigInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const MachineInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }


  // Separate impl in case ring is not commutative
  matrix operator*(ConstMatrixView M, ConstRefRingElem x)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "x*M");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)*x);
    return ans;
  }

  matrix operator*(ConstMatrixView M, const BigRat& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const BigInt& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const MachineInt& x)
  { return M * RingElem(RingOf(M),x); }


  matrix operator-(const ConstMatrixView& M)
  {
    return RingElem(RingOf(M),-1)*M;
  }


  matrix operator/(ConstMatrixView M, ConstRefRingElem x)
  {
    if (IsZeroDivisor(x)) CoCoA_ERROR(ERR::DivByZero, "Mat/x");
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "Mat/x");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)/x);
    return ans;
  }

  matrix operator/(ConstMatrixView M, const BigRat& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const BigInt& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const MachineInt& x)
  { return M / RingElem(RingOf(M),x); }


  matrix power(ConstMatrixView M, long n)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (Nrows != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "power(M,n)");
    if (n == numeric_limits<long>::min()) CoCoA_ERROR(ERR::ExpTooBig, "power(M,n)");
    if (n < 0) return power(inverse(M), -n); // cannot overflow because we have excluded n == MinLong
    if (n == 0) return NewDenseMat(IdentityMat(R,Nrows));

    // An iterative implementation of binary powering.
    long bit = 1; while (bit <= n/2) bit <<= 1;
    matrix ans = NewDenseMat(M);
    while (bit > 1)
    {
      mul(ans,ans,ans);//ans *= ans;
      bit >>= 1;
      if (n&bit) mul(ans,ans,M);//ans *= M;
    }
    return ans;
  }


  matrix power(ConstMatrixView M, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::ExpTooBig, "power(M,N)");
    return power(M, n);
  }


  // The square of the Frobenius norm of a matrix.
  RingElem FrobeniusNorm2(ConstMatrixView A)
  {
    RingElem FrNorm2 = zero(RingOf(A));
    for (long i=0; i < NumRows(A); ++i)
      for (long j=0; j < NumCols(A); ++j)
        FrNorm2 += A(i,j)*A(i,j);
    return FrNorm2;
  }


  // Compute the induced infty-norm of a matrix
  RingElem OperatorNormInfinity(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsOrderedDomain(R))
      CoCoA_ERROR(ERR::NotOrdDom, "OperatorNormInfinity(Mat)");
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    RingElem MaxRowSize = zero(R);

    for (long i=0; i < Nrows; ++i)
    {
      RingElem RowSize(zero(R));
      for (long j=0; j < Ncols; ++j) {	RowSize += abs(M(i,j)); }
      if (RowSize > MaxRowSize)
        MaxRowSize = RowSize;
    }
    return MaxRowSize;
  }

  RingElem OperatorNorm1(ConstMatrixView M)
  {
    return OperatorNormInfinity(transpose(M));
  }


  RingElem det(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det(Mat)");
    RingElem d(RingOf(M));
    M->myDet(d);
    return d;
  }


  long rk(const ConstMatrixView& M)
  {
    if (!IsIntegralDomain(RingOf(M)))
      CoCoA_ERROR(ERR::NotIntegralDomain, "rk(Mat)");
    return M->myRank();
  }


  matrix inverse(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "inverse(Mat)");
    return InverseByGauss(M);
  }


  matrix adj(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "adj(Mat)");
///JAA???return AdjDirect(M);
///JAA???return AdjByDetOfMinors(M);
    if (!IsIntegralDomain(RingOf(M)) || (IsPolyRing(RingOf(M)) && NumIndets(RingOf(M)) > 1))
      return AdjByDetOfMinors(M);
    if (IsZero(det(M))) return AdjByDetOfMinors(M);
    return AdjByInverse(M);
  }


  // Should restriction to full rank be in the name?
  matrix PseudoInverse(ConstMatrixView M)
  {
    // BUG??? Would it make sense to generalize to non fields???
    const ring R = RingOf(M);
    if (!IsField(R))
      CoCoA_ERROR(ERR::NotField, "PseudoInverse(Mat)");

    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    const long rank = rk(M);
    if (rank < Nrows && rank < Ncols) 
      CoCoA_ERROR(ERR::NYI, "PseudoInverse of non full rank matrix");

    // Easy case: a square full rank matrix
    if (Nrows == Ncols)
      return inverse(M);

    if (Nrows < Ncols)
      return transpose(M)*inverse(M*transpose(M));
    else
      return inverse(transpose(M)*M)*transpose(M);
  }


  matrix LinSolve(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R || NumRows(M) != NumRows(rhs))
      CoCoA_ERROR(ERR::BadArg, "LinSolve");
    if (IsField(R)) return LinSolveByGauss(M, rhs);
    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_ERROR(ERR::NYI, "LinSolve over non-field, non-gcddomain, non-polynomial-ring");
    return LinSolve(M,rhs); // never reached -- just to keep compiler quiet
  }


  matrix LinKer(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (IsField(R)) return LinKerByGauss(M);
    //    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    //    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_ERROR(ERR::NYI, "LinKer over non-field");
    return LinKer(M); // never reached -- just to keep compiler quiet
  }


  // WARNING!! Pivot selection strategy is simple rather than clever!
  long RankAndGauss(matrix& M, const int ToDoCols)
  {
    const ring R = RingOf(M);
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "gauss");  
    if (ToDoCols > NumCols(M)) CoCoA_ERROR(ERR::BadColIndex, "gauss");  
    const long Mrows = NumRows(M);

    long rank = 0;
    for (long j=0; j < ToDoCols; ++j)
    {
      // Look for a pivot in col j.
      long PivotRow=rank;
      while (PivotRow < Mrows && M(PivotRow, j) == 0)
        ++PivotRow;
      if (PivotRow == Mrows) continue; // col was zero, try next col.
      if (PivotRow != rank)  SwapRows(M, rank, PivotRow);
      M->myRowMul(rank, 1/M(rank, j));  // make pivot entry = 1
      for (long i=0; i < Mrows; ++i)
      {
        if (i == rank) continue;
        if (M(i, j) == 0) continue;
        M->myAddRowMul(i, rank, -M(i,j));
      }
      ++rank;
    }
    return rank;
  }


  matrix LinSolveByGauss(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R || NumRows(M) != NumRows(rhs))
      CoCoA_ERROR(ERR::BadArg, "LinSolveByGauss");
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "LinSolveByGauss");
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);
    const long RHScols = NumCols(rhs);
    matrix tmp = NewDenseMat(ConcatHor(M, rhs));

    // Do row redn by Gauss
    const long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, RHScols); // initially full of zeroes
    long col=0;
    for (long i=0; i < rank; ++i)
    {
      while (tmp(i,col) == 0) { ++col; }
      for (long j=0; j < RHScols; ++j)
        SetEntry(ans, col, j, tmp(i, j+Mcols));
    }
    for (long i=rank; i < Mrows; ++i)
    {
      for (long j=0; j < RHScols; ++j)
        {
          if (tmp(i, j+Mcols) != 0)
            return NewDenseMat(R,0,0); // to indicate that no soln exists
        }
    }
    return ans;
  }


  matrix LinKerByGauss(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "LinKerByGauss");
    matrix tmp = NewDenseMat(M);
  
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);

    // Do row redn by Gauss
    const long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, Mcols-rank); // initially full of zeroes
    long row=0;
    long anscol=0;
    vector<long> PivotCols; // the i-th pivot is in col PivotCols[i]
    ConstMatrixView Z = ZeroMat(R, std::max(Mcols-Mrows, (long)0), Mcols);
    ConstMatrixView SqTmp = ConcatVer(tmp, Z); // make it square
    for (long j=0; j < Mcols; ++j) // we consider only Mcols x Mcols
      if (SqTmp(row,j) != 0) // j-th col with pivot
      {
        PivotCols.push_back(j);
        ++row;
      }
      else // j-th col without pivot
      {
        for (long i=0; i < len(PivotCols); ++i)  // "copy" j-th column
          SetEntry(ans, PivotCols[i], anscol, -SqTmp(i, j));
        SetEntry(ans, j, anscol, one(R));
        ++anscol;
      }
    return ans;
  }


  matrix LinSolveByHNF(ConstMatrixView M, ConstMatrixView rhs)
  {
    // HNF works only for PIDs: i.e. ZZ or k[x]
    // NB Could work in k[x,y,z] if M is univariate!
    CoCoA_ERROR(ERR::NYI, "LinSolveByHNF");
    return NewDenseMat(RingOf(M),0,0); // never reached -- just to keep compiler quiet
  }

  matrix LinSolveByModuleRepr(ConstMatrixView M, ConstMatrixView rhs)
  {
    // Works in k[x,y,z] where k is a field.  Probably slow.
    CoCoA_ERROR(ERR::NYI, "LinSolveByModuleRepr");
    return NewDenseMat(RingOf(M),0,0); // never reached -- just to keep compiler quiet
  }

  /***************************************************************************/


  RingElem det2x2(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det2x2(Mat)");
    if (NumRows(M) != 2)
      CoCoA_ERROR(ERR::BadRowIndex, "det2x2(Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    return M(0,0)*M(1,1) - M(0,1)*M(1,0);
  }
  

  RingElem det3x3(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det3x3(Mat)");
    if (NumRows(M) != 3)
      CoCoA_ERROR(ERR::BadRowIndex, "det3x3(Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    return M(2,2) *(M(0,0)*M(1,1) - M(0,1)*M(1,0)) 
         - M(2,1) *(M(0,0)*M(1,2) - M(0,2)*M(1,0))
         + M(2,0) *(M(0,1)*M(1,2) - M(0,2)*M(1,1));
  }
  

  RingElem det4x4(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det4x4(Mat)");
    if (NumRows(M) != 4)
      CoCoA_ERROR(ERR::BadRowIndex, "det4x4(Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));

    // 2x2 dets for rows 2,3
    const RingElem det_col01 = M(2,0)*M(3,1) - M(2,1)*M(3,0);
    const RingElem det_col02 = M(2,0)*M(3,2) - M(2,2)*M(3,0);
    const RingElem det_col03 = M(2,0)*M(3,3) - M(2,3)*M(3,0);
    const RingElem det_col12 = M(2,1)*M(3,2) - M(2,2)*M(3,1);
    const RingElem det_col13 = M(2,1)*M(3,3) - M(2,3)*M(3,1);
    const RingElem det_col23 = M(2,2)*M(3,3) - M(2,3)*M(3,2);

    // 3x3 dets for rows 1,2,3
    const RingElem det_col012 = M(1,0)*det_col12 - M(1,1)*det_col02 + M(1,2)*det_col01;
    const RingElem det_col013 = M(1,0)*det_col13 - M(1,1)*det_col03 + M(1,3)*det_col01;
    const RingElem det_col023 = M(1,0)*det_col23 - M(1,2)*det_col03 + M(1,3)*det_col02;
    const RingElem det_col123 = M(1,1)*det_col23 - M(1,2)*det_col13 + M(1,3)*det_col12;

    // 4x4 det
    const RingElem det = M(0,0)*det_col123 - M(0,1)*det_col023 + M(0,2)*det_col013 - M(0,3)*det_col012;
    return det;
  }


  RingElem det5x5(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det5x5(Mat)");
    if (NumRows(M) != 5)
      CoCoA_ERROR(ERR::BadRowIndex, "det5x5(Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));

    // 2x2 dets for rows 3,4
    const RingElem det_col01 = M(3,0)*M(4,1) - M(3,1)*M(4,0);
    const RingElem det_col02 = M(3,0)*M(4,2) - M(3,2)*M(4,0);
    const RingElem det_col03 = M(3,0)*M(4,3) - M(3,3)*M(4,0);
    const RingElem det_col04 = M(3,0)*M(4,4) - M(3,4)*M(4,0);
    const RingElem det_col12 = M(3,1)*M(4,2) - M(3,2)*M(4,1);
    const RingElem det_col13 = M(3,1)*M(4,3) - M(3,3)*M(4,1);
    const RingElem det_col14 = M(3,1)*M(4,4) - M(3,4)*M(4,1);
    const RingElem det_col23 = M(3,2)*M(4,3) - M(3,3)*M(4,2);
    const RingElem det_col24 = M(3,2)*M(4,4) - M(3,4)*M(4,2);
    const RingElem det_col34 = M(3,3)*M(4,4) - M(3,4)*M(4,3);

    // 3x3 dets for rows 2,3,4
    const RingElem det_col012 = M(2,0)*det_col12 - M(2,1)*det_col02 + M(2,2)*det_col01;
    const RingElem det_col013 = M(2,0)*det_col13 - M(2,1)*det_col03 + M(2,3)*det_col01;
    const RingElem det_col014 = M(2,0)*det_col14 - M(2,1)*det_col04 + M(2,4)*det_col01;
    const RingElem det_col023 = M(2,0)*det_col23 - M(2,2)*det_col03 + M(2,3)*det_col02;
    const RingElem det_col024 = M(2,0)*det_col24 - M(2,2)*det_col04 + M(2,4)*det_col02;
    const RingElem det_col034 = M(2,0)*det_col34 - M(2,3)*det_col04 + M(2,4)*det_col03;
    const RingElem det_col123 = M(2,1)*det_col23 - M(2,2)*det_col13 + M(2,3)*det_col12;
    const RingElem det_col124 = M(2,1)*det_col24 - M(2,2)*det_col14 + M(2,4)*det_col12;
    const RingElem det_col134 = M(2,1)*det_col34 - M(2,3)*det_col14 + M(2,4)*det_col13;
    const RingElem det_col234 = M(2,2)*det_col34 - M(2,3)*det_col24 + M(2,4)*det_col23;

    // 4x4 dets for rows 1,2,3,4
    const RingElem det_col0123 = M(1,0)*det_col123 - M(1,1)*det_col023 + M(1,2)*det_col013 - M(1,3)*det_col012;
    const RingElem det_col0124 = M(1,0)*det_col124 - M(1,1)*det_col024 + M(1,2)*det_col014 - M(1,4)*det_col012;
    const RingElem det_col0134 = M(1,0)*det_col134 - M(1,1)*det_col034 + M(1,3)*det_col014 - M(1,4)*det_col013;
    const RingElem det_col0234 = M(1,0)*det_col234 - M(1,2)*det_col034 + M(1,3)*det_col024 - M(1,4)*det_col023;
    const RingElem det_col1234 = M(1,1)*det_col234 - M(1,2)*det_col134 + M(1,3)*det_col124 - M(1,4)*det_col123;

    // 5x5 det
    const RingElem det = M(0,0)*det_col1234 - M(0,1)*det_col0234 + M(0,2)*det_col0134 - M(0,3)*det_col0124 + M(0,4)*det_col0123;
    
    return det;
  }


  // Compute det as "alternating" sum of products
  RingElem DetDirect(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    const ring R = RingOf(M);
    CoCoA_ASSERT(IsPolyRing(R)); // !!!BUG BUG BUG!!! TEMPORARY so we can use geobucket

    vector<int> cols(n);
    for (long i=0; i < n; ++i) { cols[i] = i; }

    RingElem det(R);
    geobucket ans(R);
    do {
      RingElem term = one(R);
      for (int i=0; i < n; ++i)
        term *= M(i, cols[i]);
//      det += signature(cols)*term;
      if (IsZero(term)) continue;
      if (NumTerms(term) == 1)
        ans.myAddMulLM(term, RingElem(R,signature(cols)), 1);
      else det += signature(cols)*term;
    } while ( std::next_permutation(cols.begin(),cols.end()) );

    AddClear(det,ans);
    return det;
  }


  // Known defect: this algm is valid only over IntegralDomains
  RingElem DetByGauss(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "DetByGauss(Mat)");
    const long N = NumRows(M);
    const ring R(RingOf(M));
    // ??? this function works only within integral domains!
    if (!IsIntegralDomain(R))
      CoCoA_ERROR(ERR::NotIntegralDomain, "DetByGauss(Mat)");
// Is it worth adding code for this special case?  General code works fine.
//     // Handle 0x0 matrix specially: its det is 1.
//     if (N == 0) { d = one(R); return; }
    const ring K((IsField(R) ? R : NewFractionField(R)));
    matrix Gss(apply(CanonicalHom(R,K), M));
    RingElem c(K);
    RingElem determinant(one(K));
    for (long col=0; col < N; ++col)
    {
      // Pick a good pivot row
      long PivotRow = -1;
      for (long r=col ; r < N; ++r)
      {
        if (!IsZero(Gss(r,col))) { PivotRow = r; break; }
      }

      if (PivotRow == -1)
      {
        return zero(R);
      }
      if (PivotRow != col)
      {
        Gss->mySwapRows(PivotRow,col);
        determinant *= -1;
      }
      c = Gss(col,col);
      determinant *= c;
      for (long i=col+1; i < N; ++i)
        Gss->myAddRowMul(i, col, -Gss(i,col)/c);
    }
    if (R==K) return determinant;
    else return num(determinant)/den(determinant);
  }


  long HadamardRowScale(const ConstMatrixView& M); //fwd decl
  long HadamardColScale(const ConstMatrixView& M); //fwd decl
  BigInt DetBound_GPS(const ConstMatrixView& M);


  RingElem DetByCRT(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetByCRT");
    const ring& ZZ = RingOf(M);
    CoCoA_ASSERT(IsZZ(ZZ));
    const int n = NumRows(M);
    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));

    const HadamardRowCol H = HadamardBound(M); // square of hadamard bound (min of row and col bounds)
    const BigInt CRTBound = 2*FloorSqrt(min(H.myRowBound, H.myColBound)); // factor of 2 to allow for +/- sign
    VERBOSE(80) << "Hrow(sq) = " << FloatStr(H.myRowBound) << "  Hcol(sq) = " << FloatStr(H.myColBound) << "   GPS = " << FloatStr(DetBound_GPS(M)) << std::endl;

 const long RowScale = HadamardRowScale(M);
 const long ColScale = HadamardColScale(M);
 VERBOSE(80) << "RowScale = " << RowScale << "   ColScale = " << ColScale << std::endl;
 VERBOSE(80) << "Scaled Hrow = " << FloatStr(FloorSqrt(H.myRowBound/power(2,RowScale))) <<
   "   Scaled Hcol = " << FloatStr(FloorSqrt(H.myColBound/power(2,ColScale))) << std::endl;

 VERBOSE(80) << "Actually using CRTBound = " << FloatStr(CRTBound) << std::endl;
    int NumPrimes = 0;

    bool UsingHeuristic = false;
    PrimeSeqForCRT PSeq;  // this produces the primes we shall use for CRTing
    CRTMill CRT;
    while (true)
    {
      const bool finished = (CombinedModulus(CRT) > CRTBound);
      const long log2modulus = FloorLog2(CombinedModulus(CRT));
      if (finished || (UsingHeuristic && log2modulus > 50 && (IsZero(CombinedResidue(CRT)) || FloorLog2(CombinedResidue(CRT))+50 < log2modulus)))
      {
        VERBOSE(80) << "Returning; num iters: " << NumPrimes << std::endl;
        return RingElem(ZZ, CombinedResidue(CRT));
      }

      CheckForInterrupt("DetByCRT");
//      p = NextPrime(p);
      const SmallPrime p = NextPrime(PSeq);
      ++NumPrimes;
      VERBOSE(85) << "Iter " << NumPrimes << "  Using prime p=" << p << std::endl;
      SmallFpImpl ModP(p);
      MatrixFp Mp(ModP, n,n);
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          Mp(i,j) = ModP.myReduce(Mcopy[i][j]);

      const long dp = det(Mp);
      CRT.myAddInfo(dp,p, CRTMill::CoprimeModulus);
    }
  }

  // SQUARE of min of row and col Hadamard bounds
HadamardRowCol HadamardBound(const ConstMatrixView& M)
  {
    CoCoA_ASSERT(IsZZ(RingOf(M)));
    const int n = NumRows(M);
    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));

    BigInt RowBound(1);
    for (int i=0; i < n; ++i)
    {
      BigInt L2;
      for (int j=0; j < n; ++j)
        L2 += power(Mcopy[i][j],2);
      RowBound *= L2;
    }

    BigInt ColBound(1);
    for (int j=0; j < n; ++j)
    {
      BigInt L2;
      for (int i=0; i < n; ++i)
        L2 += power(Mcopy[i][j],2);
      ColBound *= L2;
    }

    return HadamardRowCol(RowBound, ColBound);
    // if (RowBound < ColBound)
    //   return RowBound;
    // return ColBound;
  }


  struct IPentry
  {
    IPentry(double cos2, long i1, long i2): cosine2(cos2), row1(i1), row2(i2) {}
    double cosine2;
    long row1;
    long row2;
  };

  inline bool ByCosine2(const IPentry& A, const IPentry& B)
  {
    return A.cosine2 > B.cosine2;
  }

  // Returns an integer K such that abs val of true det is less than
  // row Hadamard bound times 2^(-K).
  // Complexity is CUBIC in matrix dimension!!   ASSUMES ring is ZZ
  long HadamardRowScale(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("HadamardRowScale");
    const double t0 = CpuTime();
    CoCoA_ASSERT(IsZZ(RingOf(M)));
    const int n = NumRows(M);

    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));
    const double t1 = CpuTime();
//    VERBOSE(81) << "Copy time: " << t1-t0 << std::endl;

    vector<BigInt> L2(n);
    for (int i=0; i < n; ++i)
    {
      BigInt L2row;
      for (int j=0; j < n; ++j)
        L2row += power(Mcopy[i][j],2);
      L2[i] = L2row;
    }
    const double t2 = CpuTime();
//    VERBOSE(81) << "L2 time: " << t2-t1 << std::endl;

    vector< vector<double> > Mdbl(n, vector<double>(n));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        long mij_exp;
        long len_exp;
        double mij_mant = mpz_get_d_2exp(&mij_exp, mpzref(Mcopy[i][j]));
        double len_mant = mpz_get_d_2exp(&len_exp, mpzref(L2[i]));
        if (IsOdd(len_exp)) { ++len_exp; len_mant /= 2.0; }
        len_mant = sqrt(len_mant); len_exp /= 2; // exact division!
        if (len_exp - mij_exp > 50) { Mdbl[i][j] = 0.0; continue; }
        const double shift = ldexp(1.0, mij_exp - len_exp);
        Mdbl[i][j] = (mij_mant/len_mant)*shift;
//        VERBOSE(89) << "Entry (" << i << "," << j << ") = " << Mdbl[i][j] << std::endl;
      }
    const double t3 = CpuTime();
//    VERBOSE(81) << "Mdbl time: " << t3-t2 << std::endl;

    vector<IPentry> tbl; tbl.reserve((n*n-n)/2);
    for (int i1=0; i1 < n; ++i1)
    {
      double BiggestInnerProd = 0.0;
      for (int i2=i1+1; i2 < n; ++i2)
      {
        double InnerProd = 0.0;
        for (int j=0; j < n; ++j)
          InnerProd += Mdbl[i1][j]*Mdbl[i2][j];
        tbl.push_back(IPentry(pow(InnerProd,2),i1,i2));
        if (std::abs(InnerProd) > BiggestInnerProd) BiggestInnerProd = std::abs(InnerProd);
      }
      if (i1 == 0)
      {
//        VERBOSE(81) << "BiggestInnerProd = " << BiggestInnerProd << std::endl;
        
        if (BiggestInnerProd < 0.2) { /*VERBOSE(80) << "Not worth it" << std::endl;*/ return 0;} // heuristic - not worth computing redn factor
      }
    }
    const double t4 = CpuTime();
//    VERBOSE(81) << "IP tbl time: " << t4-t3 << std::endl;
    std::sort(tbl.begin(), tbl.end(), ByCosine2);
    const double t5 = CpuTime();
//    VERBOSE(81) << "sort time: " << t5-t4 << std::endl;
    // If all pairs are practically orthog, just return 0 = log(1).
    if (tbl[0].cosine2 < 1.0/n) { /*VERBOSE(80) << "Practically orthog" << std::endl;*/ return 0; }
    
    vector<int> TreeIndex(n); for (int i=0; i < n; ++i) TreeIndex[i] = i;
    int NumTrees = n;
    double RednFactor = 1.0;
    const double eps = 1.0/(1048576.0*1048576.0);
    const double Shift256 = pow(2.0, 256);
    const double LowerLimit = 1.0/Shift256;
    long exp2 = 0;
    int NumBigJumps = 0;
    for (int i=0; i < len(tbl); ++i)
    {
//      if (i < 10) VERBOSE(81) << "cos2 = " << tbl[i].cosine2 << "   (i1,i2)=("<<tbl[i].row1<<","<<tbl[i].row2<<")"<<std:: endl;
      if (tbl[i].cosine2 < 0.01) break;
      const int index1 =TreeIndex[tbl[i].row1];
      const int index2 =TreeIndex[tbl[i].row2];
      if (index1 == index2) continue;
      if (index1 < index2) TreeIndex[tbl[i].row2] = index1;
      else TreeIndex[tbl[i].row1] = index2;
      if (tbl[i].cosine2 < 1.0-eps)
      RednFactor *= (1-tbl[i].cosine2);
      else { RednFactor *= eps; ++NumBigJumps; }
      if (RednFactor < LowerLimit) { RednFactor *= Shift256; exp2 += 256; }
      --NumTrees;
      if (NumTrees == 1) break;
    }
    int exp2rest = 0;
    if (RednFactor != 1) frexp(RednFactor, &exp2rest);
    VERBOSE(81) << "NumBigJumps = " << NumBigJumps << std::endl;
//    VERBOSE(81) << "RednFactor = " << RednFactor << "  exp2rest = " << exp2rest << std::endl;
//    VERBOSE(81) << "Returning " << exp2-exp2rest << std::endl;
//    VERBOSE(81) << "TIME: " << CpuTime()-t0 << std::endl;
    return exp2-exp2rest;  // subtract 1??
  }


  // Lazy man's impl -- at least it's "obviously correct" (right?)
  long HadamardColScale(const ConstMatrixView& M)
  {
    return HadamardRowScale(transpose(M));
  }

  
  // ASSUMES ring is ZZ (or QQ?)
  // SHOULD PREPROCESS MATRIX:
  // (1) make sure that all row sums and all col sums are >= 0
  //     (by changing signs of a row/col --> does not affect abs(det))
  // (2) rescale rows/cols so that most entries are about the same size
  //     (changes det by a known factor).  JAA is not sure how to do this well.
  BigInt DetBound_GPS(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetBound_GPS");
    CoCoA_ASSERT(IsZZ(RingOf(M)));
    const int n = NumRows(M);

    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    vector<BigInt> RowSum(n);
    vector<BigInt> ColSum(n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));
        RowSum[i] += Mcopy[i][j];
        ColSum[j] += Mcopy[i][j];
      }

    BigInt SumAll;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "(INPUT MAT) Entry sum = " << FloatStr(SumAll) << std:: endl;

    // Make sure all RowSums and ColSums are >= 0
    while (true)
    {
      int MinRow = 0; // index of min value
      for (int i=1; i < n; ++i)
        if (RowSum[i] < RowSum[MinRow]) MinRow = i;
      int MinCol = 0; // index of min value
      for (int j=1; j < n; ++j)
        if (ColSum[j] < ColSum[MinCol]) MinCol = j;
      if (RowSum[MinRow] >= 0 && ColSum[MinCol] >= 0) break;
      if (RowSum[MinRow] <= ColSum[MinCol])
      {
        // Negate MinRow
        RowSum[MinRow] = -RowSum[MinRow];
        for (int j=0; j < n; ++j)
        {
          Mcopy[MinRow][j] = -Mcopy[MinRow][j];
          ColSum[j] += 2*Mcopy[MinRow][j];
        }
      }
      else
      {
        // negate MinCol
        ColSum[MinCol] = -ColSum[MinCol];
        for (int i=0; i < n; ++i)
        {
          Mcopy[i][MinCol] = -Mcopy[i][MinCol];
          RowSum[i] += 2*Mcopy[i][MinCol];
        }
      }
    }

    VERBOSE(150) << "After first stage" << std::endl;
    SumAll = 0;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "Entry sum = " << FloatStr(SumAll) << std:: endl;

    // Now check for row-col pairs
    bool changed = true;
    while (changed)
    {
      changed = false;
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
        {
          if (RowSum[i]+ColSum[j] >= 2*Mcopy[i][j]) continue;
          changed = true;
//          VERBOSE(22) << "Pivot " << i << "," << j << std::endl;
          // negate col j
          ColSum[j] = -ColSum[j];
          for (int ii=0; ii < n; ++ii)
          {
            Mcopy[ii][j] = -Mcopy[ii][j];
            RowSum[ii] += 2*Mcopy[ii][j];
          }
          // negate row i
          RowSum[i] = -RowSum[i];
          for (int jj=0; jj < n; ++jj)
          {
            Mcopy[i][jj] = -Mcopy[i][jj];
            ColSum[jj] += 2*Mcopy[i][jj];
          }
        }
    }
    
    VERBOSE(150) << "After second stage" << std::endl;
    SumAll = 0;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "Entry sum = " << FloatStr(SumAll) << std:: endl;

    BigInt SumEntries;
    BigInt SumSquares;
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        SumEntries += Mcopy[i][j];
        SumSquares += power(Mcopy[i][j],2);
      }

    const BigInt delta = power(SumEntries,2) - n*SumSquares;
    VERBOSE(150) << "alpha = " << SumEntries << "  over " << n << std::endl;
    VERBOSE(150) << "beta = " << SumSquares << "  over " << n << std::endl;
    VERBOSE(150) << "delta = " << delta << "  over " << n*n*(n-1) << std::endl;
    if (delta <= 0)
    {
      // Always worse than Hadanard's bound
      if (IsEven(n)) return power(SumSquares,n/2);
      return FloorSqrt(power(SumSquares,n));
    }
    VERBOSE(150) << "Good case" << std::endl;
    if (IsOdd(n))
      return (abs(SumEntries)*power((n*(n-1))*SumSquares-delta, (n-1)/2))/(power(n,n)*power((n-1),(n-1)/2)); // int division
    return FloorSqrt(power(SumEntries,2)*power((n*(n-1))*SumSquares-delta,n-1)/power(n-1,n-1))/power(n,n); // integer division (twice)
  }


  long RankByGauss(std::vector<long>& IndepRows, ConstMatrixView M)
  {
    // ??? this function works only within integral domains!
    const ring R(RingOf(M));
    if (!IsIntegralDomain(R))
      CoCoA_ERROR(ERR::NotIntegralDomain, "RankByGauss(v,Mat)  over non-integral domain");
    const ring K((IsField(R) ? R : NewFractionField(R)));
    const RingHom R2K(R==K ? IdentityHom(K) : EmbeddingHom(K));
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix Gss(NewDenseMat(K, Nrows, Ncols));

    // below copy M via RingHom into Gss
    {
      // This will eventually become a separate function
      for (long row=0; row < Nrows; ++row)
        for (long col=0; col < Ncols; ++col)
          SetEntry(Gss, row, col, R2K(M(row,col)));
    }
    IndepRows.clear();
    RingElem pivot(K);
    long row=0;
    for (long col=0; col < Ncols; ++col)
    {
      if (IsZero(Gss(row,col)))
      {
        long i=row+1;
        for ( ; i < Nrows; ++i)
          if (!IsZero(Gss(i,col)))
          {
            Gss->mySwapRows(i,row);
            break;
          }
        if (i==Nrows) continue;
      }
      IndepRows.push_back(row);
      pivot = Gss(row,col);
      for (long i=row+1; i < Nrows; ++i)
        Gss->myAddRowMul(i, row, -Gss(i,col)/pivot);
      ++row;
      if (row == Nrows) return row;
    }
    return row;
  }


  matrix InverseByGauss(ConstMatrixView M)
  {
    // this code works only if the base ring is an integral domain
    CoCoA_ASSERT(IsIntegralDomain(RingOf(M)));
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "InverseByGauss(Mat)");
    const ring R(RingOf(M));
    if (!IsIntegralDomain(R))
      CoCoA_ERROR(ERR::NotIntegralDomain, "InverseByGauss(Mat)  over non-integral domain");
    const ring    K(IsField(R) ?  R : NewFractionField(R));
    const RingHom R2K(R==K ? IdentityHom(K) : EmbeddingHom(K));
    const long N = NumRows(M);
    matrix Gss(apply(CanonicalHom(R,K), M));
    matrix inv = NewDenseMat(IdentityMat(K, N));
    RingElem c(K);
    for (long j=0; j < N; ++j)
    {
      if (IsZero(Gss(j,j)))
      {
        long i=j+1;
        for ( ; i < N; ++i)
          if (!IsZero(Gss(i,j))) break;
        if (i == N) CoCoA_ERROR(ERR::NotInvMatrix, "InverseByGauss(Mat)");
        Gss->mySwapRows(i,j);
        inv->mySwapRows(i,j);
      }
      c = 1/Gss(j,j);
      Gss->myRowMul(j, c);
      inv->myRowMul(j, c);
      for (long i=0; i < N; ++i)
        if (i != j)
        {
          c = -Gss(i,j);
          Gss->myAddRowMul(i, j, c); // AddRowMul(Gss, i, j, c);
          inv->myAddRowMul(i, j, c);
        }
    }
    if (R==K) return inv;
    matrix inv_R = NewDenseMat(R,N,N);
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        SetEntry(inv_R, i, j, num(inv(i,j))/den(inv(i,j)));
    return inv_R;
  }


//////////////////////////////////////////////////////////////////
// Bareiss

  namespace /*anonymous*/
  {
    RingElem DetByBareiss_ZZ(const ConstMatrixView& M)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)));
      const int n = NumRows(M);
      CoCoA_ASSERT(NumCols(M) == n);
      const ring R = RingOf(M);
      if (n == 0) return one(R);
      if (n == 1) return M(0,0);
      vector< vector< BigInt > > M2(n, vector<BigInt>(n));
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          IsInteger(M2[i][j], M(i,j)); // ignore return value (must be true)
      BigInt d(1);
      int sign = 1;
      for (int k=0; k < n-1; ++k)
      {
        // Find a non-zero pivot in k-th column
        int row = -1;
        for (int i=k; i < n; ++i)
          if (!IsZero(M2[i][k])) { row = i; break; }
        if (row == -1) return zero(R);
        if (row != k) { swap(M2[k], M2[row]); sign = -sign; }
        BigInt tmp; // temporary workspace used in inner loop below
        for (int i=k+1; i < n; ++i)
          for (int j=k+1; j < n; ++j)
          {
//  This block effectively does the following, but is usefully faster (about 4x)
//            M2[i][j] = (M2[i][j]*M2[k][k]-M2[i][k]*M2[k][j])/d;
            mpz_mul(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(M2[k][k]));
            mpz_mul(mpzref(tmp), mpzref(M2[i][k]), mpzref(M2[k][j]));
            mpz_sub(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(tmp));
            mpz_divexact(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(d));
          }
        d = M2[k][k];
      }
      return RingElem(R,sign*M2[n-1][n-1]);
    }
  } // end of anonymous namespace


  RingElem DetByBareiss(const ConstMatrixView& M)
  {
    const ring R = RingOf(M);
    if (IsZZ(R)) return DetByBareiss_ZZ(M); // a bit faster
    CoCoA_ASSERT(IsIntegralDomain(R));
    const int n = NumRows(M);
    CoCoA_ASSERT(NumCols(M) == n);
    if (n == 0) return one(R);
    if (n == 1) return M(0,0);
    matrix M2 = NewDenseMat(M);
    RingElem d = one(R);
    int sign = 1;
    for (int k=0; k < n-1; ++k)
    {
      // Find a non-zero pivot in k-th column
      int row = -1;
      for (int i=k; i < n; ++i)
        if (!IsZero(M2(i,k))) { row = i; break; }
      if (row == -1) return zero(R);
      if (row != k) { M2->mySwapRows(k, row); sign = -sign; }
      // Now use pivot row to reduce all lower rows
      for (int i=k+1; i < n; ++i)
        for (int j=k+1; j < n; ++j)
        {
          SetEntry(M2,i,j,(M2(i,j)*M2(k,k)-M2(i,k)*M2(k,j))/d );
        }
      d = M2(k,k);
    }
    if (sign == 1)
      return M2(n-1,n-1);
    return -M2(n-1,n-1);
  }


//////////////////////////////////////////////////////////////////

  matrix AdjByInverse(ConstMatrixView M)
  {
    // This code works only if the matrix is invertible!
    // It is morally equivalent to swap(lhs, inverse(M)*det(M));
    CoCoA_ASSERT(IsIntegralDomain(RingOf(M))); // needed for this method, not mathematically necessary BUG BUG!!!
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    if (IsField(RingOf(M)))
    {
      const RingElem d = det(M);
      if (IsZero(d)) CoCoA_ERROR(ERR::DivByZero, "AdjByInverse");
      return d*inverse(M);
    }
    FractionField K(NewFractionField(RingOf(M)));
    RingHom R2K(EmbeddingHom(K));
    const matrix ans_K = adj(apply(R2K, M));
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix ans = NewDenseMat(RingOf(M), Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i,j, num(ans_K(i,j)));
    return ans;
  }


  // Simple (but probably not very fast).
  matrix AdjByDetOfMinors(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    matrix adj = NewDenseMat(RingOf(M), n,n);
    vector<long> rows(n-1);
    for (long i=0; i < n-1; ++i) { rows[i] = i+1; }

    vector<long> cols(n-1);
    for (long i=0; i < n; ++i)
    {
      if (i > 0) rows[i-1] = i-1;
      for (long j=0; j < n-1; ++j) { cols[j] = j+1; }
      for (long j=0; j < n; ++j)
      {
        if (j > 0) cols[j-1] = j-1;
        SetEntry(adj, j,i, power(-1, i+j)*det(submat(M,rows,cols)));
      }
    }
    return adj;
  }



  matrix AdjDirect(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    const ring R = RingOf(M);
    matrix adj = NewDenseMat(R, n,n);

    for (int col = 0; col < n; ++col)
    {
      for (int row = 0; row < n; ++row)
      {
        vector<int> rows(n-1);
        for (long i=0; i < n-1; ++i) { rows[i] = i + (i >= row); }

        RingElem entry(R);

        do {
          RingElem term = one(R);
          if ((row+col)%2 == 1) term = -term;
          for (int i=0; i < n-1; ++i)
            term *= M(rows[i], i+(i>=col));
          entry += signature(rows)*term;
        } while ( std::next_permutation(rows.begin(),rows.end()) );

        SetEntry(adj, col,row, entry); // transposed!!
      }
    }
    return adj;
  }

  bool IsZero(const ConstMatrixView& M)
  { return M == ZeroMat(RingOf(M), NumRows(M), NumCols(M)); }


  bool IsZeroRow(const ConstMatrixView& M, long i)
  {
    M->myCheckRowIndex(i, "IsZeroRow(M)");
    return M->myIsZeroRow(i);
  }


  bool IsZeroCol(const ConstMatrixView& M, long j)
  {
    M->myCheckColIndex(j, "IsZeroCol(M)");
    return M->myIsZeroCol(j);
  }


  bool IsSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsSymmetric");
    return M->IamSymmetric();
  }


  bool IsAntiSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsAntiSymmetric");
    return M->IamAntiSymmetric();
  }


  bool IsDiagonal(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsDiagonal");
    return M->IamDiagonal();
  }

  
  bool IsUpperTriangular(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsUpperTriangular");
    return M->IamUpperTriangular();
  }

  
  bool IsLowerTriangular(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsLowerTriangular");
    return M->IamLowerTriangular();
  }

  
  bool IsMat0x0(const ConstMatrixView& M)
  {
    return NumRows(M) == 0 && NumCols(M) == 0;
  }


  bool HasNegEntry(const ConstMatrixView& M)
  {
    return M->IhaveNegEntry();
  }
  

  namespace // anonymous
  {

    BigInt FindRowAndColScales(vector<BigInt>& Rscale, vector<BigInt>& Cscale, const ConstMatrixView& M)
    {
      VerboseLog VERBOSE("FindRowAndColScales");

      CoCoA_ASSERT(IsQQ(RingOf(M)));
      CoCoA_ASSERT(NumCols(M) == NumRows(M));
      const int n = NumRows(M);
      vector< vector<BigInt> > denom(n, vector<BigInt>(n));
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          denom[i][j] = ConvertTo<BigInt>(den(M(i,j)));

//    VERBOSE(25) << "denom: " << denom << std::endl;
      vector<BigInt> RowScale(n, BigInt(1));
      vector<BigInt> ColScale(n, BigInt(1));

      // Spit notionally into k=sqrt(n) blocks;
      // compute gcd of each block, then take largest lcm of pairs of blocks
      const int sqrtn = 1+FloorSqrt(n);
      vector<BigInt> BlockGCD(sqrtn);
      // Compute preliminary row factors
      for (int i=0; i < n; ++i) // row indexer
      {
        int end=n;
        for (int k=0; k < sqrtn; ++k)
        {
          BigInt g;
          const int BlockSize = end/(sqrtn-k); // integer division!
          for (int j=end-BlockSize; j < end; ++j)
            g = gcd(g, denom[i][j]);
          BlockGCD[k] = g;
          end -= BlockSize;
        }
//      VERBOSE(25) << "BlockGCD: " << BlockGCD << std::endl;
        BigInt scale = BlockGCD[0];
        for (int k=1; k < sqrtn; ++k)
          scale = lcm(scale, BlockGCD[k]);
        RowScale[i] = scale;
        if (scale != 1)
          for (int j=0; j < n; ++j)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
      if (VerbosityLevel() >= 90)
      {
        BigInt Stage1Fac(1);
        for (int i=0; i < n; ++i)
          Stage1Fac *= RowScale[i];
        VERBOSE(90) << "Stage1 RowFac = " << FloatStr(Stage1Fac) << std::endl;
      }
      if (VerbosityLevel() >= 91)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(91) << "Stage 1 LogRowScale: " << LogRowScale << std::endl;
      }

      // Compute preliminary col factors
      for (int j=0; j < n; ++j) // col indexer
      {
        int end=n;
        for (int k=0; k < sqrtn; ++k)
        {
          BigInt g;
          const int BlockSize = end/(sqrtn-k); // integer division!
          for (int i=end-BlockSize; i < end; ++i)
            g = gcd(g, denom[i][j]);
          BlockGCD[k] = g;
          end -= BlockSize;
        }
//      VERBOSE(25) << "BlockGCD: " << BlockGCD << std::endl;
        BigInt scale = BlockGCD[0];
        for (int k=1; k < sqrtn; ++k)
          scale = lcm(scale, BlockGCD[k]);
        ColScale[j] = scale;
        if (scale != 1)
          for (int i=0; i < n; ++i)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
      if (VerbosityLevel() >= 90)
      {
        BigInt Stage1Fac(1);
        for (int j=0; j < n; ++j)
          Stage1Fac *= ColScale[j];
        VERBOSE(90) << "Stage1 ColFac = " << FloatStr(Stage1Fac) << std::endl;
      }
      if (VerbosityLevel() >= 91)
      {
        vector<long> LogColScale(n);
        for (int j=0; j < n; ++j)
          LogColScale[j] = FloorLog2(ColScale[j]);
        VERBOSE(91) << "Stage 1 LogColScale: " << LogColScale << std::endl;
      }

//    VERBOSE(25) << "denom: " << denom << std::endl;
//    VERBOSE(25) << "RowScale: " << RowScale << std::endl;

#if 0
      // First idea: gcd(lcm(first_half), lcm(second_half))
      // Works OK, but not great
      const int n2 = n/2; // integer division!
      // Fill RowScale -- phase 1
      for (int i=0; i < n; ++i)
      {
        BigInt block1(1);
        for (int j=0; j < n2; ++j)
          block1 = lcm(block1, denom[i][j]);
        BigInt block2(1);
        for (int j=n2; j < n; ++j)
          block2 = lcm(block2, denom[i][j]);
        const BigInt scale = gcd(block1, block2);
        RowScale[i] = scale;
        if (scale != 1)
          for (int j=0; j < n; ++j)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
    
      if (VerbosityLevel() >= 90)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(90) << "Stage 1 LogRowScale: " << LogRowScale << std::endl;
      }

      // Fill ColScale -- phase 1
      for (int j=0; j < n; ++j)
      {
        BigInt block1(1);
        for (int i=0; i < n2; ++i)
          block1 = lcm(block1, denom[i][j]);
        BigInt block2(1);
        for (int i=n2; i < n; ++i)
          block2 = lcm(block2, denom[i][j]);
        const BigInt scale = gcd(block1, block2);
        ColScale[j] = scale;
        if (scale != 1)
          for (int i=0; i < n; ++i)
            denom[i][j] /=  gcd(scale, denom[i][j]);
      }

      if (VerbosityLevel() >= 90)
      {
        vector<long> LogColScale(n);
        for (int i=0; i < n; ++i)
          LogColScale[i] = FloorLog2(ColScale[i]);
        VERBOSE(90) << "Stage 1 LogColScale: " << LogColScale << std::endl;
      }
#endif

// #if 0
//       THESE TWO LOOPS ARE PROBABLY NOT A GOOD IDEA!!!
//         // First reduce by rows, then by cols
//         // Uh oh: what is the complexity of this loop???
//         for (int i=0; i < n; ++i)
//       {
//         for (int j1=0; j1 < n; ++j1)
//         {
//           if (IsOne(denom[i][j1])) continue;
//           for (int j2=j1+1; j2 < n; ++j2)
//           {
//             const BigInt g = gcd(denom[i][j1], denom[i][j2]);
//             if (IsOne(g)) continue;
//             RowScale[i] *= g;
//             for (int j3=0; j3 < n; ++j3)
//               denom[i][j3] = denom[i][j3]/gcd(g, denom[i][j3]);
//           }
//         }
//       }
//       VERBOSE(30) << "RowScale: " << RowScale << std::endl;

//       // Now columns
//       for (int j=0; j < n; ++j)
//       {
//         for (int i1=0; i1 < n; ++i1)
//         {
//           if (IsOne(denom[i1][j])) continue;
//           for (int i2=i1+1; i2 < n; ++i2)
//           {
//             const BigInt g = gcd(denom[i1][j], denom[i2][j]);
//             if (IsOne(g)) continue;
//             ColScale[j] *= g;
//             for (int i3=0; i3 < n; ++i3)
//               denom[i3][j] = denom[i3][j]/gcd(g, denom[i3][j]);
//           }
//         }
//       }
//       VERBOSE(30) << "ColScale: " << ColScale << std::endl;
//       END OF NOT A GOOD IDEA
// #endif
    
      // Compute row-by-row LCMs
      vector<BigInt> LCMrow(n);
      BigInt ByRowsFac(1);
      for (int i=0; i < n; ++i)
      {
        BigInt RowFac(1);
        for (int j=0; j < n; ++j)
          RowFac = lcm(RowFac, denom[i][j]);
        LCMrow[i] = RowFac;
        ByRowsFac *= RowFac;
      }

      // Compute col-by-col LCMs
      vector<BigInt> LCMcol(n);
      BigInt ByColsFac(1);
      for (int j=0; j < n; ++j)
      {
        BigInt ColFac(1);
        for (int i=0; i < n; ++i)
          ColFac = lcm(ColFac, denom[i][j]);
        LCMcol[j] = ColFac;
        ByColsFac *= ColFac;
      }
      VERBOSE(85) << "ByRowsFac = " << FloatStr(ByRowsFac) << "   ByColsFac = " << FloatStr(ByColsFac) << std::endl;

      // Now pick whichever is smaller from row-by-row or col-by-col
      if (ByRowsFac <= ByColsFac)
      {
        for (int i=0; i < n; ++i)
          RowScale[i] *= LCMrow[i];
      }
      else
      {
        for (int j=0; j < n; ++j)
          ColScale[j] *= LCMcol[j];
      }


      BigInt OverallFactor(1);
      for (int i=0; i <n; ++i) OverallFactor *= RowScale[i];
      for (int j=0; j <n; ++j) OverallFactor *= ColScale[j];

      if (VerbosityLevel() >= 90)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(90) << "Final LogRowScale: " << LogRowScale << std::endl;
        vector<long> LogColScale(n);
        for (int j=0; j < n; ++j)
          LogColScale[j] = FloorLog2(ColScale[j]);
        VERBOSE(90) << "Final LogColScale: " << LogColScale << std::endl;
      }
      VERBOSE(80) << "Overall factor " << FloatStr(OverallFactor) << std::endl;

      swap(Rscale, RowScale);  // really assignment
      swap(Cscale, ColScale);  // really assignment
      return OverallFactor;
    }


    matrix ScaleRowsAndCols(const ConstMatrixView& M, const vector<BigInt>& RowScale, const vector<BigInt>& ColScale)
    {
      const int n = NumRows(M);
      matrix ans = NewDenseMat(RingZZ(), n,n);
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
        {
          const RingElem entry = (RowScale[i]*ColScale[j])*M(i,j);
          CoCoA_ASSERT(IsOne(den(entry)));
          SetEntry(ans,i,j, num(entry));
        }

      return ans;
    }

  } // end of namespace anonymous
  
  RingElem DetOverQQ(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetOverQQ");
    CoCoA_ASSERT(IsQQ(RingOf(M)));
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const int n = NumRows(M);
    if (n == 0) return one(RingOf(M));
    if (n == 1) return M(0,0);
    VERBOSE(80) << "UNTRANSPOSED" << std::endl;
    vector<BigInt> Rscale;
    vector<BigInt> Cscale;
    const BigInt OverallFactor = FindRowAndColScales(Rscale, Cscale, M);
    VERBOSE(80) << std::endl;
    VERBOSE(80) << "TRANSPOSED" << std::endl;
    vector<BigInt> trRscale;
    vector<BigInt> trCscale;
    const BigInt trOverallFactor = FindRowAndColScales(trRscale, trCscale, transpose(M));

    if (OverallFactor <= trOverallFactor)
    {
      return RingElem(RingOf(M), BigRat(ConvertTo<BigInt>(det(ScaleRowsAndCols(M, Rscale, Cscale))), OverallFactor));
    }
    else
      return RingElem(RingOf(M), BigRat(ConvertTo<BigInt>(det(ScaleRowsAndCols(M, trCscale, trRscale))), trOverallFactor));
      
    // BigInt DetZZ = ConvertTo<BigInt>(det(MoverZZ));
    // BigInt CombinedDenom(1);
    // for (int i=0; i < n; ++i)
    //   CombinedDenom *= Rscale[i]*Cscale[i];
    // return RingElem(RingOf(M), BigRat(DetZZ, CombinedDenom));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps.C,v 1.4 2018/06/26 10:01:10 abbott Exp $
// $Log: MatrixOps.C,v $
// Revision 1.4  2018/06/26 10:01:10  abbott
// Summary: Increased VerbosityLevel needed to print blurb
//
// Revision 1.3  2018/06/12 13:55:33  abbott
// Summary: Some cleaning
//
// Revision 1.2  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2018/05/17 15:25:53  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.24  2018/04/18 14:27:15  abbott
// Summary: Minor cleaning
//
// Revision 1.23  2018/03/02 13:44:11  abbott
// Summary: Major revision to conversion mat(QQ) --> mat(ZZ) for determinant
//
// Revision 1.22  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.21  2018/02/27 10:58:59  abbott
// Summary: Added include NumTheory_prime; added first impl for converting det over QQ to det over ZZ
//
// Revision 1.20  2018/02/15 17:03:53  abbott
// Summary: Removed some verbosity
//
// Revision 1.19  2018/02/12 14:52:11  abbott
// Summary: Revised det2x2, det3x3, DetByGauss.  Added det4x4, det5x5, and DetByCRT
//
// Revision 1.18  2017/04/07 14:20:43  bigatti
// -- 2 small changes to unused lines to keep compiler quiet
//
// Revision 1.17  2016/10/27 14:05:39  abbott
// Summary: Added code for det of 0x0 and 1x1 matrices; see issue 956
//
// Revision 1.16  2016/06/10 12:01:12  bigatti
// -- just a change of sign to the result of LinKer
//
// Revision 1.15  2015/12/11 15:51:23  bigatti
// -- added IsLowerTriangular
//
// Revision 1.14  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.13  2015/06/26 14:59:02  abbott
// Summary: Corrected assertion inside InverseByGauss
// Author: JAA
//
// Revision 1.12  2015/05/20 13:27:18  abbott
// Summary: adj now calls AdjByDetOfMinors if det is zero
// Author: JAA
//
// Revision 1.11  2015/05/19 07:26:36  abbott
// Summary: Minor improvement to DetDirect
// Author: JAA
//
// Revision 1.10  2015/05/15 15:59:04  bigatti
// -- fixed swapped arguments in DetDirect
//
// Revision 1.9  2015/05/11 15:49:03  bigatti
// -- now InverseByGauss works also if matrix is over IntegralDomain (and invertible)
//
// Revision 1.8  2015/04/27 10:08:48  bigatti
// Summary: changed myAddMul --> myAddMulLM
//
// Revision 1.7  2015/04/13 16:16:12  abbott
// Summary: Changed "rank" --> "rk"; adjoint" --> "adj"; added "AdjDirect", "DetDirect"
// Author: JAA
//
// Revision 1.6  2014/08/26 12:55:58  abbott
// Summary: Cleaned up DetByGauss; added DetByBareiss
// Author: JAA
//
// Revision 1.5  2014/07/30 14:06:24  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.4  2014/07/08 08:35:16  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.3  2014/07/07 12:23:22  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.2  2014/04/17 13:38:47  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.1  2014/04/11 15:42:37  abbott
// Summary: Renamed from MatrixArith
// Author: JAA
//
// Revision 1.46  2014/04/08 15:40:50  abbott
// Summary: Replaced test for IsZero by IsZeroDivisor
// Author: JAA
//
// Revision 1.45  2014/01/16 16:10:56  abbott
// Removed a blank line.
//
// Revision 1.44  2012/11/23 17:28:26  abbott
// Modified LinSolve so that it returns a 0x0 matrix when no soln exists
// (previously it threw an exception).
//
// Revision 1.43  2012/10/16 09:45:44  abbott
// Replaced RefRingElem by RingElem&.
//
// Revision 1.42  2012/10/03 15:25:05  abbott
// Replaced swap by assignment in DetByGauss; new impl of swap did not work
// in that instance (since one value was a temporary).
//
// Revision 1.41  2012/07/10 12:59:34  bigatti
// -- added two lines to keep compiler quiet
//
// Revision 1.40  2012/07/10 09:48:41  bigatti
// -- fixes of some naive errors
//
// Revision 1.39  2012/07/10 09:23:30  bigatti
// -- separated gauss code from LinSolveByGauss
// -- added LinKerByGauss
//
// Revision 1.38  2012/06/19 15:43:27  abbott
// Added division of a matrix by a scalar.
//
// Revision 1.37  2012/06/11 08:20:33  abbott
// Added multiplication on the right by a scalar.
//
// Revision 1.36  2012/06/10 22:57:31  abbott
// Added negation for matrices -- same as doing (-1)*M.
//
// Revision 1.35  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.34  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.33  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.32  2012/05/04 15:39:29  abbott
// Corrected a comment.
//
// Revision 1.31  2012/04/27 14:49:33  abbott
// Added LinSolve family (incl. LinSolveByGauss, LinSolveByHNF, LinSolveByModuleRepr).
//
// Revision 1.30  2012/04/16 09:21:19  abbott
// Added missing include directive.
//
// Revision 1.29  2012/04/13 16:24:35  abbott
// Added solve and SolveByGauss.
//
// Revision 1.28  2012/04/11 14:03:24  abbott
// Very minor change: slight improvement to readability.
//
// Revision 1.27  2012/03/16 14:42:46  bigatti
// -- fixed AdjointByDetOfMinors
// -- fixed adjoint (field + det(M)=0)
//
// Revision 1.26  2012/02/10 10:26:39  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.25  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.24  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.23  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.22  2011/05/13 16:47:20  abbott
// Added power fn for matrices: partial impl, cannot yet handle negative powers.
//
// Revision 1.21  2011/05/03 13:48:02  abbott
// Now using CanonicalHom inside DetByGauss.
// Cleaner and avoids a mysterious compiler warning.
//
// Revision 1.20  2011/03/22 16:44:19  bigatti
// -- fixed check in det
//
// Revision 1.19  2011/03/16 15:41:06  bigatti
// -- minor cleaning
//
// Revision 1.18  2011/03/10 11:25:46  bigatti
// -- now using long instead of size_t, and len(v) instead of v.size()
//
// Revision 1.17  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.16  2011/03/01 14:13:24  bigatti
// -- added f*M
//
// Revision 1.15  2011/02/28 14:08:49  bigatti
// -- added det3x3
// -- using apply mapping matrix (in DetByGauss)
//
// Revision 1.14  2011/02/10 15:27:06  bigatti
// -- commented #include vector  (included in MatrixArith.H)
//
// Revision 1.13  2011/02/09 16:48:27  bigatti
// -- added + and - for matrices
//
// Revision 1.12  2009/12/23 18:55:16  abbott
// Removed some useless comments.
//
// Revision 1.11  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.10  2009/06/25 16:59:42  abbott
// Minor improvement to some error messages (better coherence & comprehensibility).
//
// Revision 1.9  2008/07/09 16:09:11  abbott
// Removed pointless bogus function declarations.
//
// Revision 1.8  2008/04/22 14:42:03  abbott
// Two very minor changes.
//
// Revision 1.7  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.6  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.5  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.4  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.3  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/30 15:54:15  bigatti
// -- fixed index too big in RankByGauss(ConstMatrix M)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.8  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.5  2006/11/27 16:18:32  cocoa
// -- moved classes declarations from .H to .C (DenseMatrix, DiagMatrix,
//    FieldIdeal, SpecialMatrix)
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/08/17 09:39:07  cocoa
// -- added: elimination ordering matrix for non-homogeneous input
//
// Revision 1.2  2006/07/17 16:58:05  cocoa
// -- added: NewMatrixElim(size_t NumIndets, std::vector<size_t> IndetsToElim)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/02 14:39:20  cocoa
// -- Changed "not" into "!" becuase of M$Windoze (by M.Abshoff)
//
// Revision 1.7  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.6  2006/04/10 13:20:43  cocoa
// -- fixed buglets for Elimination orderings
//
// Revision 1.5  2006/04/05 16:45:29  cocoa
// -- added comment and CoCoA_ASSERT in NewPositiveMatrix
// -- added IsPositiveGrading
//
// Revision 1.4  2006/04/05 14:49:20  cocoa
// -- fixed: NewPositiveMatrix (tested and used in OrdvArith.C)
//
// Revision 1.3  2006/01/19 15:48:49  cocoa
// -- fixed RankByGauss by Stefan Kaspar
//
// Revision 1.2  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.1  2005/03/11 18:38:32  cocoa
// -- first import
//
