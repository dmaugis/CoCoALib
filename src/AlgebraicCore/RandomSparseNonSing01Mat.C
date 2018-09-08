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

#include "CoCoA/RandomSparseNonSing01Mat.H"

#include "CoCoA/MachineInt.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/random.H"

#include <cmath>
using std::log;

namespace CoCoA
{

  matrix RandomSparseNonSing01Mat(const ring& R, const MachineInt& N)
  {
    if (IsNegative(N) || !IsSignedLong(N))
      CoCoA_ERROR(ERR::NotPositive, "RandomSparseNonSing01Mat");
    const long n = AsSignedLong(N);
    
    const double MagicFactor = 0.80;
    const double prob = MagicFactor*(std::log(n)/n);
    matrix M = NewDenseMat(R, n,n);
    while (true)
    {
      for (int i=0; i < n; ++i)
      {
        bool RowHasNZEntry = false;
        do
        {
          for (int j=0; j < n; ++j)
          {
            if (RandomBiasedBool(prob))
            {
              RowHasNZEntry = true;
              SetEntry(M, i,j, one(R));
            }
          }
        } while (!RowHasNZEntry);
      }
      if (!IsZero(det(M))) return M;
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          if (!IsZero(M(i,j))) SetEntry(M, i,j, zero(R));
    }
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RandomSparseNonSing01Mat.C,v 1.2 2018/05/17 15:40:20 bigatti Exp $
// $Log: RandomSparseNonSing01Mat.C,v $
// Revision 1.2  2018/05/17 15:40:20  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.1  2017/11/14 14:59:35  abbott
// Summary: New fns: RandomSmallPrime, RandomSparseNonSing01Matrix
//
//
