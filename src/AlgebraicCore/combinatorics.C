//   Copyright (c)  2015  John Abbott

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


#include "CoCoA/combinatorics.H"
#include "CoCoA/error.H"
#include "CoCoA/random.H"

#include <algorithm>
using std::sort;
//#include <vector>
using std::vector;

namespace CoCoA
{

  std::vector<long> RandomSubsetIndices(const MachineInt& n)
  {
    if (IsNegative(n)) CoCoA_ERROR(ERR::NotNonNegative, "RandomSubsetIndices");
    const long N = AsSignedLong(n);
    vector<long> ans;
    for (int i=0; i < N; ++i)
      if (RandomBool())
        ans.push_back(i);
    return ans;
  }


  // Adapted from Wikipedia entry "Reservoir sorting"
  std::vector<long> RandomSubsetIndices(const MachineInt& n, const MachineInt& r)
  {
    if (IsNegative(n) || IsNegative(r)) CoCoA_ERROR(ERR::NotNonNegative, "RandomSubsetIndices");
    const long N = AsSignedLong(n);
    const long R = AsSignedLong(r);
    if (R > N) CoCoA_ERROR(ERR::BadArg, "RandomSubsetIndices");
    vector<long> ans(R);
    for (long i=0; i < R; ++i)
      ans[i] = i;
    for (long k=R; k < N; ++k)
    {
      const long i = RandomLong(0,k); // both extremes included
      if (i < R) ans[i] = k;
    }
    sort(ans.begin(), ans.end());
    return ans;
  }


  std::vector<long> RandomTupleIndices(const MachineInt& n, const MachineInt& r)
  {
    if (IsNegative(n) || IsNegative(r)) CoCoA_ERROR(ERR::NotNonNegative, "RandomSubsetIndices");
    const long N = AsSignedLong(n);
    const long R = AsSignedLong(r);
    vector<long> ans(R);
    for (long i=0; i < R; ++i)
      ans[i] = RandomLong(0,N-1); // both extremes included
    return ans;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/combinatorics.C,v 1.1 2015/05/20 13:38:32 abbott Exp $
// $Log: combinatorics.C,v $
// Revision 1.1  2015/05/20 13:38:32  abbott
// Summary: New files for combinatorial fns: RandomSubset, RandomTuple
// Author: JAA
//
//
