//   Copyright (c)  2016  John Abbott

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

#include "CoCoA/VectorOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/utils.H"

#include<list>
using std::list;
#include<vector>
using std::vector;

namespace CoCoA
{

  // This should be a template fn?
  BigInt BinaryProduct(vector<BigInt>& v)
  {
    while (true)
    {
      const int n = len(v);
      if (n < 4) break;
      if (n&1) v[n-2] *= v[n-1];
      const int n2 = n/2; // integer division!
      for (int i=0; i < n2; ++i)
        v[i] = v[2*i]*v[2*i+1];
      v.resize(n2);
    }
    const int n = len(v);
    if (n == 1) return v[0];
    if (n == 2) return v[1]*v[1];
    return v[0]*v[1]*v[2];
  }

  BigInt BinaryProduct(const list<BigInt>& l)
  {
    if (l.empty()) CoCoA_ERROR(ERR::Empty, "BinaryProduct");
    if (len(l) == 1) return l.front();
    if (len(l) == 2) return l.front()*l.back();
    list<BigInt> tmp;
    for (list<BigInt>::const_iterator it = l.begin(); it != l.end(); ++it)
    {
      const list<BigInt>::const_iterator prev = it;
      ++it;
      if (it == l.end()) { tmp.back() *= *it; break; }
      tmp.push_back((*prev)*(*it));
    }
    return BinaryProduct(tmp);
  }
  
} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/VectorOps.C,v 1.3 2018/05/22 14:16:40 abbott Exp $
// $Log: VectorOps.C,v $
// Revision 1.3  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.2  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2018/05/17 15:24:30  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.1  2017/02/01 12:17:57  abbott
// Summary: Added new fn BinaryProduct -- never called currently
//
//
