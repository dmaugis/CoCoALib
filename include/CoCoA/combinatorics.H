#ifndef CoCoA_combinatorics_H
#define CoCoA_combinatorics_H

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
//   -------------------------------------------------------------------------


#include "CoCoA/MachineInt.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  std::vector<long> RandomSubsetIndices(const MachineInt& n);
  std::vector<long> RandomSubsetIndices(const MachineInt& n, const MachineInt& r);

  std::vector<long> RandomTupleIndices(const MachineInt& n, const MachineInt& r);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/combinatorics.H,v 1.1 2015/05/20 13:38:32 abbott Exp $
// $Log: combinatorics.H,v $
// Revision 1.1  2015/05/20 13:38:32  abbott
// Summary: New files for combinatorial fns: RandomSubset, RandomTuple
// Author: JAA
//
//

#endif
