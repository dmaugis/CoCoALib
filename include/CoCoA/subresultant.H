#ifndef CoCoA_subresultant_H
#define CoCoA_subresultant_H

//   Copyright (c)  2017  John Abbott

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

#include "CoCoA/ring.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

//  class RingElem;  // fwd decl -- defined in ring.H
//  class ConstRefRingElem;  // fwd decl -- defined in ring.H

  RingElem resultant(ConstRefRingElem f, ConstRefRingElem g, long x);
  std::vector<RingElem> SubresultantSeq(ConstRefRingElem f, ConstRefRingElem g, long x);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/subresultant.H,v 1.2 2017/04/05 14:30:42 abbott Exp $
// $Log: subresultant.H,v $
// Revision 1.2  2017/04/05 14:30:42  abbott
// Summary: Correct fn name
//
// Revision 1.1  2017/02/20 12:08:53  abbott
// Summary: Added new files MatrixFp.C and subresultant.C
//
//

#endif
