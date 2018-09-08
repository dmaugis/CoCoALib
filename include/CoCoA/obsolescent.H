#ifndef CoCoA_obsolescent_H
#define CoCoA_obsolescent_H

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
//   -------------------------------------------------------------------------


// ***** !!! IMPORTANT !!! *****
// This file is not included automatically in CoCoA/library.H;
// you must use a separate #include "CoCoA/obsolescent.H" directive.

// ***** OBSOLESCENT FUNCTIONS *****

namespace CoCoA
{

  class ideal;                // fwd decl -- defined in ideal.H
  class ConstRefPPMonoidElem; // fwd decl -- defined in PPMonoid.H
  class MachineInt;           // fwd decl -- defined in MachineInt.H
  class PPOrdering;           // fwd decl -- defined in PPOrdering.H

  bool IsRadical(ConstRefPPMonoidElem pp);  // RENAMED to IsSqFree

  bool AreGensSquareFreeMonomial(const ideal& I);  // renamed to AreGensSqFreeMonomial

  PPOrdering NewLexOrdering(const MachineInt& NumIndets);          // use pseudo-ctor lex
  PPOrdering NewStdDegLexOrdering(const MachineInt& NumIndets);    // use pseudo-ctor StdDegLex
  PPOrdering NewStdDegRevLexOrdering(const MachineInt& NumIndets); // use pseudo-ctor StdDegRevLex

  ideal minimalize(const ideal& I);       // use IdealOfMinGens
  FGModule minimalize(const FGModule& M); // use SubmoduleOfMinGens
  

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/obsolescent.H,v 1.5 2017/11/20 20:38:27 bigatti Exp $
// $Log: obsolescent.H,v $
// Revision 1.5  2017/11/20 20:38:27  bigatti
// -- added minimalized
//
// Revision 1.4  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.3  2016/11/07 14:16:51  bigatti
// -- added AreGensSquareFreeMonomial
//
// Revision 1.2  2016/11/04 20:40:47  abbott
// Summary: Cleaning
//
// Revision 1.1  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
//

#endif
