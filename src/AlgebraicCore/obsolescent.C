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


#include "CoCoA/library.H"

namespace CoCoA
{

  namespace // anonymous for file local fn
  {
    // Procedure to give error if obsolescent fns are forbidden, and otherwise print out a warning on CoCoA::LogStream.
    void LogObsolescentFn(const char* const FnName, const char* const UsefulAdvice)
    {
      if (!IsAllowedObsolescentFnCall())
        CoCoA_ERROR(ERR::OBSOLESCENT, FnName + string(" -- ") + UsefulAdvice);
      LogStream() << "WARNING: called obsolescent fn `" << FnName << "' -- " << UsefulAdvice << endl;
    }

  } // end of namespace anonymous

  
  ///////////////////////////////////////////////////////
  // The obsolescent fns below are ALWAYS defined.
  
  bool IsRadical(ConstRefPPMonoidElem pp)  // RENAMED to IsSqFree
  {
    LogObsolescentFn("IsRadical(ConstRefPPMonoidElem)", "renamed to IsSqFree");
    return IsSqFree(pp);
  }


  bool AreGensSquareFreeMonomial(const ideal& I)
  {
    LogObsolescentFn("AreGensSquareFreeMonomial(ideal)",
                     "renamed to AreGensSqFreeMonomial");
    return AreGensSqFreeMonomial(I);
  }

  
  PPOrdering NewLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewLexOrdering", "use pseudo-ctor `lex'");
    return lex(NumIndets);
  }

  PPOrdering NewStdDegLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewStdDegLexOrdering", "use pseudo-ctor `StdDegLex'");
    return StdDegLex(NumIndets);
  }
  
  PPOrdering NewStdDegRevLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewStdDegRevLexOrdering", "use pseudo-ctor `StdDegRevLex'");
    return StdDegRevLex(NumIndets);
  }

  ideal minimalize(const ideal& I)
  {
    LogObsolescentFn("minimalize", "use \"IdealOfMinGens\"");
    return IdealOfMinGens(I);
  }

  FGModule minimalize(const FGModule& M)
  {
    LogObsolescentFn("minimalize", "use \"SubmoduleOfMinGens\"");
    return SubmoduleOfMinGens(M);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/obsolescent.C,v 1.8 2017/11/20 20:38:27 bigatti Exp $
// $Log: obsolescent.C,v $
// Revision 1.8  2017/11/20 20:38:27  bigatti
// -- added minimalized
//
// Revision 1.7  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.6  2017/03/29 15:40:40  abbott
// Summary: Now prints UsefulAdvice in log message
//
// Revision 1.5  2017/01/25 13:01:44  abbott
// Summary: Warning message now output to CoCoA::LogStream (instead of clog)
//
// Revision 1.4  2016/11/07 14:16:50  bigatti
// -- added AreGensSquareFreeMonomial
//
// Revision 1.3  2016/11/05 16:34:17  abbott
// Summary: Put LogObsolescentFn into anon namespace
//
// Revision 1.2  2016/11/04 20:44:07  abbott
// Summary: Cleaned and simplified
//
// Revision 1.1  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
//
