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


#include "CoCoA/exception.H"

#include <iostream>
#include <sstream>
using std::ostringstream;
//#include <string>
using std::string;

namespace CoCoA
{

  // Must define this because of the virtual fns
  exception::~exception()
  {}
  

  std::ostream& exception::myOutputSelf(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CoCoA::exception(\"" << myMessage << '"';
    if (!myContext.empty()) out << ", context=\"" << myContext << '"';
    out << ")";
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const exception& exc)
  {
    return exc.myOutputSelf(out);
  }

  void PrintInFrame(std::ostream& out, const exception& exc)
  {
    ostringstream buffer;
    buffer << ">>>> " << exc << " <<<<";
    const long n = buffer.str().size();
    const string HorizontalLine(n, '-');
    using std::endl;
    out << endl
        << HorizontalLine << endl
        << buffer.str() << endl
        << HorizontalLine << endl;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/exception.C,v 1.3 2017/07/22 12:59:38 abbott Exp $
// $Log: exception.C,v $
// Revision 1.3  2017/07/22 12:59:38  abbott
// Summary: Added virt mem fn myOutputSelf; added PrintInFrame.
//
// Revision 1.2  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.1  2015/06/26 14:56:25  abbott
// Summary: Created new class "exception"
// Author: JAA
//
//
