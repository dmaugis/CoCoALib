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


#include "CoCoA/LogStream.H"

#include <iostream>

namespace CoCoA
{

  std::ostream* LogStreamForThisBlock::ourActiveLogStreamPtr = &std::cout;

  
  std::ostream& operator<<(std::ostream& out, const LogStreamForThisBlock& /*dummy*/)
  {
    return out << "LogStreamForThisBlock";
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/LogStream.C,v 1.2 2017/04/05 14:23:47 abbott Exp $
// $Log: LogStream.C,v $
// Revision 1.2  2017/04/05 14:23:47  abbott
// Summary: Revised impl of LogStream
//
// Revision 1.1  2017/01/25 13:03:48  abbott
// Summary: Added new global CoCoA::LogStream
//
//
