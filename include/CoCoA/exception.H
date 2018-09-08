#ifndef CoCoA_exception_H
#define CoCoA_exception_H

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

#include <string>
// using std::string;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class exception
  {
  public:
    explicit exception(const std::string& message): myMessage(message), myContext() {};
    exception(const std::string& message, const std::string& context): myMessage(message), myContext(context) {};
///???    const char* what() const throw() { return myMessage.c_str(); }
    virtual std::ostream& myOutputSelf(std::ostream& out) const;
    virtual ~exception();
    friend const std::string& message(const exception& e);
    friend const std::string& context(const exception& e);
  protected: // data members
    std::string myMessage;
    std::string myContext;
  };


  std::ostream& operator<<(std::ostream& out, const exception& exc);
  void PrintInFrame(std::ostream& out, const exception& exc);

  inline const std::string& message(const exception& e)
  {
    return e.myMessage;
  }

  inline const std::string& context(const exception& e)
  {
    return e.myContext;
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/exception.H,v 1.3 2017/07/22 12:59:37 abbott Exp $
// $Log: exception.H,v $
// Revision 1.3  2017/07/22 12:59:37  abbott
// Summary: Added virt mem fn myOutputSelf; added PrintInFrame.
//
// Revision 1.2  2017/07/15 15:16:33  abbott
// Summary: Removed useless include directive
//
// Revision 1.1  2015/06/26 14:56:25  abbott
// Summary: Created new class "exception"
// Author: JAA
//
//

#endif
