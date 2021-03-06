// Copyright (c) 2017  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some very basic C++:  \n"
  "basic types and how to print.                    \n";

const string LongDescription =
  "This is an example showing some very basic C++:       \n"
  "basic types and how to print.  See also ex-c++-bool.C \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // This is a comment... up to the end of the line --->

    /* This is a comment...
       ... up to the end-of-comment marker --> */

    // some BASIC TYPES in C++
    int a = 1;  // an integer value (of limited range)
    long b = 2; // an integer value (of greater range)
    string str = "a succession of characters";
    double pi = 3.14159; // a floating-point value
    bool flag = true; // boolean value: either "true" or "false"
    
    //-- PRINTING -------------------------------------------------
    // `cout' is the computer screen (usually)
    // the operator `<<' means "send to" or "print on"
    // `endl' means to end the current line and start a new one

    cout << "Hello, World!" << endl;

    cout << "The string variable `str' has value " << str << endl;

    cout << "pi is approximately " << pi << endl;
    cout << endl;  // just produces an empty line
    
    cout << a << "+" << b << " = " << a+b << endl;

    cout << "The boolean variable `flag' has value " << flag << endl;
  }

} // end of namespace CoCoA

// IGNORE THE STUFF BELOW (at least for now)

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-basic.C,v 1.4 2017/02/24 16:58:00 abbott Exp $
// $Log: ex-c++-basic.C,v $
// Revision 1.4  2017/02/24 16:58:00  abbott
// Summary: Removed an empty line
//
// Revision 1.3  2017/02/15 12:21:53  abbott
// Summary: Added bool as a basic type
//
// Revision 1.2  2017/02/10 16:40:35  abbott
// Summary: Minor improvement
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
