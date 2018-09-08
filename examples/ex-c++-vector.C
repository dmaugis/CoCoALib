// Copyright (c) 2017  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic use of C++ vectors. \n";

const string LongDescription =
  "This is an example showing some basic use of C++ vectors. \n"
  "C++ has few built-in types, but the STL (Standard Template\n"
  "Library) contains many useful extensions.  std::vector is \n"
  "is one of these extensions.  It is a homogeneous array: a \n"
  "collection of many objects of the same type, and each one \n"
  "may be accessed directly by an integer index.             \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // An array of 10 "long" integers; initially all 0.
    vector<long> v(10);
    int n = len(v);  // number of entries in v
    // CAREFUL!  Indices start at 0 and go to n-1
    for (int i=0; i < n; ++i)
    {
      v[i] = i; // index goes inside SQUARE brackets
    }
    cout << "v = " << v << endl; // CoCoALib can print a vector

    // Vector of strings; initially empty.
    vector<string> mesg;
    mesg.push_back("Hello"); // append new value at the end
    mesg.push_back("World"); // (ditto)
    cout << mesg << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-vector.C,v 1.1 2017/02/15 12:22:42 abbott Exp $
// $Log: ex-c++-vector.C,v $
// Revision 1.1  2017/02/15 12:22:42  abbott
// Summary: New C++ examples
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
