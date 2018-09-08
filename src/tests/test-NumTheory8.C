// Copyright (c) 2018  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/GlobalManager.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/BuildInfo.H"

#include "CoCoA/NumTheory-prime.H"

#include <iostream>
using namespace std;

namespace CoCoA
{

  // Test the various "prime sequence" implementations

  void program()
  {
    GlobalManager CoCoAFoundations;

    // FastFinitePrimeSeq is really just table-lookup.
    // Check the sequence is coherent: starts with 2, and each
    // successive entry is NextPrime of previous entry.
    FastFinitePrimeSeq seq3;
    long p3 = *seq3;
    CoCoA_ASSERT_ALWAYS(p3 == 2);
    while (!IsEnded(++seq3))
    {
      CoCoA_ASSERT_ALWAYS(*seq3 == NextPrime(p3));
      p3 = *seq3;
    }


    // FastMostlyPrimeSeq is same as FastFinitePrimeSeq up to index 82024
    // then it generates some composite numbers (but with no factor < 29)
    FastMostlyPrimeSeq seq0;
    for (int i=0; i < 100000; ++i)
    {
      const long n = *seq0;
      if (n < 25) { CoCoA_ASSERT_ALWAYS(IsPrime(n)); ++seq0; continue; }
      CoCoA_ASSERT_ALWAYS(n%2 != 0);
      CoCoA_ASSERT_ALWAYS(n%3 != 0);
      CoCoA_ASSERT_ALWAYS(n%5 != 0);
      CoCoA_ASSERT_ALWAYS(n%7 != 0);
      CoCoA_ASSERT_ALWAYS(n%11 != 0);
      CoCoA_ASSERT_ALWAYS(n%13 != 0);
      CoCoA_ASSERT_ALWAYS(n%17 != 0);
      CoCoA_ASSERT_ALWAYS(n%19 != 0);
      CoCoA_ASSERT_ALWAYS(n%23 != 0);
      ++seq0;
    }


    // PrimeSeq generates the sequence of primes starting from 2.
    // Check the sequence is coherent: starts with 2, and each
    // successive entry is NextPrime of previous entry.
    // Internally it "changes gear" at around 82024; check for
    // correct behaviour around this point.
    PrimeSeq seq1;
    CoCoA_ASSERT_ALWAYS(*seq1 == 2);
    for (int i=0; i < 82000; ++i)
      ++seq1;
    long p1 = *seq1;
    for (int i=0; i < 100; ++i)
    {
      ++seq1;
      CoCoA_ASSERT_ALWAYS(*seq1 == NextPrime(p1));
      p1 = *seq1;
    }

    // PrimeSeqForCRT generates "large" small primes.
    // It uses a table for the first 50000ish entries then
    // generates on the fly.  We check for a smooth transition.
    PrimeSeqForCRT seq2;
    for (int i=0; i < 50000; ++i)
      ++seq2;
    long p2 = *seq2;
    for (int i=0; i < 100; ++i)
    {
      ++seq2;
      CoCoA_ASSERT_ALWAYS(*seq2 == NextPrime(p2));
      p2 = *seq2;
    }


  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-NumTheory8.C,v 1.2 2018/03/09 14:34:34 abbott Exp $
// $Log: test-NumTheory8.C,v $
// Revision 1.2  2018/03/09 14:34:34  abbott
// Summary: Corrected minor bug
//
// Revision 1.1  2018/03/02 13:46:29  abbott
// Summary: Test for prime generation (via seqs)
//
//
