#ifndef CoCoA_NumTheory_GCDFreeBasis_H
#define CoCoA_NumTheory_GCDFreeBasis_H

//   Copyright (c)  2018  John Abbott, Anna M. Bigatti

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


#include "CoCoA/assert.H"
#include "CoCoA/BigInt.H"


#include <iosfwd>
//using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H

  class GCDFreeBasis_BigInt
  {
  public:
    GCDFreeBasis_BigInt() {}
    void myAddInfo(const MachineInt& n);
    void myAddInfo(const BigInt& N);
    void myAddInfo(const std::vector<long>& v);
    void myAddInfo(const std::vector<BigInt>& v);
    friend const std::vector<BigInt>& FactorBase(const GCDFreeBasis_BigInt& GFB);
  private: // implementation detail
    struct LCR
    {
      LCR(const BigInt& L, const std::vector<BigInt>& C, const BigInt& R):
          myL(L), myR(R), myC(C) {}
    public: // data members
      BigInt myL;
      BigInt myR;
      std::vector<BigInt> myC;
    };
  private: // data members
    std::vector<BigInt> myCoprimeBasis;
    void myRefineBasis(BigInt N);
    struct LCR myLCR(BigInt A, BigInt B) const;
  };


  std::ostream& operator<<(std::ostream& out, const GCDFreeBasis_BigInt& GFB);


  inline const std::vector<BigInt>& FactorBase(const GCDFreeBasis_BigInt& GFB)
  { return GFB.myCoprimeBasis; }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-GCDFreeBasis.H,v 1.1 2018/06/25 15:31:34 abbott Exp $
// $Log: NumTheory-GCDFreeBasis.H,v $
// Revision 1.1  2018/06/25 15:31:34  abbott
// Summary: GCDFreeBasis for integers
//
//

#endif
