#ifndef CoCoA_GCDFreeBasis_H
#define CoCoA_GCDFreeBasis_H

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

#include "CoCoA/ring.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  class RingElem;  // fwd decl -- defined in ring.H

  class GCDFreeBasis_RingElem
  {
  public:
    GCDFreeBasis_RingElem() {}
    void myAddInfo(const RingElem& X);
    void myAddInfo(const std::vector<RingElem>& v);
    friend const std::vector<RingElem>& FactorBase(const GCDFreeBasis_RingElem& GFB);
  private: // implementation detail
    struct LCR
    {
      LCR(const RingElem& L, const std::vector<RingElem>& C, const RingElem& R):
          myL(L), myR(R), myC(C) {}
    public: // data members
      RingElem myL;
      RingElem myR;
      std::vector<RingElem> myC;
    };
  private: // data members
    std::vector<RingElem> myCoprimeBasis;
    void myRefineBasis(RingElem X);
    struct LCR myLCR(RingElem A, RingElem B) const;
  };


  std::ostream& operator<<(std::ostream& out, const GCDFreeBasis_RingElem& GFB);


  inline const std::vector<RingElem>& FactorBase(const GCDFreeBasis_RingElem& GFB)
  { return GFB.myCoprimeBasis; }



/////  std::vector<RingElem> GCDFreeBasis(const std::vector<RingElem>& L);
/////  RingElem RefineGCDFreeBasis(std::vector<RingElem>& basis, RingElem N);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/GCDFreeBasis.H,v 1.2 2018/06/25 12:32:10 abbott Exp $
// $Log: GCDFreeBasis.H,v $
// Revision 1.2  2018/06/25 12:32:10  abbott
// Summary: Restructured impl
//
// Revision 1.1  2017/02/01 10:37:35  abbott
// Summary: GCDFreeBasis fns
//
//

#endif
