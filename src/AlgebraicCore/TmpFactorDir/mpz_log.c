//   Copyright (c)  1997-2006,2010  John Abbott

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


#include "mpz_log.h"
#include <math.h>

/* It would be good if GMP offered this fn. */

/* It returns a close approximation to log(abs(x)) */
/* BUG: ungraceful if x=0 */
double mpz_log(const mpz_t x)
{
  /*static*/ const double log2 = log(2.0); // should be static, but not allowed in C++03 according to clang :-(
  long exponent;
  double mantissa = mpz_get_d_2exp(&exponent, x); // GMP BUG: exponent could overflow
  return exponent*log2 + log(fabs(mantissa));
}


double mpq_log(const mpq_t q)
{
  return mpz_log(mpq_numref(q)) - mpz_log(mpq_denref(q));
}
