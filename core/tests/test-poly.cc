// Tests polynomes, fractions rationnelles, polynomes de Laurent

#include "tsd/filtrage/frat.hpp"

using namespace tsd;
using namespace tsd::filtrage;

void test_poly()
{
  soit z  = FRat<float>::z(),
       F2 = z * 0.5,
       F3 = z - 0.5,
       F1 = (z /** 0.5*/) / (z);//  - 0.5);
  msg("F2 = {}", F2);
  msg("F3 = {}", F3);
  msg("F1 = {}", F1);
  assertion(F1.horner(1.0f) == 1.0f);
}

