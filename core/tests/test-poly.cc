// Tests polynomes, fractions rationnelles, polynomes de Laurent

#include "tsd/filtrage/frat.hpp"

using namespace tsd;
using namespace tsd::filtrage;

entier test_poly()
{
  soit z = FRat<float>::z();


  soit F2 = z * 0.5;
  msg("F2 = {}", F2);

  soit F3 = z - 0.5;
  msg("F3 = {}", F3);

  soit F1 = (z /** 0.5*/) / (z);//  - 0.5);


  msg("F1 = {}", F1);

  tsd_assert(F1.horner(1.0f) == 1.0f);



  retourne 0;
}

