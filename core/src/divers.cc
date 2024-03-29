#include "tsd/divers.hpp"

namespace tsd {


float sinc(float T, float f)
{
  soit a = π_f * T * f;
  si(abs(a) < 1e-7f)
    retourne T;
  retourne sin(a) / (π_f * f);
}

float sinc(float t)
{
  // Fréquence de coupure = 1/2 <-> largeur = 1
  retourne sinc(1, t);
}


double Dirichlet(entier N, double Ω)
{
  soit ε = 1e-5;
  Ω = modulo_2π(Ω);
  si((Ω < ε) || ((2 * π - Ω) < ε))
    retourne 2 * N  + 1;

  si(abs(sin(Ω/2)) < 1e-7)
  {
    msg_avert("Dirichlet : N = {}, omega = {}, pb.", N, Ω);
  }

  retourne sin((N+0.5)*Ω) / sin(Ω/2);
}

Poly<float> Chebychev_TU(entier n, bouléen espèce)
{
  /** @brief polynome simple */
  soit z = Poly<float>::z;
  soit T0 = Poly<float>::one();
  soit T1 = (espèce == 1) ? z : 2.0f * z;

  si(n == 0)
    retourne T0;
  sinon si(n == 1)
    retourne T1;

  pour(auto k = 2; k <= n; k++)
  {
    soit T2 = 2.0f * z * T1 - T0;
    T0 = T1;
    T1 = T2;
  }

  retourne T1;
}

Poly<float> Chebychev_T(entier n)
{
  retourne Chebychev_TU(n, 1);
}

Poly<float> Chebychev_U(entier n)
{
  retourne Chebychev_TU(n, 2);
}

/*
// https://en.wikipedia.org/wiki/Bessel_polynomials
Poly<float> BesselRev(entier n)
{
  soit z  = Poly<float>::z;
  soit T0 = Poly<float>::one();
  soit T1 = z + 1;

  si(n == 0)
    retourne T0;
  sinon si(n == 1)
    retourne T1;

  pour(auto k = 2; k <= n; k++)
  {
    soit T2 = (2*n-1) * T1 + z * z * T0;
    T0 = T1;
    T1 = T2;
  }

  retourne T1;
}*/

}
