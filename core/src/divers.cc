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

//   retourne sin((N+0.5)*x)/sin(x/2);

// Implémentation naïve (coût quadratique)
static Poly<float> Chebychev_T_rec(entier n)
{
  /** @brief polynome simple */
  soit z = Poly<float>::z;

  si(n == 0)
    retourne Poly<float>::one();
  sinon si(n == 1)
  retourne z;

  // Attention : pas très fin -> cout quadratique
  retourne 2.0f * z * Chebychev_T(n-1) - Chebychev_T(n-2);
}


Poly<float> Chebychev_T(entier n)
{
  /** @brief polynome simple */
  soit z = Poly<float>::z;
  soit T0 = Poly<float>::one();
  soit T1 = z;

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

Poly<float> Chebychev_U(entier n)
{
  /** @brief polynome simple */
  soit z = Poly<float>::z;

  si(n == 0)
    retourne Poly<float>::one();
  sinon si(n == 1)
    retourne 2.0f * z;

  // Attention : pas très fin -> cout quadratique
  retourne 2.0f * z * Chebychev_U(n-1) - Chebychev_U(n-2);
}



}
