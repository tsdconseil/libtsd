/** @file divers.cc
 *  @brief desc
 */

#include "tsd/divers.hpp"

using namespace std;

namespace tsd {


float sinc(float T, float f)
{
  auto a = π_f * T * f;
  if(abs(a) < 1e-7f)
    return T;
  return sin(a) / (π_f * f);
}

float sinc(float t)
{
  // Fréquence de coupure = 1/2 <-> largeur = 1
  return sinc(1, t);
}


/*float sinc2(float t, float fc)
{
  return sinc(2 * fc, t);
}*/


double Dirichlet(int N, double Ω)
{
  auto ε = 1e-5;
  Ω = modulo_2π(Ω);
  if((Ω < ε) || ((2 * π - Ω) < ε))
    return 2 * N  + 1;

  if(abs(sin(Ω/2)) < 1e-7)
  {
    msg_avert("Dirichlet : N = {}, omega = {}, pb.", N, Ω);
  }

  return sin((N+0.5)*Ω) / sin(Ω/2);
}

//   return sin((N+0.5)*x)/sin(x/2);

/** @brief polynome simple */
auto z = Poly<float>::z;

Poly<float> Chebychev_T(int n)
{
  if(n == 0)
    return Poly<float>::one();
  else if(n == 1)
    return z;

  // Attention : pas très fin -> cout quadratique
  return 2.0f * z * Chebychev_T(n-1) - Chebychev_T(n-2);
}

Poly<float> Chebychev_U(int n)
{
  if(n == 0)
    return Poly<float>::one();
  else if(n == 1)
    return 2.0f * z;

  // Attention : pas très fin -> cout quadratique
  return 2.0f * z * Chebychev_U(n-1) - Chebychev_U(n-2);
}



}
