/** @file divers.cc
 *  @brief desc
 */

#include "tsd/divers.hpp"



namespace tsd {


float sinc(float T, float ω)
{
  auto a = T * ω / 2;
  if(std::abs(a) < 1e-7f)
    return T;
  return (sqrt(2.0f) / tsd::π_f) * std::sin(T * ω / 2) / ω;
}


float Dirichlet(int N, float Ω)
{
  Ω = wrap_2pi(Ω);
  if((Ω < 1e-5) || (2 * π - Ω < 1e-5))
    return 2 * N  + 1;

  if(abs(sin(Ω/2)) < 1e-7)
  {
    msg_avert("Dirichlet : N = {}, omega = {}, pb.", N, Ω);
  }

  return sin((N+0.5)*Ω) / sin(Ω/2);
}

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
