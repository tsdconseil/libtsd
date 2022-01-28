#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/tests.hpp"
#include <cstdio>

using namespace tsd;
using namespace tsd::fourier;
using namespace tsd::vue;

int test_filtre_fft()
{
  msg("Test filtre FFT...");

  int BS = 512;
  FiltreFFTConfig config;
  config.nb_zeros_min             = 0;
  config.dim_blocs_temporel   = BS;
  config.avec_fenetrage       = true;
  config.traitement_freq = [&](ArrayXcf &X)
  {
  };

  auto [ola, NOLA] = filtre_fft(config);

  int M = 1024, N = 4096;

  ArrayXcf x = ArrayXcf::Zero(N), motif = ArrayXcf::Zero(M);

  for(auto i = 0; i < M; i++)
  {
    /*float fs = 100.0f;
    float Omega = 1.0f / fs;
    x(i) = std::exp(2.0f*π*cfloat(0,1)*((float)i)*Omega);
    x(i) *= ((float) i) / N;*/


    float a = std::abs(((float) i) - ((float) M)/2);
    a = (M/2 - a) / (M/2);
    motif(i) = a - cfloat(0,1) * a;
  }

  motif.head(M/4).setConstant(cfloat(1,0));
  motif.segment(M/4, M/4).setConstant(cfloat(3,0));
  motif.segment(2*M/4, M/4).setConstant(cfloat(2,0));
  motif.segment(3*M/4, M/4).setConstant(cfloat(-3.0,0));

  x.segment(N/2,M) = motif;

  ArrayXcf y;
  ola->step(x, y);


  return 0;
}





