#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"

entier test_filtre_fft()
{
  msg("Test filtre FFT...");

  soit BS = 512;
  FiltreFFTConfig config;
  config.nb_zeros_min         = 0;
  config.dim_blocs_temporel   = BS;
  config.avec_fenetrage       = oui;
  config.traitement_freq = [&](Veccf &X)
  {
    // Ne fait rien
  };

  soit [ola, NOLA] = filtre_fft(config);

  soit M = 1024, N = 4096;

  soit x     = Veccf::zeros(N),
       motif = Veccf::zeros(M);

  pour(auto i = 0; i < M; i++)
  {
    /*soit fs = 100.0f;
    soit Omega = 1.0f / fs;
    x(i) = exp(2.0f*π*cfloat(0,1)*((float)i)*Omega);
    x(i) *= ((float) i) / N;*/


    soit a = abs(((float) i) - ((float) M)/2);
    a = (M/2 - a) / (M/2);
    motif(i) = a - cfloat(0,1) * a;
  }

  motif.head(M/4).setConstant(cfloat(1,0));
  motif.segment(M/4, M/4).setConstant(cfloat(3,0));
  motif.segment(2*M/4, M/4).setConstant(cfloat(2,0));
  motif.segment(3*M/4, M/4).setConstant(cfloat(-3.0,0));

  x.segment(N/2,M) = motif;

  soit y = ola->step(x);

  {
    Figures f;
    f.subplot().plot(x, "", "x");
    f.subplot().plot(y, "", "y");
    f.afficher();
  }

  retourne 0;
}





