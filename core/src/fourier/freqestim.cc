#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"

namespace tsd::fourier {


tuple<Vecf, Vecf> psd_welch(const Veccf &x, entier N, cstring fen)
{
  soit freqs = psd_freqs(N),
       S     = Vecf::zeros(N),
       f     = tsd::filtrage::fenetre(fen, N, non);

  pour(auto i = 0; i + N < x.rows(); i += N/2)
  {
    soit xp = x.segment(i, N) * f;
    S += fftshift(abs2(fft(xp)));
  }

  retourne {freqs, pow2db(S)};
}

Vecf tfd_freqs(entier n, bouléen avec_shift)
{
  si(avec_shift)
  {
    si(est_pair(n))
      retourne linspace(-0.5,0.5-1.0/n,n);
    retourne linspace(-0.5+1.0/n, 0.5, n);
  }

  Vecf f(n);
  si(est_pair(n))
  {
    f.head(n/2) = linspace(0, 0.5-1.0/n, n/2);
    f.tail(n/2) = linspace(-0.5, -1.0/n, n/2);
  }
  sinon
  {
    f.head(n/2+1) = linspace(0, 0.5-0.5/n, n/2+1);
    f.tail(n/2)   = linspace(-0.5, -1.0/n, n/2);
  }
  retourne f;
}

Vecf psd_freqs(entier n, bouléen complexe)
{
  /*
   * si((n & 1) == 0)
  {
    res.tail(n/2) = X.head(n/2);
    res.head(n/2) = X.tail(n/2);

    // Avant shift :
    // k = 0      <=> f = 0
    // k = N/2    <=> f = 1/2
    // k = N/2+1  <=> f = 1/2+1/N <=> -1/2+1/N
    // k = N-1  <=> f = -1/N

    // Après shift :
    // j = 0     <=> k = N/2   <=> f = 1/2 <=> -1/2
    // j = N/2-1 <=> k = N-1   <=> f = -1/N
    // j = N/2   <=> k = 0     <=> f = 0
    // j = N/2+1 <=> k = 1     <=> f = 1
    // j = N-1   <=> k = N/2-1 <=> f = 1/2-1/N
  }
  sinon
  {
    res.tail(1+n/2) = X.head(1+n/2);
    res.head(n/2)   = X.tail(n/2);
     *   X' = \left[X_{N/2+1}, X_{N/2+2}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2}}\right]
  }
   */
  si(complexe)
  {
    soit t0 = -0.5, t1 = 0.5;
    si(est_pair(n))
      t1 -= 1.0/n;
    sinon
      t0 += 1.0/n;

    retourne linspace(t0,t1,n);
  }
  sinon
  {
    soit t0 = 0.0,
         t1 = 0.5;
    si(est_impair(n))
      t1 -= 1.0/n;
    retourne linspace(t0,t1,n/2);
  }
}


static float fréq_estim_candan2(const Veccf &x)
{
  soit n = x.rows();
  soit X = fft(x); // TODO : move fft as a template to avoid this

  soit k = abs(X).index_max();

  soit X2 = X(k),
       X1 = X((k-1+n)%n),
       X3 = X((k+1)%n);

  soit p = real((X1 - X3) / (2.0f*X2-X3-X1));
  p *= std::tan(π/n)/(π/n);
  p = std::atan(p*π/n) / (π/n);

  retourne (k+p) / n;
}

static float fréq_estim_quadratique(const Veccf &x)
{
  soit n = x.rows();
  soit X = abs(fft(x));
  soit k = X.index_max();

  soit AX2 = X(k),
       AX1 = X((k-1+n)%n),
       AX3 = X((k+1)%n);

  soit p = 0.5 * (AX3 - AX1) / (2*AX2-AX3-AX1);
  retourne (k+p) / n;
}

float freqestim(const Veccf &x, FreqEstimMethode m)
{
  float res = 0;

  si(m == FreqEstimMethode::FFT)
  {
    soit n = x.rows();
    soit k = (abs(fft(x))).index_max();

    // k = n/2 <=> n-n/2 = n/2 = -n/2 ?
    // 0 <=> 0
    // 1 <=> n-1
    //
    //si(k >= n/2)
      //k = -(n - k);
    res = ((float) k) / n;
  }
  sinon si(m == FreqEstimMethode::QUADRATIC)
  {
    res = fréq_estim_quadratique(x);
  }
  sinon si(m == FreqEstimMethode::CANDAN2)
  {
    res = fréq_estim_candan2(x);
  }
  si((res > 0.5) && (res <= 1.0))
    res -= 1;
  retourne res;
}

}
