#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"

namespace tsd::fourier {


std::tuple<ArrayXf, ArrayXf> psd_welch(const ArrayXcf &x, int N, const std::string fen)
{
  ArrayXf freqs = psd_freqs(N);
  ArrayXf S = ArrayXf::Zero(N);
  ArrayXf f = tsd::filtrage::fenetre(fen, N, false);

  for(auto i = 0; i + N < x.rows(); i += N/2)
  {
    ArrayXcf xp = x.segment(i, N) * f;
    S += fftshift(fft(xp).abs2());
  }

  return {freqs, 10.0f * S.log10()};
}

ArrayXf psd_freqs(int n, bool complexe)
{
  /*
   * if((n & 1) == 0)
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
  else
  {
    res.tail(1+n/2) = X.head(1+n/2);
    res.head(n/2)   = X.tail(n/2);
     *   X' = \left[X_{N/2+1}, X_{N/2+2}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2}}\right]
  }
   */
  if(complexe)
  {
    if((n & 1) == 0)
      return linspace(-0.5,0.5-1.0/n,n);
    else
      return linspace(-0.5+1.0/n, 0.5, n);
  }
  else
  {
    if((n & 1) == 0)
      return linspace(0,0.5,n/2);
    else
      return linspace(0,0.5-1.0/n,n/2);
  }
}


static float fréq_estim_candan2(IArrayXcf x)
{
  auto n = x.rows();
  //ArrayXcf xp = x;
  ArrayXcf X = fft(x); // TODO : move fft as a template to avoid this

  int k;
  X.abs().maxCoeff(&k);

  auto X2 = X(k);
  auto X1 = X((k-1+n)%n);
  auto X3 = X((k+1)%n);

  auto p = ((X1 - X3) / (2.0f*X2-X3-X1)).real();
  p *= std::tan(π/n)/(π/n);
  p = std::atan(p*π/n) / (π/n);

  return (k+p) / n;
}

static float fréq_estim_quadratique(IArrayXcf x)
{
  auto n = x.rows();
  ArrayXcf X = fft(x);

  int k;
  X.abs().maxCoeff(&k);

  auto AX2 = abs(X(k));
  auto AX1 = abs(X((k-1+n)%n));
  auto AX3 = abs(X((k+1)%n));

  float p = 0.5 * (AX3 - AX1) / (2*AX2-AX3-AX1);
  return (k+p) / n;
}

float freqestim(IArrayXcf x, FreqEstimMethode m)
{
  float res = 0;

  if(m == FreqEstimMethode::FFT)
  {
    ArrayXcf X = fft(x);
    int k;
    int n = X.rows();
    X.abs().maxCoeff(&k);

    // k = n/2 <=> n-n/2 = n/2 = -n/2 ?
    // 0 <=> 0
    // 1 <=> n-1
    //
    //if(k >= n/2)
      //k = -(n - k);
    res = ((float) k) / n;
  }
  else if(m == FreqEstimMethode::QUADRATIC)
  {
    res = fréq_estim_quadratique(x);
  }
  else if(m == FreqEstimMethode::CANDAN2)
  {
    res = fréq_estim_candan2(x);
  }
  if((res > 0.5) && (res <= 1.0))
  {
    res -= 1;
  }
  return res;
}

}
