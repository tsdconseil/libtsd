#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"


#define DEBUG_DELAY_EST 0




namespace tsd::fourier {


// ???
// Estimation (à l'échantillon près) de délais
class EstimateurDelais
{
public:
  EstimateurDelais();
  int configure(int N);
  // Score : entre 0 et 1 (1 <=> signaux identiques, au délais et à la phase prés)
  int step(IArrayXcf x, IArrayXcf y, int &delais, float &score);
private:
  ArrayXf fenetre;
  //ArrayXcf x2, y2;
  int N;
  //tsd::fourier::FFTCorrelateur correlateur;
};

EstimateurDelais::EstimateurDelais()
{
  //configure(1);
}

int EstimateurDelais::configure(int N)
{
  this->N = N;
  //correlateur.configure(N);
  //x2.resize(N, 1);
  //y2.resize(N, 1);
  this->fenetre = tsd::filtrage::fenetre("hann", N);
  return 0;
}

int estimation_delais_entier(IArrayXcf x, IArrayXcf y, float &score)
{
  auto [xp, yp] = pad_zeros(x, y);
  EstimateurDelais est;
  est.configure(xp.rows());
  int delais;
  est.step(xp, yp, delais, score);
  return delais;
}


//std::tuple<ArrayXcf, ArrayXcf, int, float> aligne_entier(IArrayXcf x, IArrayXcf y)

template<typename T>
  std::tuple<Vecteur<T>, Vecteur<T>, int, float> aligne_entier(const Vecteur<T> &x1, const Vecteur<T> &y1)
{
  float score;

  ArrayXcf x = x1, y = y1;

  // Zéro-padding pour avoir la même dimension et qu'elle soit une puissance de 2
  auto [xp, yp] = pad_zeros(x, y, true);


  //infos("aligne entier : #x=%d, #y=%d, #xp=%d, #yp=%d.", (int) x.rows(), (int) y.rows(), (int) xp.rows(), (int) yp.rows());

  int d = estimation_delais_entier(xp, yp, score);

  ArrayXcf xa, ya;
  if(d == 0)
  {
    xa = x1;
    ya = y1;
    //xa = xp;
    //ya = yp;
  }
  else if(d < 0)
  {
    xa = x1.tail(x1.rows() + d);
    ya = y1;

    //xa = x1.tail(x1.rows()+d);
    //ya = y1;//.head(x1.rows());
    //xa = xp.tail(xp.rows()+d);
    //ya = yp.head(xa.rows());
  }
  else
  {
    xa = x1;
    ya = y1.tail(y1.rows() - d);
    //xa = y1.tail(y1.rows()-d);
    //ya = x1.head(ya.rows());
    //ya = yp.tail(yp.rows()-d);
    //xa = xp.head(ya.rows());
  }

  if(xa.rows() > ya.rows())
    xa = xa.head(ya.rows()).eval();
  else if(xa.rows() < ya.rows())
    ya = ya.head(xa.rows()).eval();

  //infos("aligne entier : d = %d, xa.rows() = %d, ya.rows() = %d", d, (int) xa.rows(), (int) ya.rows());

  tsd_assert(xa.rows() == ya.rows());

  if constexpr(est_complexe<T>())
    return {xa, ya, d, score};
  else
    return {xa.real(), ya.real(), d, score};
}


// Score : entre 0 et 1 (1 <=> signaux identiques, au d�lais et � la phase pr�s)
int EstimateurDelais::step(IArrayXcf x, IArrayXcf y, int &delais, float &score)
{
  tsd_assert(x.rows() == N);
  tsd_assert(y.rows() == N);

  ArrayXcf x2 = x * fenetre;
  ArrayXcf y2 = y * fenetre;

  auto e1 = std::sqrt(x2.abs2().mean());
  auto e2 = std::sqrt(y2.abs2().mean());


  auto [lags, corr] = ccorr(x2, y2);

  //ArrayXcf corr;
  //correlateur.step(corr, x2, y2);
  int index;
  score = corr.abs().maxCoeff(&index);
  //score = std::sqrt((float) N) * score / (e1 * e2 + 1e-50);
  score /= (e1 * e2 + 1e-50);

  // sortie du corrélateur = corrélation non normalisée


  delais = index;
  if(delais > N / 2)
    delais = -(N - delais);

  //infos("ed: index = %d (délais %d), score = %f", index, delais, score);

  // Si energie nulle, met le score à zéro
  //score *= (e1 * e2) / (1.0 + e1 * e2);
  //score *= (e1 * e2) / (1.0-20 + e1 * e2);

  // index=0:   lag = 0
  // index=1:   lag = 1
  // ...
  // index=n-1: lag = n-1 == -1 from the periodicity hypothesis.


  return 0;
}




#if 0
static int verifie_delais(int n, int d)
{
  EstimateurDelais ed;
  ed.configure(n);
  int delais;
  float score;
  ArrayXcf x(n), y(n);
  ArrayXf w = tsd::filtrage::fenetre("hann", n);

  //w = w.cwiseProduct(w);

  /*for(auto i = 0; i < n; i++)
  {
    x(i) = std::sin((2*pi*i*60)/n) * w(i);

    if((i+d<n) && (i+d>=0))
      y(i) = std::sin((2*pi*(i+d)*60)/n) * w(i+d);
    else
      y(i) = 0;
  }*/

  x.setZero(n);
  y.setZero(n);

  x(n/2) = 1.0;
  y((n/2)+d) = 1.0;


  ed.step(x, y, delais, score);
  //infos("n=%d, delais=%d->%d, score=%.2f", n, d, delais, score);
  if(delais != d)
    return -1;
  return 0;
}
#endif



class TAD
{
public:

  int num_test = 0;
  int reception_ok = 0;

  void gere_buffer(const ArrayXXcf &x)
  {
    tsd_assert(x.cols() == 2);
    tsd_assert(x.rows() == 1024);


    for(auto i = 0; i < 10; i++)
    {
      auto a = x(i,0), b = x(i,1);
      msg("{}: x = {}:{}, y ={}:{}", i, a.real(), a.imag(), b.real(), b.imag());
    }

    for(auto i = 0; i < 1024; i++)
    {
      cfloat fi(i, 0);
      if(num_test == 0)
      {
        tsd_assert(x(i,0) == fi);
        tsd_assert(x(i,1) == -fi);
      }
      else if(num_test == 1)
      {
        if(i < 5)
        {
          tsd_assert(x(i,0) == fi);
          tsd_assert(x(i,1) == -(fi+5.0f));
        }
        else
        {
          tsd_assert(x(i,0) == fi+5.0f);
          tsd_assert(x(i,1) == -(fi+5.0f));
        }
      }
      else if(num_test == 2)
      {
        tsd_assert(x(i,0) == fi + 5.0f);
        tsd_assert(x(i,1) == -fi);
      }
      else if(num_test == 3)
      {
        if(i < 750)
        {
          tsd_assert(x(i,0) == fi);
          tsd_assert(x(i,1) == -fi);
        }
        else
        {
          tsd_assert(x(i,0) == fi+5.0f);
          tsd_assert(x(i,1) == -fi);
        }
      }
    }

    reception_ok = 1;
  }

  int test_alignement_delais()
  {
    AlignementSignal as;

    // TODO
    //as.sortie.add_listener(this, &TAD::gere_buffer);
    int n = 1500;
    ArrayXXcf x(n,1), y(n,1);
    for(auto i = 0; i < n; i++)
    {
      x(i) = i;
      y(i) = - (float) i;
    }

    msg("TEST 0");
    as.configure(1024); // Taille des paquets de sortie
    num_test = 0;
    as.step(x, y, 0);
    tsd_assert(reception_ok);

    msg("TEST 1");
    as.configure(1024); // Taille des paquets de sortie
    num_test = 1;
    as.step(x.block(0, 0, 5, 1), y.block(0, 0, 5, 1), 5);
    as.step(x.block(5, 0, n-5, 1), y.block(5, 0, n-5, 1), 0);
    tsd_assert(reception_ok);

    msg("TEST 2");
    as.configure(1024); // Taille des paquets de sortie
    num_test = 2;
    as.step(x.block(0, 0, n/2, 1), y.block(0, 0, n/2, 1), -5);
    as.step(x.block(n/2, 0, n/2, 1), y.block(n/2, 0, n/2, 1), -5);
    tsd_assert(reception_ok);


    msg("TEST 3");
    as.configure(1024); // Taille des paquets de sortie
    num_test = 3;
    as.step(x.block(0, 0, n/2, 1), y.block(0, 0, n/2, 1), 0);
    as.step(x.block(n/2, 0, n/2, 1), y.block(n/2, 0, n/2, 1), -5);
    tsd_assert(reception_ok);

    return reception_ok;
  }
};


int test_alignement_delais()
{
  TAD tad;
  return tad.test_alignement_delais() ? 0 : -1;
}

template
  std::tuple<Vecteur<float>, Vecteur<float>, int, float>
    aligne_entier<float>(const Vecteur<float> &x1, const Vecteur<float> &y1);

template
  std::tuple<Vecteur<cfloat>, Vecteur<cfloat>, int, float>
    aligne_entier<cfloat>(const Vecteur<cfloat> &x1, const Vecteur<cfloat> &y1);




#if 0
int DelayEstimation::setup(uint32_t n)
{
  this->n = n;
  rfft_correlator.setup(n);
  fft_correlator.configure(n);
  corr.resize(n);
  x2.resize(n);
  return 0;
}

int DelayEstimation::step(const ArrayXf &x0, const ArrayXf &x1, float &res)
{
# if 0
  int ipart;
  uint32_t i;

# if DEBUG_DELAY_EST
  SciFile scifile("./build/out.bin");
  scifile.add_matrix<Tsample,1>(x0, "x0");
  scifile.add_matrix<Tsample,1>(x1, "x1");
# endif

  // (1) Cross-correlation to retrieve the integer part of delay
  rfft_correlator.step(corr, x0, x1);

# if DEBUG_DELAY_EST
  scifile.add_matrix<Tsample,1>(corr, "cr");
# endif

  // TODO
  // argabsmax(corr, ipart);

  // Negative lags
  if(ipart >= (int) n/2)
    ipart = -(n - ipart);

  //dsp_trace("%s: ipart = %d.", __func__, ipart);

# if DEBUG_DELAY_EST
  // Check with FFT correlator that the result is the same
  {
    int ipartc;
    Matrix<std::complex<Tsample>,1> x0c(n), x1c(n), corrc(n);
    for(i = 0; i < n; i++)
    {
      x0c(i) = x0(i);
      x1c(i) = x1(i);
    }
    fft_correlator.step(corrc, x0c, x1c);

    scifile.add_matrix<Tsample,1>(corrc, "crc");

    cargabsmax(corrc, ipartc);

    dsp_trace("%s: ipartc = %d, n = %d.", __func__, ipartc, n);

    // Negative lags
    if(ipartc >= (int) n/2)
      ipartc = -(n - ipartc);

    dsp_trace("%s: ipartc = %d.", __func__, ipartc);

    if(ipartc != ipart)
    {
      dsp_warning("%s: ipartc (%d) != ipart (%d).", __func__, ipartc, ipart);
      //return -1;
      ipart = ipartc;
    }
  }
# endif



  // (2) Shift x1 by ipart ==> x2
  for(i = 0; i < n; i++)
  {
    int i2 = i - ipart;

    if((i2 >= 0) && (i2 < (int) n))
      x2(i) = x1(i2);
    else
      x2(i) = 0;
  }

# if DEBUG_DELAY_EST
  {
    // Check if integer part is now 0
    int ipart2;
    rfft_correlator.step(corr, x0, x2);
    argabsmax(corr, ipart2);

    if(ipart2)
      dsp_anomaly("%s: Failed to zero the ipart of delay (residual: %d).", __func__, ipart2);
  }
# endif


  // (3) Compute the slope
  Vecteur<cfloat> X0(n/2), X1(n/2), X2(n/2), scratch(n/2), twiddles(n/2), itwiddles(n/2);

  transforms::fft_twiddles(twiddles, n/2);
  transforms::ifft_twiddles(itwiddles, n/2);

  // (1) Compute Real FFT on x0 and x1
  transforms::rfft_compute(X0, x0, scratch, twiddles);
  transforms::rfft_compute(X1, x2, scratch, twiddles);

  X0(0).imag(0);
  X1(0).imag(0);

  for(i = 0; i < n/2; i++)
    X2(i) = X0(i) * conj(X1(i));

# if DEBUG_DELAY_EST
  scifile.add_matrix<Tsample,1>(X2, "xc2");
  scifile.add_matrix<Tsample,1>(X0, "xc0");
  scifile.add_matrix<Tsample,1>(X1, "xc1");

  // Check with FFT that the result is the same
  {
    Matrix<std::complex<Tsample>,1> x0c(n), x1c(n), corrc(n);
    for(i = 0; i < n; i++)
    {
      x0c(i) = x0(i);
      x1c(i) = x2(i);
    }

    Matrix<std::complex<Tsample>,1> X0(n), X1(n), X2(n/2),
                                         scratch(n), ftwiddles(n);
    transforms::fft_twiddles(ftwiddles, n);

    // (1) Compute Real FFT on x0 and x1
    transforms::fft_compute(X0, x0c, scratch, ftwiddles);
    transforms::fft_compute(X1, x1c, scratch, ftwiddles);

    for(i = 0; i < n/2; i++)
      X2(i) = X0(i) * conj(X1(i));

    int ismax;
    cargabsmax(X2, ismax);
    float arg = std::arg(X2(ismax));
    float del = - (arg * n) / (2 * pi * ismax);

    dsp_trace("%s: with fft: ismax = %d, arg = %f, del = %f samples",
          __func__, ismax, arg, del);
  }
# endif

  // (2) Compute D
  float D = 0, S = 0;
  for(i = 1; i < n/2; i++)
  {
    float mag = norm(X2(i));
    float arg = std::arg(X2(i));
    D += i * arg * mag;
    S += mag * i * i;
  }

  //int ismax;
  //cargabsmax(X2, ismax);

  //float arg = std::arg(X2(ismax));

  // XD(i) = sum x(k+1) e^{-2pi j i k / n)
  //       = sum x(k) e^{-2pi j i (k-1) / n)
  //       = sum x(k) e^{-2pi j i k / n) e^{2pi j i / n)
  //       = X(i) * e^{2pi j i / n)
  // => X_tau(i) = X(i) * e^{2pi j tau * i / n)

  // => Arg(X_tau(i) * X(i)*) = 2 pi tau * i / n
  // => tau = arg * n / (2 * pi * i)

  //float del = - (arg * n) / (2 * pi * ismax);

  //dsp_trace("%s: ismax = %d, arg = %f, del = %f samples",
    //    __func__, ismax, arg, del);


  D = - (D / S) * n / (2 * π);

  res = -(D + ipart);

  //dsp_trace("%s: result = %f samples (%f + %d).", __func__, res, D, ipart);


# endif
  return 0;
}
#endif





}
