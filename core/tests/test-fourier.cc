#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/tests.hpp"


using namespace std;
using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;



void test_fftshift(int n)
{
  ArrayXf x = linspace(0,n-1,n);
  ArrayXf y = fftshift(x);

  tsd_assert(y.rows() == n);

  ArrayXf yref(n);

  int m = n/2; // n = 2 * m + 1

  if(n & 1)
  {
    // n = 2 * m + 1
    yref.head(m) = linspace(m+1, n-1, m);
    yref.tail(m+1) = linspace(0, m, m+1);
  }
  else
  {
    // n = 2 * m
    yref.head(m) = linspace(m, n-1, m);
    yref.tail(m) = linspace(0, m-1, m);
  }

  tsd_assert(y.isApprox(yref));
}


void test_fftshift()
{
  msg("test fftshift...");
  test_fftshift(15);
  test_fftshift(16);
}

void test_goertzel()
{
  {
    int n = 100;
    ArrayXf x = randn(n);

    msg_majeur("Test Goertzel (sur un buffer)...");

    for(auto i = 0; i < n / 2; i++)
    {
      float freq = (i * 1.0f) / n;
      float s1 = goertzel(x, freq);
      ArrayXcf X = fft(x);
      float s2 = 2 * norm(X(i)) / x.square().sum();
      float err = abs(s2 - s1);
      if(err > 1.5e-5)
        echec("Echec test Goertzel, i = {}, s1 = {}, s2 = {}, err = {:e}.", i, s1, s2, err);
    }
  }

  {
    msg_majeur("Test Goertzel (au fil de l'eau)...");
    int N = 8;
    int n = N * 8;
    ArrayXf x = randn(n);
    float freq = 5.0f/N;
    auto f = filtre_goertzel(freq, N);
    ArrayXf y = f->step(x);
    tsd_assert_msg(y.rows() == x.rows() / N, "Filtre Goertzel : pb décimation.");

    //msg("y = {}", y.transpose());

    float errmax = 0;
    for(auto i = 0; i + N <= n; i += N)
    {
      float s0  = y(i / N);
      float s1  = goertzel(x.segment(i, N), freq);
      float err = abs(s1 - s0);

      //msg("theo(i = {} ... {}) = {}, calc = {}, err = {:e}", i, i+N-1, s1, s0, err);
      errmax = max(errmax, err);
      if(err > 5e-7)
        echec("Echec test Goertzel fil de l'eau, i = {}, s0 = {}, s1 = {}, err = {:e}", i, s0, s1, err);
    }
    msg("errmax = {}", errmax);
  }
}

static void test_resample()
{
  int n = 16;

  msg_majeur("Test resample...");
  ArrayXf x = linspace(0, 1 - 1.0f/n, n);
  msg("  fourier...");
  ArrayXf x1 = tsd::fourier::reechan_freq(x, 2);

  tsd_assert(x1.rows() == 2 * n);

  msg("  rif...");
  ArrayXf x2 = tsd::reechan(x, 2);


  //msg("x  : [{} - {}]", x.minCoeff(), x.maxCoeff());
  //msg("x1 : [{} - {}]", x1.minCoeff(), x1.maxCoeff());
  //msg("x2 : [{} - {}]", x2.minCoeff(), x2.maxCoeff());

  tsd_assert(!x1.hasNaN());

  ArrayXf t = x;
  ArrayXf t1 = linspace(0, 1 - 1.0f/(2*n), x1.rows());
  ArrayXf t2 = linspace(0, 1, x2.rows());

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(t,  x, "bo", "x");
    f.plot(t1, x1, "g-", "resample fft");
    f.plot(t2, x2, "m-", "resample rif");
    f.afficher("Resample");
  }

  tsd_assert_msg(x1.rows() == 2 * n, "Resample de facteur 2 : nombre d'échantillons invalide");
  ArrayXf x1b = sousech(x1, 2);
  tsd_assert_msg(x1b.rows() == n, "sousech : nombre d'échantillons invalide");

  auto err = (x1b - x).abs().maxCoeff();

  msg("Erreur max resample : {}", err);

  tsd_assert_msg(err < 1e-5, "Resample : erreur trop importante ({})", err);

  msg("ok.");
}


// TODO: à déplacer dans fourier.hpp
template<typename T>
Vecteur<T> tfd(const Vecteur<T> &x, bool inv = false)
{
  int n = x.rows();
  Vecteur<T> Y(n);
  for(auto i = 0; i < n; i++)
  {
    cdouble s = 0;
    double signe = inv ? 1 : -1;
    for(auto k = 0; k < n; k++)
      s += ((cdouble) x(k)) * polar(1.0, signe*(2*π*k*i)/n);
    s /= sqrt((float) n);
    Y(i) = s;
  }
  return Y;
}


template<typename Tin, typename Tout>
static int test_fft_valide(int n, bool alea, bool inv)
{

  msg("Vérification validité FFT (n={}, alea={}, inv={})...", n, alea, inv);

  Vecteur<Tin> x;

  if(alea)
    x = Vecteur<Tin>::Random(n);
  else
    x = Vecteur<Tin>::Ones(n);


  Vecteur<Tout> X;

  if constexpr(est_complexe<Tin>())
  {
    msg("  - calcul FFT...");
    if(inv)
      X = ifft(x);
    else
      X = fft(x);
  }
  else
  {
    msg("  - calcul RFFT...");
    X = rfft(x);
  }

  tsd_assert(X.rows() == n);


  // TODO : utiliser TFD
  Vecteur<Tout> Y(n);
  msg("  - calcul DFT naive (pour comparaison)...");
  for(auto i = 0; i < n; i++)
  {
    cfloat s = 0;
    float signe = inv ? 1 : -1;
    for(auto k = 0; k < n; k++)
      s += x(k) * polar(1.0f, signe*(2*π_f*k*i)/n);
    s /= sqrt((float) n);
    Y(i) = s;
  }


  float err = (X - Y).abs().maxCoeff();

  //fmt::print("FFT = \n{}\n, DFT = \n{}\n", X, Y);

  msg("  - Erreur FFT (n = {}) = {}", n, err);

  if((err >= 1e-2) && (n < 100))
  {
    msg("TFR = {}", X);
    msg("TFD = {}", Y);
  }

  tsd_assert_msg(err < 1e-2, "Erreur FFT");

  return 0;
}

static int test_fft()
{
  msg_majeur("Tests FFT...");

  for(auto i : {16, 1, 2, 3, 4, 5, 8, 10, 17, 128, 129, 1024})
  {
    if(test_fft_valide<cfloat,cfloat>(i, false, false))
      return -1;
    if(test_fft_valide<cfloat,cfloat>(i, true, false))
      return -1;
    if(test_fft_valide<cfloat,cfloat>(i, true, true))
      return -1;

    // Test RFFT
    if(test_fft_valide<float,cfloat>(i, false, false))
      return -1;
    if(test_fft_valide<float,cfloat>(i, true, false))
      return -1;

  }


  int n = 1024;
  ArrayXf x = cos(linspace(0,8*2*π,n));

  if(tests_debug_actif)
  {
    ArrayXcf X = ifft(x);
    Figure f;
    f.plot(X.abs());
    f.afficher("ifft");
  }


  ArrayXcf X = fft(x);
  ArrayXf x2 = ifft(X).real();
  float err = sqrt((x2 - x).square().mean());

  ArrayXcf xc = x;
  ArrayXcf Xc = fft(xc);
  ArrayXf x2c = ifft(Xc).real();
  float errc = sqrt((x2c - xc).abs2().mean());


  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "", "x");
    f.subplot().plot(X.abs(), "", "X");
    auto s = f.subplot();
    s.plot(x2, "b-", "x2");
    s.plot(ifft(X).imag(), "g-", "ifft");
    f.afficher("fftifft");
  }

  msg("erreur fft + ifft = {}", errc);
  msg("erreur rfft + ifft = {}", err);

  if(err > 5e-6)
    return -1;

  return 0;
}

static ArrayXf test_signal(int n = 15*1024)
{

  ArrayXf x0(n);
  ArrayXf wnd = fenetre("hn", n / 2);

  /** Ratio between signal freq and sampling freq, between 0 and 0.5 */
  float fe = 60E6;
  float fs = 100e3;
  float normalized_fs = fs / fe;
  float periode = 1 / normalized_fs;
  //float es = 0.0;

  for(auto i = 0; i < n; i++)
  {
    float w = 0;
    if((i >= n/4) && (i < n - n / 4))
      w = (float) wnd(i - n / 4);
    x0(i) = w * sin((2.0*π*i)/periode);
    //x0(i) *= 10000.0;
  }
  return x0;
}


// Test de delais
template<typename T>
static int test_delais_fractionnaire(float d)
{
  Vecteur<T> x0 = test_signal();


  if constexpr (std::is_same<T, cfloat>::value)
  {
    x0 *= std::polar(1.0f, -π_f/4);
  }


  Vecteur<T> x1 = tsd::fourier::delais(x0, d);
  int n = x0.rows();
  tsd_assert(x0.rows() == x1.rows());

  ArrayXf err;

  if(d >= 0)
    err = (x1.tail(n-d) - x0.head(n-d)).abs();
  else
    err = (x1.head(n+d) - x0.tail(n+d)).abs();
  auto em = err.maxCoeff();

  string tps = "float";
  if(tests_debug_actif)
  {

    Figures f;

    if constexpr (is_same<T, cfloat>::value)
    {
      auto s = f.subplot(221);
      s.plot(x0.real(), "b-", "x0");
      s.plot(x1.real(), "g-", "x1");
      s.plot(x1.imag(), "r-", "x1 (imag)");
      s.titre(fmt::format("x0, x1 = x0(n-{})", d));
      tps = "cfloat";

      s = f.subplot(222);
      s.plot(x1.imag(), "r-", "x1 (imag)");
      //f.plot(x0.segment(n/2-500, 1000).real(), "b-", "x0");
      //f.plot(x1.segment(n/2-500+250, 1000).real(), "g-", "x1");
    }
    else
    {
      auto s = f.subplot(221);
      s.plot(x0, "b-", "x0");
      s.plot(x1, "g-", "x1");
      s.titre(fmt::format("x0, x1 = x0(n-{})", d));
    }
    //f.titre(fmt::format("Délais fractionnaire - {} = {}", tps, d));
    f.subplot(223).plot(err, "r-", fmt::format("erreur (max = {})", em));

    auto s = f.subplot(224);
    s.plot(linspace(n/2,n/2+19,20), x0.segment(n/2, 20).real(), "b-", "x0");
    s.plot(linspace(n/2,n/2+19,20), x1.segment(n/2, 20).real(), "g-", "x0");
    s.titre("Zoom t = [n/2,n/2+19]");

    f.afficher(fmt::format("test-delais-frac-{}-{}.png", tps, d));
  }

  msg("Délais fractionnaire - {} = {}, erreur max = {}.", tps, d, em);

  if(floor(d) == d)
  {
    tsd_assert(em < 1e-2);
  }
  else
  {
    tsd_assert(em < 1e-1);
  }


  return 0;
}



// Test de delais + estimation_delais_entier
static int test_delais_unitaire(float refdelay, float snr_db = 40.0)
{
  //float delay;

  msg("test with delay = {}...", refdelay);

  //tsd::transforms::DelayEstimation<T, float, float> estimator;
  //assert(estimator.setup(n) == 0);

  ArrayXf x0 = test_signal();
  //int n = x0.rows();

  float es = x0.square().mean();

  // Apply some delay




  ArrayXf x1 = tsd::fourier::delais(x0, refdelay);
  //tsd::transforms::intdelay<T>(x1, x0, round(refdelay));

  // snr_db = 10 log_10(E[S]/E[N])
  // E[S]/E[N] = 10^(snr_db / 10)
  float sn = pow(10.0, snr_db / 10.0);
  // E[N] = sigma�

  float σ = sqrt((1.0 / sn) * es);

  msg("SNR = {},  dB, σ = {}", snr_db, σ);

  // Apply some noise
  //x0 += sigma * tsd::randn(n);
  //x1 += sigma * tsd::randn(n);


  // Pour vérifier que le score est indépendant de la norme des signaux
  x1 *= 4;
  x0 *= 7;


  // USE FFT to apply the delay
  float score;
  ArrayXcf X0 = x0, X1 = x1;
  int d = tsd::fourier::estimation_delais_entier(X0, X1, score);
  //assert(estimator.process(x0, x1, delay) == 0);

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x0, "b-", "x0");
    f.plot(x1, "g-", "x1");
    f.titre(fmt::format("Délais réf={:.2g}, détecté={}, score={}", refdelay, d, score));
    f.afficher(fmt::format("test-delais-{:.2f}.png", refdelay));
  }

  //if(refdelay != 0)
    //assert(cutils::TestUtil::verifie_valeur(delay, refdelay, 0.144, "delay") == 0);

  auto err = d - refdelay;

  msg("Délais vrai = {}, estimé = {}, erreur = {}, score = {}", refdelay, d, err, score);

  if(abs(err) >= 1)
  {
    msg_erreur("Erreur trop importante.");
    return -1;
  }

  if(abs(score - 1) >= 1e-4)
  {
    msg_erreur("test_delais_unitaire : le score devrait être proche de 1 (il vaut {}).", score);
    return -1;
  }


  return 0;
}



std::tuple<ArrayXf, ArrayXcf> xcorr_ref(const ArrayXcf &x, const ArrayXcf &y, int m = -1, bool biais = false)
{
  int n = x.rows();
  if(m == -1)
    m = n;
  tsd_assert(x.rows() == y.rows());
  ArrayXf lags = linspace(-m-1, m-1, 2*m-1);
  ArrayXcf c(2*m-1);

  // Délais nul et négatifs (-(n-1) ... 0)
  for(auto i = 0; i < m; i++)
  {
    cfloat s = 0;
    int ns = 0;
    for(auto k = 0; k + i < n; k++)
    {
      s += x(k + i) * conj(y(k));
      ns++;
    }
    if(!biais)
      s /= ns;
    else
      s /= n;
    c(m-1-i) = s;
  }

  // Délais positifs (1..n-1) -> index = n, ..., 2*n-2

  // n = m = 3, i = m-1:    k+2+1 < 3
  for(auto i = 1; i < m; i++)
  {
    cfloat s = 0;
    int ns = 0;
    for(auto k = 0; k + i < n; k++)
    {
      s += x(k) * conj(y(k + i));
      ns++;
    }
    if(!biais)
      s /= ns;
    else
      s /= n;
    c(i+m-1) = s;
  }

  return {lags, c};
}




int test_xcorr()
{
  msg_majeur("Tests xcorr...");
  for(int n : {1, 2, 3, 10, 15, 16, 21, 32})
  {
    msg("Test xcorr..., n = {}", n);
    ArrayXcf a1 = ArrayXf::Random(n);
    ArrayXcf a2 = ArrayXf::Random(n);
    auto [lags, c] = xcorr(a1, a2);

    // n=1 : 0
    // n=2 : -1 0 1
    // n=3 : -2 -1 0 1 2
    // ...

    auto [lags_ref, cref] = xcorr_ref(a1, a2, n, false);

    //msg("n = {}, lags =\n {}", n, lags.transpose());

    tsd_assert_msg(lags.rows() == 2*n-1, "test xcorr : dim lags invalide.");
    tsd_assert_msg((lags_ref - lags).maxCoeff() == 0, "test xcorr : lags invalides.");

    auto err = (cref - c).abs().maxCoeff();
    msg("Erreur = {:g}", err);
    if(err > 1e-5)
    {
      msg_avert("Echec test xcorr: n={}, err =\n{}", n, (cref - c).transpose());
      msg("lags = {}", lags.transpose());
      msg("xcorr_ref = {}", cref.transpose());
      msg("xcorr     = {}", c.transpose());
    }
    tsd_assert_msg(err < 1e-5, "Erreur xcorr");
    msg("ok.");
  }
  return 0;
}

int test_ccorr()
{
  msg_majeur("Tests ccorr...");
  for(int n : {1, 2, 3, 10, 15, 16, 21, 32})
  {
    msg("Test ccorr..., n = {}", n);
    ArrayXcf a1 = ArrayXf::Random(n);
    ArrayXcf a2 = ArrayXf::Random(n);
    auto [lags, c] = ccorr(a1, a2);
    ArrayXcf cref(n);

    for(auto i = 0; i < n; i++)
    {
      cfloat s = 0;
      for(auto k = 0; k < n; k++)
      {
        s += a1(k) * conj(a2((k + i) % n));
      }
      cref(i) = s / ((float) n);
      //fmt::print("corr[{}] = {}, xcorr = {}\n", i, s, c(i));
    }

    //msg("cref =\n{}", cref);
    //msg("ccor =\n{}", c);

    auto err = (cref - c).abs().maxCoeff();
    msg("Erreur = {:g}", err);
    tsd_assert_msg(err < 1e-5, "Erreur ccorr");

    msg("ok.");
  }
  return 0;
}

int test_align_entier()
{
  msg_majeur("Tests de la fonction 'aligne_entier'...");

  for(int dref = -50; dref <= 50; dref += 25)
  {
    msg("Test align entier, dref = {}...", dref);
    ArrayXcf x = test_signal();
    int n = x.rows();
    ArrayXcf y = ArrayXcf::Zero(n);

    if(dref >= 0)
      y.tail(n-dref) = x.head(n-dref);
    else
      y.head(n+dref) = x.tail(n+dref);

    auto [x2, y2, d, s] = aligne_entier(x, y);
    tsd_assert_msg(x2.rows() == y2.rows(), "");
    float err = (x2 - y2).abs().maxCoeff();

    msg("  aligne_entier : d={}, s={}, err={}", d, s, err);

    tsd_assert_msg(d == dref, "aligne_entier : délais invalide.");
    tsd_assert_msg(abs(s - 1) < 1e-3, "aligne_entier : score invalide.");
    tsd_assert_msg(err < 1e-3, "aligne_entier : vecteur invalide.");
  }

  return 0;
}

int test_csym(int n)
{
  ArrayXcf X = ArrayXcf::Random(n);

  //if(n <= 5)
  //  msg("X avant csym = \n{}", X);

  force_csym(X);

  //if(n <= 5)
  //  msg("X après csym = \n{}", X);

  ArrayXcf x1 = ifft(X);              // moins précis
  ArrayXcf x2 = tfd<cfloat>(X, true); // plus précis

  auto err1 = x1.imag().abs().maxCoeff();
  auto err2 = x2.imag().abs().maxCoeff();

  msg("test csym(n={:04d}) : err max partie imaginaire : avec FFT: {}, avec TFD: {}", n, err1, err2);

  tsd_assert_msg(err1 < 1e-3, "Erreur csym trop importante");
  tsd_assert_msg(err2 < 1e-7, "Erreur csym trop importante");
  return 0;
}

int test_fourier()
{
  test_goertzel();

  test_fftshift();

  {
    int n = 16*1024*2;
    float Omega = - 250.0f * 2 * π / n;

    cfloat f1 = polar(1.0f, 54 * Omega);
    cfloat f2 = polar(1.0f, (n-54) * Omega);

    float err = abs(f1 - conj(f2));

    fmt::print("f1 = {}, f2 = {}\n", f1, f2);
    fmt::print("Erreur = {}\n", err);
    //exit(0);
  }



  if(test_fft())
    return -1;

  for(auto n: {3, 4, 5, 63, 64, 511, 512, 1000, 1001})
    test_csym(n);

  test_resample();


  msg_majeur("Test delais...");
  for(float f: {0.f, 250.f, 1.f, 10.f, -10.f, 0.5f, 0.1f, 1.5f})
  {
    test_delais_fractionnaire<cfloat>(f);
    test_delais_fractionnaire<float>(f);
  }

  if(test_ccorr())
    return -1;

  if(test_xcorr())
    return -1;



  for(float f: {0.f, 1.f, 10.f, 20.f, 30.f, 40.f, -50.f, 11.f, 1.1f})
    if(test_delais_unitaire(f))
      return -1;

  test_align_entier();

  return 0;
}


