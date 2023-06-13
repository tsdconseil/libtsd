#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"

//using namespace std;

static void test_fftplan(entier n)
{
  soit x    = Veccf::random(n);
  soit plan = tfrplan_création();
  soit y1   = plan->step(x),
       y2   = plan->step(x, non);

  soit err1 = abs(y1 - fft(x)).valeur_max(),
       err2 = abs(y2 - ifft(x)).valeur_max();

  msg("Test fft plan[n={}] : erreur fft = {}, erreur ifft = {}", n, err1, err2);
  assertion((err1 < 1e-6) && (err2 < 1e-6));
}


static void test_fftplan()
{
  for(auto i : {8, 16, 18, 19, 101})
    test_fftplan(i);
}

static void test_rfftplan()
{
  soit x    = randu(101);
  soit plan = rtfrplan_création();
  soit y    = plan->step(x);
  soit err  = abs(y - rfft(x)).valeur_max();

  msg("Test rfftplan: erreur = {}", err);

  assertion(err < 1e-6);
}

static void test_fftshift(entier n)
{
  soit x = linspace(0,n-1,n),
       y = fftshift(x);

  assertion(y.rows() == n);

  Vecf yref(n);

  soit m = n/2; // n = 2 * m + 1

  si(n & 1)
  {
    // n = 2 * m + 1
    yref.head(m)   = linspace(m+1, n-1, m);
    yref.tail(m+1) = linspace(0, m, m+1);
  }
  sinon
  {
    // n = 2 * m
    yref.head(m) = linspace(m, n-1, m);
    yref.tail(m) = linspace(0, m-1, m);
  }

  assertion(y.est_approx(yref));
}


static void test_fftshift()
{
  msg("test fftshift...");
  test_fftshift(15);
  test_fftshift(16);
}

static void test_goertzel()
{
  {
    soit n = 100;
    soit x = randn(n);

    msg_majeur("Test Goertzel (sur un buffer)...");

    pour(auto i = 0; i < n / 2; i++)
    {
      soit freq = (i * 1.0f) / n;
      soit s1   = goertzel(x, freq);
      soit X    = fft(x);
      soit s2   = 2 * std::norm(X(i)) / square(x).somme();
      soit err  = abs(s2 - s1);
      si(err > 1.5e-5)
        échec("Echec test Goertzel, i = {}, s1 = {}, s2 = {}, err = {:e}.", i, s1, s2, err);
    }
  }

  {
    msg_majeur("Test Goertzel (au fil de l'eau)...");
    soit N = 8,
         n = N * 8;
    soit freq = 5.0f/N;
    soit f = filtre_goertzel(freq, N);
    soit x = randn(n),
         y = f->step(x);
    assertion_msg(y.rows() == x.rows() / N, "Filtre Goertzel : pb décimation.");

    //msg("y = {}", y.transpose());

    soit errmax = 0.0f;
    pour(auto i = 0; i + N <= n; i += N)
    {
      soit s0  = y(i / N),
           s1  = goertzel(x.segment(i, N), freq),
           err = abs(s1 - s0);

      //msg("theo(i = {} ... {}) = {}, calc = {}, err = {:e}", i, i+N-1, s1, s0, err);
      errmax = max(errmax, err);
      si(err > 5e-7)
        échec("Echec test Goertzel fil de l'eau, i = {}, s0 = {}, s1 = {}, err = {:e}", i, s0, s1, err);
    }
    msg("errmax = {}", errmax);
  }
}

static void test_reechan()
{
  soit n = 16;

  msg_majeur("Test resample...");

  soit x  = linspace(0, 1 - 1.0f/n, n),
       x1 = tsd::fourier::rééchan_freq(x, 2),
       x2 = tsd::rééchan(x, 2);

  assertion(x1.rows() == 2 * n);
  assertion(!x1.hasNaN());

  soit t  = x.clone(),
       t1 = linspace(0, 1 - 1.0f/(2*n), x1.rows()),
       t2 = linspace(0, 1, x2.rows());

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(t,  x, "bo", "x");
    f.plot(t1, x1, "g-", "resample fft");
    f.plot(t2, x2, "m-", "resample rif");
    f.afficher("Resample");
  }

  assertion_msg(x1.rows() == 2 * n, "Resample de facteur 2 : nombre d'échantillons invalide");
  soit x1b = sousech(x1, 2);
  assertion_msg(x1b.rows() == n, "sousech : nombre d'échantillons invalide");

  soit err = abs(x1b - x).valeur_max();

  msg("Erreur max resample : {}", err);

  assertion_msg(err < 1e-5, "Resample : erreur trop importante ({})", err);

  msg("ok.");
}


// TODO: à déplacer dans fourier.hpp
template<typename T>
Vecteur<T> tfd(const Vecteur<T> &x, bouléen inv = non)
{
  soit n = x.rows();
  Vecteur<T> Y(n);
  pour(auto i = 0; i < n; i++)
  {
    cdouble s = 0;
    double signe = inv ? 1 : -1;
    pour(auto k = 0; k < n; k++)
      s += ((cdouble) x(k)) * polar(1.0, signe*(2*π*k*i)/n);
    s /= sqrt((float) n);
    Y(i) = s;
  }
  retourne Y;
}


template<typename Tin, typename Tout>
static void test_fft_valide(entier n, bouléen alea, bouléen inv)
{

  msg("Vérification validité FFT (n={}, alea={}, inv={})...", n, alea, inv);

  Vecteur<Tin> x;

  si(alea)
    x = Vecteur<Tin>::random(n);
  sinon
    x = Vecteur<Tin>::ones(n);


  Vecteur<Tout> X;

  si constexpr(est_complexe<Tin>())
  {
    msg("  - calcul FFT...");
    si(inv)
      X = ifft(x);
    sinon
      X = fft(x);
  }
  sinon
  {
    msg("  - calcul RFFT...");
    X = rfft(x);
  }

  assertion(X.rows() == n);


  // TODO : utiliser TFD
  Vecteur<Tout> Y(n);
  msg("  - calcul DFT naive (pour comparaison)...");
  pour(auto i = 0; i < n; i++)
  {
    cfloat s = 0;
    float signe = inv ? 1 : -1;
    pour(auto k = 0; k < n; k++)
      s += x(k) * polar(1.0f, signe*(2*π_f*k*i)/n);
    s /= sqrt((float) n);
    Y(i) = s;
  }

  msg("Dim(X) = {}, dim(Y) = {}", X.rows(), Y.rows());
  msg("X = {}", X);
  msg("Y = {}", Y);
  msg("X-Y = {}", X-Y);

  soit T = Veccf::valeurs({1.0f});
  msg("T = {}", T);
  msg("sq(T) = {}", square(T));

  soit V = Vecf::valeurs({1.0f});
  msg("V = {}", V);
  msg("abs(V) = {}", abs(V));

  msg("abs(T) = {}", abs(T));

  // Plante
  msg("abs(X-Y) = {}", abs(X-Y));
  soit err = abs(X - Y).valeur_max();

  //fmt::print("FFT = \n{}\n, DFT = \n{}\n", X, Y);

  msg("  - Erreur FFT (n = {}) = {}", n, err);

  si((err >= 1e-2) && (n < 100))
  {
    msg("TFR = {}", X);
    msg("TFD = {}", Y);
  }

  assertion_msg(err < 1e-2, "Erreur FFT");
}

static void test_fft()
{
  msg_majeur("Tests FFT...");

  pour(auto i : {16, 1, 2, 3, 4, 5, 8, 10, 17, 128, 129, 1024})
  {
    test_fft_valide<cfloat,cfloat>(i, non, non);
    test_fft_valide<cfloat,cfloat>(i, oui, non);
    test_fft_valide<cfloat,cfloat>(i, oui, oui);

    // Test RFFT
    test_fft_valide<float,cfloat>(i, non, non);
    test_fft_valide<float,cfloat>(i, oui, non);
  }


  soit n = 1024;
  soit x = cos(linspace(0,8*2*π,n));

  si(tests_debug_actif)
  {
    soit X = ifft(x);
    Figure f;
    f.plot(abs(X));
    f.afficher("ifft");
  }


  soit X   = fft(x);
  soit x2  = real(ifft(X));
  soit err = sqrt(square(x2 - x).moyenne());

  soit xc   = x.clone();
  soit Xc   = fft(xc);
  soit x2c  = real(ifft(Xc));
  soit errc = sqrt(abs2(x2c - xc).moyenne());


  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "", "x");
    f.subplot().plot(abs(X), "", "X");
    soit s = f.subplot();
    s.plot(x2, "b-", "x2");
    s.plot(imag(ifft(X)), "g-", "ifft");
    f.afficher("fftifft");
  }

  msg("erreur fft + ifft = {}", errc);
  msg("erreur rfft + ifft = {}", err);

  si(err > 5e-6)
    échec("Test fft : erreur trop importante.");
}

static Vecf test_signal(entier n = 15*1024)
{
  Vecf x0(n);
  soit fen = fenetre("hn", n / 2);

  /** Ratio between signal freq and sampling freq, between 0 and 0.5 */
  soit fe             = 60e6f,
       fs             = 100e3f,
       fs_normalisée  = fs / fe,
       periode        = 1 / fs_normalisée;

  pour(auto i = 0; i < n; i++)
  {
    soit w = 0.0f;
    si((i >= n/4) && (i < n - n / 4))
      w = fen(i - n / 4);
    x0(i) = w * sin((2.0*π*i)/periode);
  }
  retourne x0;
}


// Test de delais
template<typename T>
static void test_delais_fractionnaire(float d)
{
  soit x0 = test_signal().as<cfloat>();

  si constexpr (std::is_same<T, cfloat>::value)
  {
    x0 *= polar(1.0f, -π_f/4);
  }

  soit x1 = tsd::fourier::délais(x0, d);
  soit n  = x0.rows();
  assertion(x0.rows() == x1.rows());

  Vecf err;

  si(d >= 0)
    err = abs(x1.tail(n-d) - x0.head(n-d));
  sinon
    err = abs(x1.head(n+d) - x0.tail(n+d));

  soit em = err.valeur_max();

  string tps = "float";
  si(tests_debug_actif)
  {
    Figures f;

    si constexpr (std::is_same<T, cfloat>::value)
    {
      soit s = f.subplot(221);
      s.plot(real(x0), "b-", "x0");
      s.plot(real(x1), "g-", "x1");
      s.plot(imag(x1), "r-", "x1 (imag)");
      s.titre(sformat("x0, x1 = x0(n-{})", d));
      tps = "cfloat";

      s = f.subplot(222);
      s.plot(imag(x1), "r-", "x1 (imag)");
      //f.plot(x0.segment(n/2-500, 1000).real(), "b-", "x0");
      //f.plot(x1.segment(n/2-500+250, 1000).real(), "g-", "x1");
    }
    sinon
    {
      soit s = f.subplot(221);
      s.plot(x0, "b-", "x0");
      s.plot(x1, "g-", "x1");
      s.titre(sformat("x0, x1 = x0(n-{})", d));
    }
    //f.titre(sformat("Délais fractionnaire - {} = {}", tps, d));
    f.subplot(223).plot(err, "r-", sformat("erreur (max = {})", em));

    soit s = f.subplot(224);
    s.plot(linspace(n/2,n/2+19,20), real(x0.segment(n/2, 20)), "b-", "x0");
    s.plot(linspace(n/2,n/2+19,20), real(x1.segment(n/2, 20)), "g-", "x0");
    s.titre("Zoom t = [n/2,n/2+19]");

    f.afficher(sformat("test-delais-frac-{}-{}.png", tps, d));
  }

  msg("Délais fractionnaire - {} = {}, erreur max = {}.", tps, d, em);

  si(floor(d) == d)
    assertion(em < 1e-2);
  sinon
    assertion(em < 1e-1);
}



// Test délais + estimation_délais
static void test_delais_unitaire(float délais_vrai, float snr_db, entier type_signal, entier N)
{
  msg("test avec délais = {}, type signal = {}, N = {}...",
      délais_vrai, (type_signal == 0) ? "impulsion gaussienne" : "randn", N);

  Veccf x0;

  si(type_signal == 0)
    x0 = test_signal(N).as<cfloat>();
  sinon
    x0 = randn(N);


  si(floor(délais_vrai) == délais_vrai)
    msg(" -> délais entier");
  sinon
    msg(" -> délais fractionnaire");

  // Apply some delay
  soit x1 = délais(x0, délais_vrai);

  //float es = x0.square().mean();
  // snr_db = 10 log_10(E[S]/E[N])
  // E[S]/E[N] = 10^(snr_db / 10)
  //float sn = pow(10.0, snr_db / 10.0);
  // E[N] = sigma²
  //float σ = sqrt((1.0 / sn) * es);
  //msg("SNR = {},  dB, σ = {}", snr_db, σ);

  // pour vérifier que le score est indépendant de la norme des signaux
  x0 *= 7;
  x1 *= 4;

  Veccf X0 = x0, X1 = x1;
  soit [d, score] = estimation_délais(X0, X1);

  si(tests_debug_actif)
  {
    msg("Figures...");
    Figure f;
    f.plot(x0/7, "b-", "x0");
    f.plot(x1/4, "g-", "x1");
    f.titre(sformat("Délais réf={:.2g}, détecté={}, score={}", délais_vrai, d, score));
    f.afficher(sformat("test-delais-{:.2f}.png", délais_vrai));

    {
      Figure f;
      soit [lags, c] = xcorrb(x0, x1);
      f.plot(lags, c);
      f.afficher("XCORRB");
    }
  }

  soit err_pos    = abs(d - délais_vrai),
       err_score  = abs(score - 1),
       tol_pos = (N == 32) ? 0.1f : 0.02f,
       tol_score = 0.4f;

  msg("Délais vrai = {}, estimé = {} (erreur = \033[33m{}\033[0m), score = \033[33m{}\033[0m", délais_vrai, d, err_pos, score);



  assertion_msg(err_pos < tol_pos/*0.02*/, "test_delais_unitaire : erreur délais trop importante : {}.", err_pos);
  assertion_msg(err_score < tol_score/*1e-4*/, "test_delais_unitaire : le score devrait être proche de 1 (il vaut {}).", score);
}



static tuple<Vecf, Veccf> xcorr_ref(const Veccf &x, const Veccf &y, entier m = -1, bouléen biais = non)
{
  soit n = x.rows();
  si(m == -1)
    m = n;
  assertion(x.rows() == y.rows());

  soit lags = linspace(-(m-1), m-1, 2*m-1);
  Veccf c(2*m-1);

  // Délais nul et négatifs (-(n-1) ... 0)
  pour(auto i = 0; i < m; i++)
  {
    cfloat s = 0;
    soit ns = 0;
    pour(auto k = 0; k + i < n; k++)
    {
      s += x(k + i) * conj(y(k));
      ns++;
    }
    si(!biais)
      s /= ns;
    sinon
      s /= n;
    c(m-1-i) = s;
  }

  // Délais positifs (1..n-1) -> index = n, ..., 2*n-2

  // n = m = 3, i = m-1:    k+2+1 < 3
  pour(auto i = 1; i < m; i++)
  {
    cfloat s = 0;
    soit ns = 0;
    pour(auto k = 0; k + i < n; k++)
    {
      s += x(k) * conj(y(k + i));
      ns++;
    }
    si(!biais)
      s /= ns;
    sinon
      s /= n;
    c(i+m-1) = s;
  }

  retourne {lags, c};
}




static void test_xcorr(bouléen biaisé)
{
  msg_majeur("Tests xcorr{}...", biaisé ? "b" : "");
  pour(entier n : {1, 2, 3, 10, 15, 16, 21, 32})
  {
    msg("Test xcorr{}..., n = {}", biaisé ? "b" : "", n);
    soit a1 = Veccf::random(n),
         a2 = Veccf::random(n);

    Vecf lags;
    Veccf c;

    si(biaisé)
      tie(lags, c) = xcorrb(a1, a2);
    sinon
      tie(lags, c) = xcorr(a1, a2);

    // n=1 : 0
    // n=2 : -1 0 1
    // n=3 : -2 -1 0 1 2
    // ...

    soit [lags_ref, cref] = xcorr_ref(a1, a2, n, biaisé);

    //msg("n = {}, lags =\n {}", n, lags.transpose());

    assertion_msg(lags.rows() == 2*n-1, "test xcorr : dim lags invalide.");

    msg("lags     = {}", lags);
    msg("lags ref = {}", lags_ref);
    assertion_msg((lags_ref - lags).valeur_max() == 0, "test xcorr : lags invalides.");

    soit err = abs(cref - c).valeur_max();
    msg("Erreur = {:g}", err);
    si(err > 1e-5)
    {
      msg_avert("Echec test xcorr: n={}, err =\n{}", n, cref - c);
      msg("lags = {}", lags);
      msg("xcorr_ref = {}", cref);
      msg("xcorr     = {}", c);
    }
    assertion_msg(err < 1e-5, "Erreur xcorr");
    msg("ok.");
  }
}

static void test_ccorr()
{
  msg_majeur("Tests ccorr...");
  pour(entier n : {1, 2, 3, 10, 15, 16, 21, 32})
  {
    msg("Test ccorr..., n = {}", n);
    soit a1 = Veccf::random(n),
         a2 = Veccf::random(n);
    soit [lags, c] = ccorr(a1, a2);
    Veccf cref(n);

    //msg("a1 conj = {}", a1.conjugate());

    pour(auto i = 0; i < n; i++)
    {
      cfloat s = 0;
      pour(auto k = 0; k < n; k++)
        s += a1(k) * conj(a2((k + i) % n));
      cref(i) = s / ((float) n);
    }

    soit err = abs(cref - c).valeur_max();
    msg("Erreur = {:g}", err);

    msg("a1 = {}, a2 = {}", a1, a2);
    msg("cref = {}, c = {}", cref, c);

    msg("Erreur ccorr = {}", err);
    assertion_msg(err < 1e-5, "Erreur ccorr");

    msg("ok.");
  }
}

static void test_align_entier()
{
  msg_majeur("Tests de la fonction 'aligne_entier'...");

  pour(entier dref = -50; dref <= 50; dref += 25)
  {
    msg("Test align entier, dref = {}...", dref);
    soit x = test_signal();
    soit n = x.rows();
    soit y = Vecf::zeros(n);

    si(dref >= 0)
      y.tail(n-dref) = x.head(n-dref);
    sinon
      y.head(n+dref) = x.tail(n+dref);

    soit [x2, y2, d, s] = aligne_entier(x, y);
    assertion_msg(x2.rows() == y2.rows(), "");
    soit err = abs(x2 - y2).valeur_max();

    msg("  aligne_entier : d={}, s={}, err={}", d, s, err);

    assertion_msg(d == dref, "aligne_entier : délais invalide.");
    assertion_msg(abs(s - 1) < 1e-3, "aligne_entier : score invalide.");
    assertion_msg(err < 1e-3, "aligne_entier : vecteur invalide.");
  }
}


// TODO: à utiliser ou supprimer
static float csym_erreur(const Veccf &X)
{
  soit n = X.rows();

  Vecf verr;
  si((n % 2) == 0)
    verr = abs(X.segment(1, n/2 - 1) - X.tail(n/2 - 1).reverse().conjugate());
  sinon
    verr = abs(X.segment(1, n/2) - X.tail(n/2).reverse().conjugate());

  soit [err, idx] = verr.max();
  idx++;

  err = max(err, abs(imag(X(0))));
  si((n % 2) == 0)
    err = max(err, abs(imag(X(n/2))));

  retourne err;
}

static void test_csym(entier n)
{
  soit X = Veccf::random(n);

  csym_forçage(X);

  soit x1 = ifft(X),             // moins précis
       x2 = tfd<cfloat>(X, oui); // plus précis

  soit err1 = abs(imag(x1)).valeur_max(),
       err2 = abs(imag(x2)).valeur_max();

  msg("test csym(n={:04d}) : err max partie imaginaire : avec FFT: {}, avec TFD: {}", n, err1, err2);

  assertion_msg(err1 < 1e-3, "Erreur csym trop importante");
  assertion_msg(err2 < 1e-7, "Erreur csym trop importante");
}

void test_fourier()
{


  test_fftshift();

  {
    soit n     = 16*1024*2;
    soit Omega = - 250.0f * 2 * π_f / n;

    soit f1 = polar(1.0f, 54 * Omega),
         f2 = polar(1.0f, (n-54) * Omega);

    soit err = abs(f1 - conj(f2));

    fmt::print("f1 = {}, f2 = {}\n", f1, f2);
    fmt::print("Erreur = {}\n", err);
  }

  test_fft();

  pour(auto n: {3, 4, 5, 63, 64, 511, 512, 1000, 1001})
    test_csym(n);

  test_fftplan();
  test_rfftplan();
  test_goertzel();

  test_reechan();


  msg_majeur("Test délais...");
  pour(float f: {0.f, 250.f, 1.f, 10.f, -10.f, 0.5f, 0.1f, 1.5f})
  {
    test_delais_fractionnaire<cfloat>(f);
    test_delais_fractionnaire<float>(f);
  }

  test_ccorr();
  test_xcorr(non);
  test_xcorr(oui);


  pour(auto N: {32, 1024, 15*1024})
    pour(auto type_signal : {0})//), 1})
      pour(float f: {0.f, 1.f, 10.f, 20.f, 30.f, 40.f, -50.f, 11.f, 1.1f})
      {
        si((N == 32) && (type_signal == 0))
          continue;
        si((N == 32) && (abs(f) > 10))
          continue;
        test_delais_unitaire(f, 40, type_signal, N);
      }

  test_align_entier();


  msg("test ola complexité...");

  float C;
  entier Nf, Nz, Ne;
  soit M = 2560;
  ola_complexité_optimise(M, C, Nf, Nz, Ne);
  msg("M={} -> Nf={}, Nz={}, Ne={}", M, Nf, Nz, Ne);
  assertion(Ne < 100000);
}


