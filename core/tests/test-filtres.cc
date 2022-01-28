#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/filtrage/frat.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/tests.hpp"

using namespace std;
using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::vue;

/** @brief Appel un filtre, non pas directement sur l'ensemble d'un vecteur,
 *  mais en coupant le vecteur en tronçons de longueur BS.
 *  (ceci n'a pas d'intérêt pour un utilisateur, uniquement pour les tests). */
template<typename T>
Vecteur<T> filtre_par_bloc(sptr<FiltreGen<T,T>> filtre, const Vecteur<T> &x, int BS)
{
  vector<Vecteur<T>> lst;
  int N = x.rows(), offset = 0, n = 0;
  while(offset < N)
  {
    auto nl = min(BS, N - offset);
    Vecteur<T> xp = x.segment(offset, nl);
    offset += nl;
    Vecteur<T> yp = filtre->step(xp);
    n += yp.rows();
    lst.push_back(yp);
  }
  Vecteur<T> y(n);

  offset = 0;
  for(auto &v: lst)
  {
    y.segment(offset, v.rows()) = v;
    offset += v.rows();
  }
  return y;
}

static void test_design_rif_prod(int n1, int n2)
{
  msg_majeur("  design_rif_prod: n1={},n2={}", n1, n2);

  ArrayXf h1 = randn(n1);
  ArrayXf h2 = randn(n2);

  ArrayXf hp = design_rif_prod(h1, h2);
  tsd_assert(hp.rows() == (n1+n2-1));

  ArrayXf x = ArrayXf::Zero(n1+n2);
  x(0) = 1;

  ArrayXf y1 = filtrer(h1, x);
  ArrayXf y2 = filtrer(h2, y1);
  ArrayXf yp = filtrer(hp, x);

  //msg("x:{}, y1:{}, y2:{}, yp:{}", x.rows(), y1.rows(), y2.rows(), yp.rows());

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(y2, "", "y2");
    f.subplot().plot(yp, "", "yp");
    f.subplot().plot(y2-yp, "r-o", "Erreur");
    f.afficher();
  }

  auto err = (y2-yp).abs().maxCoeff();
  msg("  erreur = {}", err);
  tsd_assert_msg(err < 1e-5, "Erreur rif prod trop importante.");


}

static void test_design_rif_prod()
{
  msg_majeur("Test design_rif_prod...");

  test_design_rif_prod(11, 11);
  test_design_rif_prod(11, 15);
  test_design_rif_prod(15, 11);

  test_design_rif_prod(10, 20);
  test_design_rif_prod(20, 10);
  test_design_rif_prod(20, 20);

  msg("ok.");
}

// Ne marche pas !
ArrayXf design_rif_rcs_alt(int nc, float bet, float fc)
{
  ArrayXf h_cs    = design_rif_cs (nc, bet, fc);

  ArrayXcf HC = tsd::fourier::fft(h_cs);

  HC = HC.sqrt().eval();
  tsd::fourier::force_csym(HC);
  ArrayXf h2 = tsd::fourier::ifft(HC).real();



  return h2;//tsd::fourier::fftshift(h2);
}

static void test_rcs1()
{
  ArrayXf h1 = design_rif_rcs1(15, 0.2, 4);
  ArrayXf h2 = design_rif_rcs1(15, 0.2, 3);

  if(tests_debug_actif)
  {
    {
      Figure f;
      f.plot(h1, "", "RCS1, OSF = 4");
      f.afficher();
    }
    {
      Figure f;
      f.plot(h2, "", "RCS2, OSF = 3");
      f.afficher();
    }
  }

}

static void test_rcs(int nc)
{
  msg_majeur("Test filtre RCS² VS CS (ncoefs = {})...", nc);

  ArrayXf h_cs    = design_rif_cs (nc, 0.25, 0.25);
  ArrayXf h_rcs   = design_rif_rcs(nc, 0.25, 0.25);
  ArrayXf h_p     = design_rif_prod(h_rcs, h_rcs);


  {
    auto f = affiche_filtre(h_cs);
    f.afficher("CS");
    f = affiche_filtre(h_rcs);
    f.afficher("RCS");
    f = affiche_filtre(h_p);
    f.afficher("PROD");
  }

  ArrayXf h_p2 = h_p.segment(nc/2,nc);

  ArrayXf t1 = linspace(0, nc-1, nc);
  ArrayXf t2 = linspace(-nc/2, nc-1+nc/2, 2*nc-1);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(t1, h_cs, "|bo", "CS");
    f.gcf().plot(t2, h_p, "|go", "RCS*2");
    f.subplot().plot(t1, h_p2 - h_cs, "|or", "Erreur");
    //f.subplot().plot()
    f.afficher(fmt::format("RCS : nc = {}", nc));
  }

  float err = (h_p2 - h_cs).abs().maxCoeff();

  msg("Erreur test RCS : nc = {}, err = {}", nc, err);

  auto err_max = 6e-2;
  if(nc >= 15)
    err_max = 1e-2;
  if(nc >= 63)
    err_max = 7e-4;

  tsd_assert_msg(err <= err_max, "Erreur RCS trop importante.");

  msg("ok");
}

static void test_rcs()
{
  for(auto nc: {5, 15, 31, 63, 127})
    test_rcs(nc);
}


static void test_gaussien_unit(int n, float BT, int osf)
{
  // TODO !!!
  // filtre de comparaison = concat NRZ + Gaussien
  ArrayXf h1 = design_rif_gaussien(n, 0.2);//BT, osf);
  ArrayXf h2 = design_rif_gaussien_telecom(n, BT, osf);

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(h1, "-bo", "RIF Gaussien");
    f.plot(h2, "-go", "RIF Gaussien - telecom");
    f.afficher(fmt::format("RIF Gaussien - BT = {}", BT));
  }

  msg("h1 = {}", h1.transpose());
  msg("h2 = {}", h2.transpose());

  tsd_assert(!h1.hasNaN());
  tsd_assert(!h2.hasNaN());
}


static void test_gaussien()
{
  msg_majeur("Test filtre Gaussien...");
  test_gaussien_unit(15, 0.8, 4);

  test_gaussien_unit(15, 10, 4);

  test_gaussien_unit(15, 100, 4);
  msg("ok.");
}


void test_decimateur()
{
  msg_majeur("Test décimateur...");
  int R = 3;
  auto dec = decimateur<float>(R);

  ArrayXf x = linspace(0, 89, 90);
  ArrayXf y = filtre_par_bloc(dec, x, 4);

  ArrayXf xr = linspace(0, 87, 30);

  tsd_assert(y.rows() == 30);

  msg("xr =\n{}", xr.transpose());
  msg("y =\n{}", y.transpose());

  tsd_assert((y - xr).abs2().sum() == 0);

  msg("ok.");
}


static ArrayXf signal_test()
{
  int n = 5000;
  ArrayXf t = linspace(0,n-1,n);
  return 0.1 * randn(n) + (t*(2*π/n)*20).sin() * (-((t- n/2)/(n/8)).abs()).exp();
}

void test_ligne_a_retard(int δ)
{
  msg_majeur("Test ligne à retard, δ={}...", δ);
  auto lar = ligne_a_retard<float>(δ);

  ArrayXf x = signal_test();
  int n = x.rows();
  ArrayXf y = filtre_par_bloc(lar, x, 311);

  tsd_assert_msg(y.rows() == n, "ligne à retard : pb dim");

  auto err = (y.tail(n-δ) - x.head(n-δ)).maxCoeff();

  msg("Err = {}", err);
  tsd_assert_msg(err == 0, "Ligne à retard : erreur trop importante.");

  msg("ok.");
}

void test_ligne_a_retard()
{
  test_ligne_a_retard(0);
  test_ligne_a_retard(70);
}

void test_filtre_mg(int R)
{
  msg_majeur("Test filtre mg, R={}...", R);
  int n = 1000;
  auto f = filtre_mg<float,double>(R);
  ArrayXf x = randn(n);

  ArrayXf y = filtre_par_bloc(f, x, 80);
  //ArrayXf y = f->step(x);

  tsd_assert_msg(y.rows() == n, "filtre mg : pb dim");

  ArrayXf yref(n);
  for(auto i = 0; i < n; i++)
  {
    int imin  = max(0,i-(R-1));
    int nelem = (i - imin) + 1;
    yref(i) = x.segment(imin, nelem).sum() / R;
  }

  // msg("yref : {}", yref.transpose());
  // msg("y    : {}", y.transpose());
  // msg("err  : {}", (y - yref).transpose());

  float err = (y - yref).abs().maxCoeff();
  tsd_assert_msg(err < 5e-7, "Echec filtre MG : err = {}.", err);
}

void test_filtre_mg()
{
  test_filtre_mg(4);
  test_filtre_mg(11);
}





int test_filtfilt()
{
  msg("Test filtfilt...");
  int n = 500;
  ArrayXf x(n);
  x.head(n/2) = linspace(0,n/2-1,n/2);
  x.tail(n/2) = x.head(n/2).reverse();

  x += randn(n);

  ArrayXf h = design_rif_fen(63, "lp", 0.05);

  ArrayXf y1 = filtrer(h, x);

  ArrayXf y = filtfilt(h, x);

  msg_avert("TODO : automatiser test filtfilt.");

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x, "b-", "x");
    auto c = f.plot(y1, "r-", "filt(x)");
    c.def_epaisseur(3);
    c = f.plot(y, "g-", "filtfilt(x)");
    c.def_epaisseur(3);
    f.afficher("filtfilt");
  }

  return 0;
}


void test_design_biquad()
{
  {
    auto h = design_biquad({.type = BiquadSpec::PASSE_BAS, .f = 0.25, .Q = 0.5});
    analyse_filtre(h).afficher("biquad lp");
  }
  {
    auto h = design_biquad({.type = BiquadSpec::PASSE_HAUT, .f = 0.25, .Q = 0.5});
    analyse_filtre(h).afficher("biquad hp");
  }
  {
    auto h = design_biquad({.type = BiquadSpec::PASSE_BANDE, .f = 0.25, .Q = 0.5});
    analyse_filtre(h).afficher("biquad bp");
  }
  {
    auto h = design_biquad({.type = BiquadSpec::COUPE_BANDE, .f = 0.25, .Q = 0.5});
    analyse_filtre(h).afficher("biquad notch");
  }
}


template<typename T>
void test_design_passe_bas(const T &h, float atten_min = -1, float err_max = 0.1)
{
  //analyse_filtre(h).afficher("Vérification design passe-bas");

  int npts = 2048;
  auto [fr,xm] = frmag(h, npts);

  //int n = h.rows();

  float emax_bp = (xm.head(800) - 1.0f).abs().maxCoeff();
  float emax_bc = (xm.tail(800) - 0.0f).abs().maxCoeff();

  ArrayXf proto = ArrayXf::Zero(npts);
  proto.head(800).setOnes();

  ArrayXf err = proto - xm;

  if(tests_debug_actif)
  {
    Figures fs;
    auto f = fs.subplot();
    f.plot(20*log10(xm), "b-", "réponse");
    f.plot(20*log10(proto), "g-", "prototype");
    fs.subplot().plot(err, "r-");
    fs.afficher("Erreur filtre");
  }

  float atten = -20 * log10(emax_bc);

  msg("Err. max BP = {}, err. max BC = {} ({:.1f} dB).", emax_bp, emax_bc, atten);

  if(max(emax_bp, emax_bc) > err_max)
    echec("Gabarit non respecté.");

  if((atten_min > 0) && (atten < atten_min))
    echec("Atténuation insuffisante ({} dB minimum attendu).", atten_min);
}



void test_filtrage_ola()
{
  msg_majeur("Test filtrage OLA (domaine fréquentiel)...");
  tsd::fourier::FiltreFFTConfig config;

  config.avec_fenetrage         = true;
  config.dim_blocs_temporel     = 512;
  config.nb_zeros_min           = 0;

  config.traitement_freq = [&](ArrayXcf &X)
    {
      int N = 512;
      X.segment(N/32, (30*N)/32).setZero();
    };


  auto [ola, N] = tsd::fourier::filtre_fft(config);


  ArrayXcf x = signal_test();
  ArrayXcf y = filtre_par_bloc<cfloat>(ola, x, 1000);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x,"b-","x");
    f.subplot().plot(y,"g-","filtrage");
    f.afficher("Filtrage OLA freq");
  }
}

static void test_retard(int d)
{
  ArrayXf h = ArrayXf::Ones(d);
  ArrayXf x = ArrayXf::Zero(20);
  x(3) = 1;

  ArrayXf y = filtrer(h, x);
  ArrayXf y2 = filtrer(h, y);

  int index;
  y2.maxCoeff(&index);
  tsd_assert(index == 3+d-1);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x,"b-o","x");
    f.subplot().plot(y,"b-o","y");
    f.subplot().plot(y2,"b-o","y2");
    f.afficher(fmt::format("Test retard - filtre mg({}) - retard {}", d, (d-1)/2.0f));
  }
}

static void test_retard()
{
  msg("Test retard filtre RIF...");
  test_retard(3);
  test_retard(4);
}


void test_filtre_rif()
{
  int nc = 31;
  int n  = nc + 50;
  ArrayXf h = linspace(1, nc, nc);
  ArrayXf x = ArrayXf::Zero(n);
  x(0) = 1;
  ArrayXf y = filtrer(h, x);

  if(y.rows() != n)
    echec("x.rows() = {}, y.rows() = {}", n, y.rows());

  ArrayXf verr = y-vconcat(h, ArrayXf::Zero(n - nc));

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot(311).plot(h, "b|o", "Coefficients");
    f.subplot(312).plot(y, "g|o", "Réponse impulsionnelle");
    f.subplot(313).plot(verr, "r|o", "Erreur");
    f.afficher();
  }



  auto err = verr.abs().maxCoeff();
  //auto err = (y.head(nc)-h).abs().maxCoeff();
  //err = max(err, y.tail(n - nc).abs().maxCoeff());

  msg("err = {}", err);
  if(err > 1e-7)
    echec("Erreur trop importante.");
}


void test_rif_vs_rif_fft()
{
  msg_majeur("Tests rif vs rif fft...");

  ArrayXf h = design_rif_fen(127, "lp", 0.02);

  auto f1 = filtre_rif<float>(h);
  auto f2 = filtre_rif_fft<float>(h);

  ArrayXf x = signal_test();

  ArrayXf y1 = filtre_par_bloc(f1, x, 1000);
  ArrayXf y2 = filtre_par_bloc(f2, x, 1000);


  msg("y1 : {} - {}", y1.minCoeff(), y1.maxCoeff());
  msg("y2 : {} - {}", y2.minCoeff(), y2.maxCoeff());

  auto [y1p,y2p,i,score] = tsd::fourier::aligne_entier(y1, y2);

  ArrayXf err = (y1p - y2p).head(y1p.rows() - 4 * h.rows());
  auto errmax = err.abs().maxCoeff();

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "b-", "x");
    f.subplot().plot(y1, "g-", "rif std");
    f.subplot().plot(y2, "b-", "rif fft");
    auto s = f.subplot(514);
    s.plot(y1p, "g-", "rif std");
    s.plot(y2p, "b-", "rif fft");
    s.titre("std vs fft (alignés)");
    f.subplot().plot(err, "r-", "Erreur (max = {})", errmax);
    f.afficher("Filtre RIF", {1400,800});
  }


  msg("Erreur max = {}", errmax);
  if(errmax > 1e-6)
    echec("Test RIF FFT : erreur trop importante ({})", errmax);
}

void test_filtre_rii()
{
  // Test simple avec un lissage exponentiel (facile de générer un vecteur de référence)
  //  y_n + b_1 y_{n-1} =  a_0 x_n + a_1 x_{n-1} + ...

  msg_majeur("Test filtre RII...");

  float α = 0.1;

  ArrayXf b(2), a(1);
  a(0) = α;
  b(0) = 1.0f;
  b(1) = -(1-α);

  auto h = FRat<float>::rii(a, b);

  msg("H = {}", h);
  msg("H(z^-1) = {}", h.eval_inv_z());

  auto f = filtre_rii<float,float>(h);

  int n = 20;
  ArrayXf x = ArrayXf::Ones(n);

  ArrayXf yref(n), y;

  y = f->step(x);

  yref(0) = α;
  for(auto i = 1; i < n; i++)
    yref(i) = yref(i-1) + α * (1 - yref(i-1));


  if(tests_debug_actif)
  {
    Figures fig;
    auto fi = fig.subplot();
    fi.plot(x, "b-o", "x");
    fi.plot(y, "m-o", "y");
    fi.plot(yref, "g-s", "yref");
    fi = fig.subplot();
    fi.plot(y - yref, "r-o", "erreur");
    fig.afficher("Filtrage RII");
  }

  auto err = (yref - y).abs().maxCoeff();

  msg("Erreur max filtre RII : {}", err);

  if(err > 1e-6f)
    echec("Erreur trop importante");


  msg("ok.");
}




void test_rif_freq()
{
  msg_majeur("Test RIF freq...");
  // Nombre de coefficients souhaités
  int n = 19;
  // Nombre de points du gabarit
  int m = (n+1)/2; // m = 10

  ArrayXf d = ArrayXf::Zero(m);
  d.head(m/2) = 1.0f;

  ArrayXf h = design_rif_freq(n, d);
  tsd_assert(h.rows() == n);
  test_design_passe_bas(h, -1, 0.2);

  auto [fr,xm] = frmag<float>(h);
  ArrayXf fr1 = design_rif_freq_freqs(n);

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(fr,xm,"b-", "Réponse obtenue");

    f.plot(fr1, d, "gs", "Points d'échantillonnage (gabarit)");
    f.titre("Réponse fréquentielle");
    f.afficher();
  }


  // Vérification réponse fréquentielle = réponse attendue
  ArrayXcf y = repfreq(h, fr1);

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(fr1, d.abs(), "bo-", "Réponse attendue");
    f.plot(fr1, y.abs(), "gs-", "Réponse réelle");
    f.plot(fr1, (d.abs()-y.abs()).abs(), "r-", "Erreur");
    f.afficher();
  }


  auto err = (y.abs() - d.abs()).abs().maxCoeff();

  msg("Erreur sur les points d'interpolation : {}", err);

  // Vérification filtre de type I (phase linéaire)
  auto err_lp = (h.head(n/2) - h.tail(n/2).reverse()).abs().maxCoeff();

  msg("Erreur phase linéaire : {}", err_lp);

  tsd_assert_msg(err + err_lp < 1e-6, "Design RIF freq invalide.");

  msg("ok.");
}


static void test_riia()
{
  msg_majeur("Test design RIIA...");
  for(auto structure : {"ellip", "butt", "cheb1", "cheb2"})
  {
    msg_majeur("Test design RIIA {}...", structure);

    // Ondulations max : 0.1 dB dans la bande passante, ou -60 dB dans la bande coupée
    //auto h = design_riia(4, "lp", structure, 0.2, 3, 30);
    auto h = design_riia(12, "lp", structure, 0.25, 0.1, 60);
    test_design_passe_bas(h);
  }
  msg("ok");
}

int test_filtres()
{


  test_retard();
  test_design_rif_prod();




  test_rif_freq();

  test_decimateur();

  test_ligne_a_retard();

  test_filtre_rii();

  test_riia();

  {

    ArrayXf h;

    msg_majeur("Test des design RIF...");

    msg("Test RCS...");
    h = design_rif_rcs(127, 0.2, 0.25);
    test_design_passe_bas(h);


    msg("Test CS...");
    h = design_rif_cs(127, 0.2, 0.25);
    test_design_passe_bas(h);


    test_rcs1();
    test_rcs();

    msg("Test RIF FEN...");
    h = design_rif_fen(127, "lp", 0.25);
    test_design_passe_bas(h);

    msg("Test RIF FEN - kaiser...");
    h = design_rif_fen_kaiser("lp", 0.25, 60, 0.01);
    test_design_passe_bas(h);

    {
      msg("test RIF freq");
      int n = 127;
      int m = 512;
      ArrayXf d(m);

      d.head(m/2) = 1.0f;
      d.tail(m/2) = 0.0f;

      //stdo.def_dossier_sortie("./test-log/rif-freq");
      h = design_rif_freq(n, d);
      //stdo.flush();
      msg("nb coefs : {}, demandés : {}", h.rows(), n);
      tsd_assert(h.rows() == n);
      test_design_passe_bas(h);
    }



    msg("Test RIF EQ...");
    // TEST TROP LONG
    for(int nc : {15, 63, 127, 255})//, 1023})
    {
      msg("  NC = {}", nc);

      // Ne fonctionne pas au-delà
      if(nc <= 127)
      {
        msg("Mode Cheby....");
        h = tsd::filtrage::design_rif_eq(nc, {{0, 0.2, 1}, {0.3, 0.5, 0}}, false);
        test_design_passe_bas(h);
      }



      msg("Mode Linf....");
      h = tsd::filtrage::design_rif_eq(nc, {{0, 0.2, 1}, {0.3, 0.5, 0}}, tests_debug_actif, false);
      test_design_passe_bas(h);
    }

    //
    //ArrayXf design_rif_eq2(int nc, IArrayXf D, IArrayXf W, bool debug)

    //msg("Test RIF EQ2...");
    /*for(int nc : {15, 63, 127})
    {
      msg("  NC = {}", nc);
      h = tsd::filtrage::design_rif_eq(nc, {{0, 0.2, 1}, {0.3, 0.5, 0}}, true);
      test_design_passe_bas(h);
    }*/


  }


  test_gaussien();


  ArrayXf h = tsd::filtrage::design_rif_rcs(63, 0.1, 1.0/(2*8));

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(h, "", "Réponse impulsionnelle");
    f.afficher();
  }

  test_filtre_rif();


  if(test_filtfilt())
    return -1;



  test_rif_vs_rif_fft();


  test_filtrage_ola();


  test_filtre_mg();

  test_design_biquad();
  return 0;
}
