#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"

using namespace std;

/** @brief Appel un filtre, non pas directement sur l'ensemble d'un vecteur,
 *  mais en coupant le vecteur en tronçons de longueur BS.
 *  (ceci n'a pas d'intérêt pour un utilisateur, uniquement pour les tests). */
template<typename T>
Vecteur<T> filtre_par_bloc(sptr<FiltreGen<T,T>> filtre, const Vecteur<T> &x, entier BS)
{
  vector<Vecteur<T>> lst;
  soit N = x.rows(), offset = 0, n = 0;
  tantque(offset < N)
  {
    soit nl = min(BS, N - offset);
    Vecteur<T> xp = x.segment(offset, nl);
    offset += nl;
    Vecteur<T> yp = filtre->step(xp);
    n += yp.rows();
    lst.push_back(yp);
  }
  Vecteur<T> y(n);

  offset = 0;
  pour(auto &v: lst)
  {
    y.segment(offset, v.rows()) = v;
    offset += v.rows();
  }
  retourne y;
}

static void test_design_rif_prod(entier n1, entier n2)
{
  msg_majeur("  design_rif_prod: n1={},n2={}", n1, n2);

  soit h1 = randn(n1);
  soit h2 = randn(n2);

  soit hp = design_rif_prod(h1, h2);
  tsd_assert(hp.rows() == (n1+n2-1));

  soit x = Vecf::zeros(n1+n2);
  x(0) = 1;

  soit y1 = filtrer(h1, x);
  soit y2 = filtrer(h2, y1);
  soit yp = filtrer(hp, x);

  //msg("x:{}, y1:{}, y2:{}, yp:{}", x.rows(), y1.rows(), y2.rows(), yp.rows());

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(y2, "", "y2");
    f.subplot().plot(yp, "", "yp");
    f.subplot().plot(y2-yp, "r-o", "Erreur");
    f.afficher();
  }

  soit err = abs(y2-yp).valeur_max();
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
Vecf design_rif_rcs_alt(entier nc, float bet, float fc)
{
  soit h_cs    = design_rif_cs (nc, bet, fc);

  soit HC = fft(h_cs);

  HC = sqrt(HC);
  csym_forçage(HC);
  retourne real(ifft(HC));
}

static void test_rcs1()
{
  soit h1 = design_rif_rcs1(15, 0.2, 4);
  soit h2 = design_rif_rcs1(15, 0.2, 3);

  si(tests_debug_actif)
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

static void test_rcs(entier nc)
{
  msg_majeur("Test filtre RCS² VS CS (ncoefs = {})...", nc);

  soit h_cs    = design_rif_cs (nc, 0.25, 0.25);
  soit h_rcs   = design_rif_rcs(nc, 0.25, 0.25);
  soit h_p     = design_rif_prod(h_rcs, h_rcs);


  {
    soit f = plot_filtre(h_cs);
    f.afficher("CS");
    f = plot_filtre(h_rcs);
    f.afficher("RCS");
    f = plot_filtre(h_p);
    f.afficher("PROD");
  }

  soit h_p2 = h_p.segment(nc/2,nc);

  soit t1 = linspace(0, nc-1, nc);
  soit t2 = linspace(-nc/2, nc-1+nc/2, 2*nc-1);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(t1, h_cs, "|bo", "CS");
    f.gcf().plot(t2, h_p, "|go", "RCS*2");
    f.subplot().plot(t1, h_p2 - h_cs, "|or", "Erreur");
    //f.subplot().plot()
    f.afficher(sformat("RCS : nc = {}", nc));
  }

  soit err = abs(h_p2 - h_cs).valeur_max();

  msg("Erreur test RCS : nc = {}, err = {}", nc, err);

  soit err_max = 6e-2;
  si(nc >= 15)
    err_max = 1e-2;
  si(nc >= 63)
    err_max = 7e-4;

  tsd_assert_msg(err <= err_max, "Erreur RCS trop importante.");

  msg("ok");
}

static void test_rcs()
{
  pour(auto nc: {5, 15, 31, 63, 127})
    test_rcs(nc);
}


static void test_gaussien_unit(entier n, float BT, entier osf)
{
  // TODO !!!
  // filtre de comparaison = concat NRZ + Gaussien
  soit h1 = design_rif_gaussien(n, 0.2);//BT, osf);
  soit h2 = design_rif_gaussien_telecom(n, BT, osf);

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(h1, "-bo", "RIF Gaussien");
    f.plot(h2, "-go", "RIF Gaussien - telecom");
    f.afficher(sformat("RIF Gaussien - BT = {}", BT));
  }

  msg("h1 = {}", h1);
  msg("h2 = {}", h2);

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
  entier R = 3;
  soit dec = decimateur<float>(R);

  soit x = linspace(0, 89, 90);
  soit y = filtre_par_bloc(dec, x, 4);

  soit xr = linspace(0, 87, 30);

  tsd_assert(y.rows() == 30);
  tsd_assert(abs2(y - xr).somme() == 0);

  msg("ok.");
}


static Vecf signal_test()
{
  soit n = 5000;
  soit t = linspace(0,n-1,n);
  retourne 0.1 * randn(n) + sin(t*(2*π/n)*20) * exp(-abs((t- n/2)/(n/8)));
}

void test_ligne_a_retard(entier δ)
{
  msg_majeur("Test ligne à retard, δ={}...", δ);
  soit lar = ligne_a_retard<float>(δ);

  soit x = signal_test();
  entier n = x.rows();
  soit y = filtre_par_bloc(lar, x, 311);

  tsd_assert_msg(y.rows() == n, "ligne à retard : pb dim");

  soit err = (y.tail(n-δ) - x.head(n-δ)).valeur_max();

  msg("Err = {}", err);
  tsd_assert_msg(err == 0, "Ligne à retard : erreur trop importante.");

  msg("ok.");
}

void test_ligne_a_retard()
{
  test_ligne_a_retard(0);
  test_ligne_a_retard(70);
}

void test_filtre_mg(entier R)
{
  msg_majeur("Test filtre mg, R={}...", R);
  entier n = 1000;
  soit f = filtre_mg<float,double>(R);
  soit x = randn(n);
  soit y = filtre_par_bloc(f, x, 80);

  tsd_assert_msg(y.rows() == n, "filtre mg : pb dim");

  Vecf yref(n);
  pour(soit i = 0; i < n; i++)
  {
    soit imin  = max(0,i-(R-1));
    soit nelem = (i - imin) + 1;
    yref(i) = x.segment(imin, nelem).somme() / R;
  }

  soit err = abs(y - yref).valeur_max();

  msg("err = {}", err);

  tsd_assert_msg(err < 5e-7, "Echec filtre MG : err = {}.", err);
}

void test_filtre_mg()
{
  test_filtre_mg(4);
  test_filtre_mg(11);
  test_filtre_mg(20);
}







entier test_filtfilt()
{
  msg("Test filtfilt...");
  entier n = 500;
  Vecf x(n);
  x.head(n/2) = linspace(0,n/2-1,n/2);
  x.tail(n/2) = x.head(n/2).reverse();

  x += randn(n);

  soit h  = design_rif_fen(63, "lp", 0.05);
  soit y1 = filtrer(h, x);
  soit y  = filtfilt(h, x);

  msg_avert("TODO : automatiser test filtfilt.");

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(x, "b-", "x");
    soit c = f.plot(y1, "r-", "filt(x)");
    c.def_epaisseur(3);
    c = f.plot(y, "g-", "filtfilt(x)");
    c.def_epaisseur(3);
    f.afficher("filtfilt");
  }

  retourne 0;
}


void test_design_biquad()
{
  {
    soit h = design_biquad({.type = BiquadSpec::PASSE_BAS, .f = 0.25, .Q = 0.5});
    plot_filtre(h).afficher("biquad lp");
  }
  {
    soit h = design_biquad({.type = BiquadSpec::PASSE_HAUT, .f = 0.25, .Q = 0.5});
    plot_filtre(h).afficher("biquad hp");
  }
  {
    soit h = design_biquad({.type = BiquadSpec::PASSE_BANDE, .f = 0.25, .Q = 0.5});
    plot_filtre(h).afficher("biquad bp");
  }
  {
    soit h = design_biquad({.type = BiquadSpec::COUPE_BANDE, .f = 0.25, .Q = 0.5});
    plot_filtre(h).afficher("biquad notch");
  }
}


struct TestDesignConfig
{
  TypeFiltre type     = TypeFiltre::PASSE_BAS;
  float atten_min     = -1;
  float err_max       = 0.1;
  float fc_3dB_théo   = -1;
  float fc_6dB_théo   = -1;
};


template<typename T>
void test_design(const T &h, const TestDesignConfig &config = TestDesignConfig())
{
  //analyse_filtre(h).afficher("Vérification design passe-bas");

  soit npts = 2048;
  soit [fr,xm] = frmag(h, npts);
  soit err_gain = 0.0f;

  soit a = analyse_LIT(h, false);

  si(config.fc_3dB_théo >= 0)
  {
    soit err_fc = abs(a.largeur_lp - config.fc_3dB_théo);
    msg("fc théo = {}, fc réelle = {} -> erreur = {}", config.fc_3dB_théo, a.largeur_lp, err_fc);
    tsd_assert_msg(err_fc < 1.0f / npts, "Fréquence de coupure à -3 dB invalide.");
  }

  si(config.fc_6dB_théo >= 0)
  {
    soit err_fc = abs(a.largeur_lp_6dB - config.fc_6dB_théo);
    msg("fc théo = {}, fc réelle = {} -> erreur = {}", config.fc_6dB_théo, a.largeur_lp_6dB, err_fc);
    tsd_assert_msg(err_fc < 1.0f / npts, "Fréquence de coupure -6 dB invalide.");
  }


  // TODO : infos d'après analyse_LIT...
  si constexpr(std::is_same<T, Vecf>())
  {
    entier n = h.rows();
    si(config.type == TypeFiltre::PASSE_BAS)
    {
      // Gain DC
      err_gain = h.somme() - 1;
    }
    sinon si(config.type == TypeFiltre::PASSE_HAUT)
    {
      // Gain Nyquist
      err_gain = (signyquist(n) * h).somme() - 1;
    }
  }



  soit gain_low = 1.0f, gain_high = 0.0f;
  si(config.type == TypeFiltre::PASSE_HAUT)
    std::swap(gain_low, gain_high);


  soit emax_bp = abs(xm.head(800) - gain_low).valeur_max(),
       emax_bc = abs(xm.tail(800) - gain_high).valeur_max();

  soit proto = Vecf::zeros(npts);
  si(config.type == TypeFiltre::PASSE_HAUT)
    proto.tail(800).setConstant(1);
  sinon
    proto.head(800).setConstant(1);

  soit err = proto - xm;

  si(tests_debug_actif)
  {
    Figures fs;
    soit f = fs.subplot();
    f.plot(mag2db(xm), "b-", "Réponse");
    f.plot(mag2db(proto), "g-", "Prototype");
    fs.subplot().plot(err, "r-", "Erreur");
    fs.afficher("Test conception filtre");
  }

  soit atten = -mag2db(emax_bc);

  msg("Err. max BP = {}, err. max BC = {} ({:.1f} dB), erreur de gain = {}.",
      emax_bp, emax_bc, atten, err_gain);

  si(abs(err_gain) > /*1e-3*/0.03)
    echec("Design passe-bas / passe-haut : gain DC / nyquist != 1 (erreur = {}).", err_gain);

  si(max(emax_bp, emax_bc) > config.err_max)
    echec("Gabarit non respecté.");

  si((config.atten_min > 0) && (atten < config.atten_min))
    echec("Atténuation insuffisante ({} dB minimum attendu).", config.atten_min);
}




template<typename T>
void test_design_passe_bas(const T &h, const TestDesignConfig &config)
{
  soit c2 = config;
  c2.type = TypeFiltre::PASSE_BAS;
  test_design(h, c2);
}



void test_filtrage_ola()
{
  msg_majeur("Test filtrage OLA (domaine fréquentiel)...");
  FiltreFFTConfig config;

  config.avec_fenetrage         = oui;
  config.dim_blocs_temporel     = 512;
  config.nb_zeros_min           = 0;

  config.traitement_freq = [&](Veccf &X)
    {
      entier N = 512;
      X.segment(N/32, (30*N)/32).setZero();
    };


  soit [ola, N] = filtre_fft(config);


  soit x = signal_test();
  soit y = filtre_par_bloc<cfloat>(ola, x, 1000);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x,"b-","x");
    f.subplot().plot(y,"g-","filtrage");
    f.afficher("Filtrage OLA freq");
  }
}

static void test_retard(entier d)
{
  msg(" ..retard = {}", d);
  soit h = Vecf::ones(d);
  soit x = Vecf::zeros(20);
  x(3) = 1;

  soit y  = filtrer(h, x);
  soit y2 = filtrer(h, y);

  soit index = y2.index_max();
  tsd_assert(index == 3+d-1);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x,"b-o","x");
    f.subplot().plot(y,"b-o","y");
    f.subplot().plot(y2,"b-o","y2");
    f.afficher(sformat("Test retard - filtre mg({}) - retard {}", d, (d-1)/2.0f));
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
  soit nc = 31;
  soit n  = nc + 50;
  soit h = linspace(1, nc, nc);

  // Note: pourrait être fait avec repimp,
  // mais ici le but est de tester filtre_rif
  soit x = sigimp(n);
  soit f = filtre_rif<float, float>(h);
  soit y = f->step(x);

  si(y.rows() != n)
    echec("n = {}, y.rows() = {}", n, y.rows());

  soit verr = y-vconcat(h, Vecf::zeros(n - nc));

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot(311).plot(h, "b|o", "Coefficients");
    f.subplot(312).plot(y, "g|o", "Réponse impulsionnelle");
    f.subplot(313).plot(verr, "r|o", "Erreur");
    f.afficher();
  }

  soit err = abs(verr).valeur_max();

  msg("err = {}", err);
  si(err > 1e-7)
    echec("Erreur trop importante.");
}


void test_rif_vs_rif_fft()
{
  msg_majeur("Tests rif vs rif fft...");

  soit h = design_rif_fen(127, "lp", 0.02);

  soit f1 = filtre_rif<float>(h);
  soit f2 = filtre_rif_fft<float>(h);

  soit x = signal_test();

  soit y1 = filtre_par_bloc(f1, x, 1000);
  soit y2 = filtre_par_bloc(f2, x, 1000);


  msg("y1 : {} - {}", y1.valeur_min(), y1.valeur_max());
  msg("y2 : {} - {}", y2.valeur_min(), y2.valeur_max());

  soit [y1p,y2p,i,score] = aligne_entier(y1, y2);

  soit err = (y1p - y2p).head(y1p.rows() - 4 * h.rows());
  soit errmax = abs(err).valeur_max();

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "b-", "x");
    f.subplot().plot(y1, "g-", "rif std");
    f.subplot().plot(y2, "b-", "rif fft");
    soit s = f.subplot();
    s.plot(y1p, "g-", "rif std");
    s.plot(y2p, "b-", "rif fft");
    s.titre("std vs fft (alignés)");
    f.subplot().plot(err, "r-", "Erreur (max = {})", errmax);
    f.afficher("Filtre RIF", {1400,800});
  }


  msg("Erreur max = {}", errmax);
  si(errmax > 1e-6)
    echec("Test RIF FFT : erreur trop importante ({})", errmax);
}

void test_filtre_rii()
{
  // Test simple avec un lissage exponentiel (facile de générer un vecteur de référence)
  //  y_n + b_1 y_{n-1} =  a_0 x_n + a_1 x_{n-1} + ...

  msg_majeur("Test filtre RII...");

  float α = 0.1;

  Vecf b(2), a(1);
  a(0) = α;
  b(0) = 1.0f;
  b(1) = -(1-α);

  soit h = FRat<float>::rii(a, b);

  msg("H = {}", h);
  msg("H(z^-1) = {}", h.eval_inv_z());

  soit f = filtre_rii<float,float>(h);

  entier n = 20;
  soit x = Vecf::ones(n);
  Vecf yref(n), y;

  y = f->step(x);

  yref(0) = α;
  pour(auto i = 1; i < n; i++)
    yref(i) = yref(i-1) + α * (1 - yref(i-1));


  si(tests_debug_actif)
  {
    Figures fig;
    soit fi = fig.subplot();
    fi.plot(x, "b-o", "x");
    fi.plot(y, "m-o", "y");
    fi.plot(yref, "g-s", "yref");
    fi = fig.subplot();
    fi.plot(y - yref, "r-o", "erreur");
    fig.afficher("Filtrage RII");
  }

  soit err = abs(yref - y).valeur_max();

  msg("Erreur max filtre RII : {}", err);

  si(err > 1e-6f)
    echec("Erreur trop importante");


  msg("ok.");
}




void test_rif_freq()
{
  msg_majeur("Test RIF freq...");
  // Nombre de coefficients souhaités
  entier n = 19;
  // Nombre de points du gabarit
  entier m = (n+1)/2; // m = 10

  soit d = Vecf::zeros(m);
  d.head(m/2).setConstant(1);

  soit h = design_rif_freq(n, d);
  tsd_assert(h.rows() == n);
  test_design(h, {TypeFiltre::PASSE_BAS, -1, 0.2});

  soit [fr,xm] = frmag<float>(h);
  soit fr1 = design_rif_freq_freqs(n);

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(fr,xm,"b-", "Réponse obtenue");

    f.plot(fr1, d, "gs", "Points d'échantillonnage (gabarit)");
    f.titre("Réponse fréquentielle");
    f.afficher();
  }


  // Vérification réponse fréquentielle = réponse attendue
  soit y = repfreq(h, fr1);

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(fr1, abs(d), "bo-", "Réponse attendue");
    f.plot(fr1, abs(y), "gs-", "Réponse réelle");
    f.plot(fr1, abs(abs(d)-abs(y)), "r-", "Erreur");
    f.afficher();
  }


  soit err = (abs(abs(y) - abs(d))).valeur_max();

  msg("Erreur sur les points d'interpolation : {}", err);

  // Vérification filtre de type I (phase linéaire)
  soit err_lp = abs(h.head(n/2) - h.tail(n/2).reverse()).valeur_max();

  msg("Erreur phase linéaire : {}", err_lp);

  tsd_assert_msg(err + err_lp < 1e-6, "Design RIF freq invalide.");

  msg("ok.");
}


static void test_riia()
{
  msg("Test designs RIIA...");
  pour(auto structure : {"ellip", "butt", "cheb1", "cheb2"})
  {
    msg_majeur("Test design RIIA {}...", structure);

    // Ondulations max : 0.1 dB dans la bande passante, ou -60 dB dans la bande coupée
    //soit h = design_riia(4, "lp", structure, 0.2, 3, 30);
    soit h = design_riia(12, "lp", structure, 0.25, 0.1, 60);
    test_design(h);
  }
  msg("ok");
}

static void test_design_lexp()
{
  msg("Test design lexp...");

  soit h = design_lexp(0.5);

  msg("h = {}", h);

  // h = z / (2z-1)
  tsd_assert(h.horner(1.0f) == 1.0f);

  msg("h(1)={}, h(2)={}", h.horner(1.0f), h.horner(2.0f));

  tsd_assert(h.horner(2.0f) == 2.0f/3.0f);

  msg("ok.");
}

entier test_filtres()
{
  msg("Test RIF FEN (PB)...");

  tsd::filtrage::debug_design = tests_debug_actif;


  test_design_lexp();

  test_design(design_rif_fen(31, "lp", 0.25),
      {.err_max     = 1,
       .fc_6dB_théo = 0.25});


  test_rif_vs_rif_fft();

  test_retard();
  test_design_rif_prod();




  test_rif_freq();

  test_decimateur();

  test_ligne_a_retard();

  test_filtre_rii();

  test_riia();

  {

    Vecf h;

    msg_majeur("Test des design RIF...");

    msg("Test RCS...");
    test_design(design_rif_rcs(127, 0.2, 0.25));


    msg("Test CS...");
    test_design(design_rif_cs(127, 0.2, 0.25));


    test_rcs1();
    test_rcs();

    msg("Test RIF FEN (PB)...");
    test_design(design_rif_fen(127, "lp", 0.25), {.fc_6dB_théo = 0.25});


    msg("Test RIF FEN (PH)...");
    test_design(design_rif_fen(127, "ph", 0.25), {TypeFiltre::PASSE_HAUT});


    msg("Test RIF FEN - kaiser...");
    test_design(design_rif_fen_kaiser("lp", 0.25, 60, 0.01));

    {
      msg("test RIF freq");
      soit n = 127, m = 512;
      Vecf d(m);
      d.head(m/2).setConstant(1.0f);
      d.tail(m/2).setConstant(0.0f);
      soit h = design_rif_freq(n, d);
      tsd_assert(h.rows() == n);
      test_design(h);
    }



    msg("Test RIF EQ...");
    // TODO: TEST TROP LONG
    pour(entier nc : {15, 63, 127, 255})//, 1023})
    {
      msg("  NC = {}", nc);
      test_design(design_rif_eq(nc, {{0, 0.2, 1}, {0.3, 0.5, 0}}));
    }
  }


  test_gaussien();


  soit h = design_rif_rcs(63, 0.1, 1.0/(2*8));

  si(tests_debug_actif)
  {
    Figure f;
    f.plot(h, "", "Réponse impulsionnelle");
    f.afficher();
  }

  test_filtre_rif();


  si(test_filtfilt())
    retourne -1;






  test_filtrage_ola();


  test_filtre_mg();

  test_design_biquad();
  retourne 0;
}
