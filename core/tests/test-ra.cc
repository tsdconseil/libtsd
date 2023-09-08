#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


struct ResultatPurete
{
  float freq_spurius,
        max_spurius_db;
};

static ResultatPurete verifie_sinus(cstring s, const Vecf x, float f)
{
  assertion((f >= 0) && (f <= 1));

  soit n = x.rows();

  soit fen  = fenêtre("hn", x.rows(), non),
       X    = abs2(fft(x * fen).head(x.rows()/2));

  soit [freqs, X_psd] = psd(x);

  soit etotal = X.somme();
  // 0.5 <=> n/2
  soit idx = (entier) (f * x.rows());
  soit ef     = X(idx);
  si(idx > 0)
    ef += X(idx-1);
  si(idx + 1 < n/2)
    ef += X(idx+1);

  si(tests_debug_actif)
  {
    Figure fig;
    fig.plot(X, "-b", "FFT");
    fig.plot((float) idx, ef, "or");
    fig.afficher("Vérif sinus");
    msg("f = {}, idx = {}, n = {}", f, idx, n);
  }

  soit score = ef / etotal;

  si(idx >= 10)
    X.segment(idx-10, 20).setConstant(0);
  entier idx_spur;
  ResultatPurete res;

  res.max_spurius_db  = 10 * log10(X.maxCoeff(&idx_spur) / ef);
  res.freq_spurius    = ((float) idx_spur) / x.rows();

  msg("{} : check fréq = {}, ratio énergie = {}, max spurius = {:.1f} dB (@ f = {} * fe)",
      s, f, score, res.max_spurius_db, res.freq_spurius);
  si(std::isnan(score) || (score < 0.8))
    échec(" Signal invalide (une sinusoide pure est attendue).");
  retourne res;
}

static void test_ra_unit(cstring nom, float ratio, sptr<FiltreGen<float>> ra, float max_spurius_dB = -50)
{
  msg("Test adaptation de rythme [{}] : ratio = {}", nom, ratio);
  soit fe = 100e3f,
       f1 = 2e6f,
       f2 = 2e3f;

  msg("Fe = {} Hz, f1 = {} Hz, f2 = {} Hz.", fe, f1, f2);
  msg("Ratio = {} --> fe' = {}", ratio, fe * ratio);

  soit t = intervalle_temporel(1000, fe),
       x = sin(t*2*π*f2),
       y = ra->step(x);

  verifie_sinus(nom + " - signal d'entrée", x, f2/fe);
  soit spy = verifie_sinus(nom + " - signal de sortie", y, f2/(ratio*fe));

  si(tests_debug_actif)
  {
    Figures fs;
    soit f = fs.subplot(221);
    f.plot(x);
    f.titres("Entrée", "Echantillons");
    fs.subplot(222).plot_psd(x, fe);
    f = fs.subplot(223);
    f.plot(y);
    f.titres("Sortie", "Echantillons");
    f = fs.subplot(224);

    soit [freqs, Y] = psd(y);
    msg("freqs: {}, Y : {}", freqs.rows(), Y.rows());

    f.plot_psd(y, fe*ratio);

    //ArrayXf freq = linspace(0, fe*ratio/2, Y.rows());
    f.plot(freqs*(fe*ratio), Y, "-g", "avec psd()");

    /*ArrayXf x1(1), y1(1);
    x1(0) = spy.freq_spurius * fe;
    y1(0) = Y(spy.freq_spurius / fe);
    f.plot(x1,y1,"ro");*/

    entier id = spy.freq_spurius * 2 * Y.rows();

    si(id < Y.rows())
    {
      f.canva().set_couleur(tsd::vue::Couleur::Rouge);
      f.canva().marqueur({(spy.freq_spurius * fe * ratio), Y(id)}, tsd::vue::Marqueur::CERCLE, 7);
    }

    fs.afficher(format("test ra [{}], ratio = {}", nom, ratio));

    //Figure f2;
    //f2.plot_psd(x, fe);
    //f2.afficher("F2");

  }
  /*{
    Figure f(format("test ra, ratio = {}, zoom", ratio));
    f.subplot(211);
    f.plot(x.head(50), "b-o", "x");
    f.subplot(212);
    f.plot(y.head(50*ratio), "g-o", "y");
    f.afficher();
  }*/

  soit err = 100.0 * abs((y.rows() - ratio * x.rows())/x.rows());

  //msg(" nb ech entrée = {}, sortie = {}, err = {} %", x.rows(), y.rows(), 100 * err);

  // TODO : test resample
  //resample(x, lom)

  assertion_msg(err < 1, "ra : nombre d'échantillons invalide.");

  soit amp1 = x.valeur_max() - x.valeur_min(),
       amp2 = y.valeur_max() - y.valeur_min(),
       erra = 100 * (amp1 - amp2) / amp1,
       rms1 = sqrt(square(x).moyenne()),
       rms2 = sqrt(square(y).moyenne()),
       errr = 100 * (rms1 - rms2) / rms1;

  msg("amp1 = {}, amp2 = {}, err = {} %", amp1, amp2, erra);
  msg("rms1 = {}, rms2 = {}, err = {} %", rms1, rms2, errr);

  assertion_msg(erra < /*0.1*/10, "ra [{}] : problème d'amplitude.", nom);
  assertion_msg(spy.max_spurius_db <= max_spurius_dB, "ra [{}] : trops de spurius ({:.1f} dB)", nom, spy.max_spurius_db);

  msg_majeur("Tests RA ok.");
}



static void test_ra_unit(float ratio)
{
  msg_majeur("Test ratio = {}, filtre_itrp (cspline)", ratio);
  soit itrp = itrp_cspline<float>();
  soit ra = filtre_itrp<float>(ratio, itrp);
  test_ra_unit("filtre_itrp/cspline", ratio, ra);

  msg_majeur("Test ratio = {}, filtre_itrp (sinc)", ratio);
  itrp = itrp_sinc<float>({127, 256, 0.5, "hn"});
  ra = filtre_itrp<float>(ratio, itrp);
  test_ra_unit("filtre_itrp/sinc", ratio, ra);

  msg_majeur("Test ratio = {}, filtre_reechan", ratio);
  ra = filtre_reechan<float>(ratio);
  test_ra_unit("filtre_reechan", ratio, ra);
}



static void test_filtre_rif_demi_bande()
{
  msg_majeur("Test filtre_rif_demi_bande");
  soit h = design_rif_fen(15, "lp", 0.25, "hn");
  soit ra = filtre_rif_demi_bande<float,float>(h);
  test_ra_unit("rif demi-bande", 0.5, ra);
}


static void test_filtre_rif_decim(entier R)
{
  msg_majeur("Test filtre_rif_decim : R = {}", R);
  soit h = design_rif_fen(15, "lp", 0.5/R, "hn");

  msg("Somme(h) = {}", h.somme());

  soit ra = filtre_rif_decim<float,float>(h, R);
  test_ra_unit("rif decim", 1.0/R, ra);
}

static void test_filtre_rif_decim()
{
  pour(auto R: {2, 3, 4, 5, 8})
    test_filtre_rif_decim(R);
}

void test_filtre_rif_ups()
{
  msg_majeur("Test filtre_rif_ups");
  soit h  = design_rif_fen(15, "lp", 0.25, "hn");
  soit ra = filtre_rif_ups<float,float>(h, 2);
  test_ra_unit("rif ups", 2.0, ra);
}





void test_ra()
{
  test_filtre_rif_decim();
  test_filtre_rif_ups();

  test_filtre_rif_demi_bande();

  //soit ratios = ;
  pour(auto ratio: {1.f, 1.5f, 0.5f, 2.f, 1.2f, π_f})
    test_ra_unit(ratio);
}


