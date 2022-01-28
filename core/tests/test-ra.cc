#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/tests.hpp"

using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

struct ResultatPurete
{
  float freq_spurius;
  float max_spurius_db;
  //float max_spurius_db_abs;
};

static ResultatPurete verifie_sinus(const std::string &s, const ArrayXf x, float f)
{
  tsd_assert((f >= 0) && (f <= 1));

  //ArrayXf X = (psd(x) * (log(10.0f)/10)).exp();


  //ArrayXf X    = rfft(x).head(x.rows()/2).abs2();

  ArrayXf fen  = tsd::filtrage::fenetre("hn", x.rows(), false);
  ArrayXf X    = fft(x * fen).head(x.rows()/2).abs2();

  auto [freqs, X_psd] = psd(x);

  // Ca ne marche pas (le produit par la fenêtre)
  //ArrayXf X    = fft(x * fenetre("hn", x.rows())).head(x.rows()/2).abs2();

/*
# if 0

# else

  ArrayXcf yc   = x;
  ArrayXf fen   = tsd::filtrage::fenetre("hn", x.rows(), false);
  ArrayXcf ycf  = yc * fen;
  ArrayXf X     = ((tsd::fourier::fft(ycf)).abs2()).head(x.rows() / 2);
# endif
*/
  float etotal = X.sum();
  // 0.5 <=> n/2
  int idx = f * x.rows();
  float ef     = X(idx) + X(idx-1) + X(idx+1);
  float score = ef / etotal;


  //X = 10*log10(X).eval();
  /*{
    Figure f;
    f.plot(10*log10(X));
    f.plot(X_psd, "g-", "X PSD");
    f.titre("verifie sinus");
    f.afficher();
  }*/

  X.segment(idx-10, 20).setZero();
  int idx_spur;
  ResultatPurete res;



  res.max_spurius_db  = 10 * log10(X.maxCoeff(&idx_spur) / ef);
  res.freq_spurius    = ((float) idx_spur) / x.rows();

  //float max_spurius = std::max(X.head(idx-10).maxCoeff(), X.tail(X.rows() - (idx+10)).maxCoeff());
  //max_spurius = 10 * log10(max_spurius / ef);

  msg("{} : check fréq = {}, ratio énergie = {}, max spurius = {:.1f} dB (@ f = {} * fe)",
      s, f, score, res.max_spurius_db, res.freq_spurius);
  if(std::isnan(score) || (score < 0.8))
    echec(" Signal invalide (une sinusoide pure est attendue).");
  return res;
}

static void test_ra_unit(const std::string &nom, float ratio, sptr<FiltreGen<float>> ra, float max_spurius_dB = -50)
{
  msg("Test adaptation de rythme [{}] : ratio = {}", nom, ratio);
  auto fe = 100e3;
  auto f1 = 2e6f, f2 = 2e3f;

  msg("Fe = {} Hz, f1 = {} Hz, f2 = {} Hz.", fe, f1, f2);
  msg("Ratio = {} --> fe' = {}", ratio, fe * ratio);

  ArrayXf t = trange(1000, fe);
  ArrayXf x = (t*2*π*f2).sin();

  verifie_sinus(nom + " - signal d'entrée", x, f2/fe);

  ArrayXf y = ra->step(x);

  auto spy = verifie_sinus(nom + " - signal de sortie", y, f2/(ratio*fe));

  if(tests_debug_actif)
  {
    Figures fs;
    auto f = fs.subplot(221);
    f.plot(x);
    f.titres("Entrée", "Echantillons");
    f = fs.subplot(222);
    f.plot_psd(x, fe);
    f = fs.subplot(223);
    f.plot(y);
    f.titres("Sortie", "Echantillons");
    f = fs.subplot(224);

    auto [freqs, Y] = psd(y);
    msg("freqs: {}, Y : {}", freqs.rows(), Y.rows());

    f.plot_psd(y, fe*ratio);

    //ArrayXf freq = linspace(0, fe*ratio/2, Y.rows());
    f.plot(freqs*(fe*ratio), Y, "-g", "avec psd()");

    /*ArrayXf x1(1), y1(1);
    x1(0) = spy.freq_spurius * fe;
    y1(0) = Y(spy.freq_spurius / fe);
    f.plot(x1,y1,"ro");*/

    int id = spy.freq_spurius * 2*Y.rows();

    if(id < Y.rows())
    {
      f.canva().set_couleur(tsd::vue::Couleur::Rouge);
      f.canva().marqueur({(float) (spy.freq_spurius * fe * ratio), Y(id)}, tsd::vue::Marqueur::CERCLE, 7);
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

  auto err = 100.0 * std::abs((y.rows() - ratio * x.rows())/x.rows());

  //msg(" nb ech entrée = {}, sortie = {}, err = {} %", x.rows(), y.rows(), 100 * err);

  // TODO : test resample
  //resample(x, lom)

  tsd_assert_msg(err < 1, "ra : nombre d'échantillons invalide.");

  float amp1 = x.maxCoeff() - x.minCoeff();
  float amp2 = y.maxCoeff() - y.minCoeff();
  float erra = 100 * (amp1 - amp2) / amp1;


  float rms1 = sqrt(x.square().mean());
  float rms2 = sqrt(y.square().mean());
  float errr = 100 * (rms1 - rms2) / rms1;

  msg("amp1 = {}, amp2 = {}, err = {} %", amp1, amp2, erra);
  msg("rms1 = {}, rms2 = {}, err = {} %", rms1, rms2, errr);

  tsd_assert_msg(erra < /*0.1*/10, "ra [{}] : problème d'amplitude.", nom);


  tsd_assert_msg(spy.max_spurius_db <= max_spurius_dB, "ra [{}] : trops de spurius ({:.1f} dB)", nom, spy.max_spurius_db);


  //auto ra2 = filtre_reechan<float>(1.0/ratio);
  //ArrayXf z =

  msg_majeur("Tests RA ok.");
}



static void test_ra_unit(float ratio)
{
  //auto itrp = itrp_sinc<float>(15, 0.25, "hn");
  //auto itrp = itrp_sinc<float>(127, 0.25, "hn");

  msg_majeur("Test ratio = {}, filtre_itrp (cspline)", ratio);
  auto itrp = itrp_cspline<float>();
  auto ra = filtre_itrp<float>(ratio, itrp);
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
  ArrayXf h = design_rif_fen(15, "lp", 0.25, "hn");
  auto ra = filtre_rif_demi_bande<float,float>(h);
  test_ra_unit("rif demi-bande", 0.5, ra);
}


static void test_filtre_rif_decim(int R)
{
  msg_majeur("Test filtre_rif_decim : R = {}", R);
  ArrayXf h = design_rif_fen(15, "lp", 0.5/R, "hn");

  msg("Somme(h) = {}", h.sum());


  auto ra = filtre_rif_decim<float,float>(h, R);
  test_ra_unit("rif decim", 1.0/R, ra);
}

static void test_filtre_rif_decim()
{
  for(auto R: {2, 3, 4, 5, 8})
    test_filtre_rif_decim(R);
}

void test_filtre_rif_ups()
{
  msg_majeur("Test filtre_rif_ups");
  ArrayXf h = design_rif_fen(15, "lp", 0.25, "hn");
  auto ra = filtre_rif_ups<float,float>(h, 2);
  test_ra_unit("rif ups", 2.0, ra);
}





int test_ra()
{






  test_filtre_rif_decim();
  test_filtre_rif_ups();

  test_filtre_rif_demi_bande();

  //auto ratios = ;
  for(auto ratio: {1.f, 1.5f, 0.5f, 2.f, 1.2f, π_f})
    test_ra_unit(ratio);

  return 0;
}


