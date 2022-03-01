#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/figure.hpp"
#include "tsd/tests.hpp"

using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

/*void test_fenetre(const std::string &nom, const ArrayXf &x)
{
  fenetre_analyse(nom, x);
}*/

void verifie_fenetre(const std::string &nom, const ArrayXf &x, const FenInfos &ref)
{
  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x, "", "Fenetre " + nom + " : vue temporelle");
    f.afficher();
  }

  if(x.hasNaN())
  {
    msg_erreur("Fenetre {} avec Nan.", nom);
  }
  //int n = x.rows();
  auto mes = fenetre_analyse(nom, x, tests_debug_actif);
  if(tests_debug_actif)
    mes.fig.afficher();
  //auto err1 = std::max(ref.atten_ls - mes.atten_pls, 0.0f);
  auto err1 = std::abs(ref.atten_ls - mes.atten_ls);
  msg("  erreur atten : {:.2f} dB", err1);
  auto err = std::abs(mes.largeur_lp - ref.largeur_lp) / ref.largeur_lp;
  msg("  erreur relative largeur lp : {:.2f} %", err*100);
  if((ref.atten_ls > 0) && (err1 > 1))
  {
    msg_erreur("Fenêtre {} : atténuation invalide ({:.1f} dB, {:.1f} dB attendu).", nom, mes.atten_ls, ref.atten_ls);
  }
  if((ref.largeur_lp > 0) && (err > 0.1))
  {
    msg_erreur("Fenêtre {} : largeur lobe principal invalide ({}, {} attendu, erreur relative = {} %).", nom, mes.largeur_lp, ref.largeur_lp, err*100);
  }
  if((ref.symetrique) && (!mes.symetrique))
  {
    msg_erreur("Fenêtre {} : symétrie attendue = {}, détectée = {}", nom, ref.symetrique, mes.symetrique);
  }
}

void test_kaiser_param(float As, float β_attendu)
{
  auto [β, n] = kaiser_param(As, 0.01);
  float err = std::abs(β - β_attendu);

  msg("test_kaiser_param : As={:.1f} dB, β attentu = {}, β calculé = {}, err={}.",
      As, β_attendu, β, err);
  //msg("  n calculé = {}", n);

  if(err > 0.1)
  {
    msg_erreur("Test kaiser param : erreur trop importante.");
    //throw;
  }
}

int test_fenetres()
{
  msg_majeur("Test des fenêtres...");

  int n = 128;

  auto nn = {127, 128};

  for(auto n : nn)
  {
    verifie_fenetre("re", fenetre("re", n, true), {13.32f, 0.88f/n, true});
    verifie_fenetre("hn", fenetre("hn", n, true), {31.47f, 1.44f/n, true});
    verifie_fenetre("hm", fenetre("hm", n, true), {42.68f, 1.30f/n, true});
    verifie_fenetre("tr", fenetre("tr", n, true), {26.50f, 1.27f/n, true});
    verifie_fenetre("bm", fenetre("bm", n, true), {-1, -1, false});
    verifie_fenetre("ch", fenetre("ch", n, true), {-1, -1, false});

    verifie_fenetre("re", fenetre("re", n, false), {13.32f, 0.88f/n, false});
    verifie_fenetre("hn", fenetre("hn", n, false), {31.47f, 1.44f/n, false});
    verifie_fenetre("hm", fenetre("hm", n, false), {42.68f, 1.30f/n, false});
    verifie_fenetre("tr", fenetre("tr", n, false), {26.50f, 1.27f/n, false});
    verifie_fenetre("bm", fenetre("bm", n, false), {-1, -1, false});
    verifie_fenetre("ch", fenetre("ch", n, false), {-1, -1, false});
  }

  verifie_fenetre("slepian", fenêtre_slepian(128, 0.1), {-1, -1});

  // Attention, ceci est faux : on aimerait plutôt spécifier une attén de filtre
  verifie_fenetre("cheby", fenêtre_chebychev(9,   60,  true), {-1, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  40,  true), {40, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  60,  true), {60, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  100, true), {100, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(128, 40,  true), {40, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(128, 60,  true), {60, -1});


  verifie_fenetre("cheby", fenêtre_chebychev(1024,     60, true), {60, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(8*1024,   60, true), {60, -1});

  // TODO
  // Ne fonctionne pas (1.6 dB d'erreur)
  //verifie_fenetre("cheby", fenêtre_chebychev(32*1024,  60, true), {60, -1});
  //verifie_fenetre("cheby", fenêtre_chebychev(128*1024, 60, true), {60, -1});
  //verifie_fenetre("cheby", fenêtre_chebychev(512*1024, 60, true), {60, -1});

  if(tests_debug_actif)
  {
    int N = 128*1024;
    auto f = fenêtre_chebychev(N, 60, true);
    ArrayXcf x1 = ArrayXcf::Ones(N) + ArrayXcf::Random(N);
    ArrayXcf x2 = x1 * f;
    Figures fig;
    fig.subplot().plot_psd(x1, 1.0, "", "psd x1");
    fig.subplot().plot_psd(x2, 1.0, "", "psd x2");
    fig.afficher("Test Cheby / N=128K");
  }




  // Note : imprécision ici
  //verifie_fenetre("cheby", fenetre_chebychev(128, 100, true), {100, -1});

  msg("Test kaiser params...");
   // [n,w,beta] = kaiserord([0,0.01],[1 0],10^(-30/20)) sous Octave

  test_kaiser_param(20, 0.0000f);
  test_kaiser_param(30, 2.1166/π);
  test_kaiser_param(40, 3.3953/π);
  test_kaiser_param(50, 4.5335/π);
  test_kaiser_param(60, 5.6533/π);
  test_kaiser_param(100, 10.061/π);

  // Spec directement du coefficient de forme
  msg("Test kaiser, β = {}", 1.25);
  verifie_fenetre("kaiser", fenêtre_kaiser1(n, 1.25, false), {29.47f, 1.19f/n});


  // Pourquoi commenté ?
  //msg_avert("Vérification")
  //verifie_fenetre("kaiser", fenetre_kaiser1(n, ),  {30, -1});
  /*verifie_fenetre("kaiser", fenetre_kaiser(30, 0.05),  {30, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(60, 0.05),  {60, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(50, 0.05),  {50, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(80, 0.05),  {80, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(100, 0.05), {100, -1});*/

  return 0;
}


