#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


void verifie_fenetre(const std::string &nom, const Vecf &x, const FenInfos &ref)
{
  si(tests_debug_actif)
  {
    Figure f;
    f.plot(x, "", "Fenetre " + nom + " : vue temporelle");
    f.afficher();
  }

  si(x.hasNaN())
  {
    msg_erreur("Fenetre {} avec Nan.", nom);
  }
  //entier n = x.rows();
  soit mes = fenetre_analyse(nom, x, tests_debug_actif);
  si(tests_debug_actif)
    mes.fig.afficher();
  //soit err1 = max(ref.atten_ls - mes.atten_pls, 0.0f);
  soit err1 = abs(ref.atten_ls - mes.atten_ls);
  msg("  erreur atten : {:.2f} dB", err1);
  soit err = abs(mes.largeur_lp - ref.largeur_lp) / ref.largeur_lp;
  msg("  erreur relative largeur lp : {:.2f} %", err*100);
  si((ref.atten_ls > 0) && (err1 > 1))
  {
    msg_erreur("Fenêtre {} : atténuation invalide ({:.1f} dB, {:.1f} dB attendu).", nom, mes.atten_ls, ref.atten_ls);
  }
  si((ref.largeur_lp > 0) && (err > 0.1))
  {
    msg_erreur("Fenêtre {} : largeur lobe principal invalide ({}, {} attendu, erreur relative = {} %).", nom, mes.largeur_lp, ref.largeur_lp, err*100);
  }
  si((ref.symetrique) && (!mes.symetrique))
  {
    msg_erreur("Fenêtre {} : symétrie attendue = {}, détectée = {}", nom, ref.symetrique, mes.symetrique);
  }
}

void test_kaiser_param(float As, float β_attendu)
{
  soit [β, n] = kaiser_param(As, 0.01);
  float err = abs(β - β_attendu);

  msg("test_kaiser_param : As={:.1f} dB, β attentu = {}, β calculé = {}, err={}.",
      As, β_attendu, β, err);
  //msg("  n calculé = {}", n);

  si(err > 0.1)
  {
    msg_erreur("Test kaiser param : erreur trop importante.");
    //throw;
  }
}

entier test_fenetres()
{
  msg_majeur("Test des fenêtres...");


  {
    // Vérification petites fenêtres
    pour(auto n: {1, 2, 3, 4, 5})
    {
      verifie_fenetre("hn", fenetre("hn", n, oui), {0, 0, oui});
      verifie_fenetre("hn", fenetre("hn", n, non), {0, 0, non});
    }
  }


  soit nn = {127, 128};

  pour(auto n : nn)
  {
    verifie_fenetre("re", fenetre("re", n, oui), {13.32f, 0.88f/n, oui});
    verifie_fenetre("hn", fenetre("hn", n, oui), {31.47f, 1.44f/n, oui});
    verifie_fenetre("hm", fenetre("hm", n, oui), {42.68f, 1.30f/n, oui});
    verifie_fenetre("tr", fenetre("tr", n, oui), {26.50f, 1.27f/n, oui});
    verifie_fenetre("bm", fenetre("bm", n, oui), {-1, -1, non});
    verifie_fenetre("ch", fenetre("ch", n, oui), {-1, -1, non});

    verifie_fenetre("re", fenetre("re", n, non), {13.32f, 0.88f/n, non});
    verifie_fenetre("hn", fenetre("hn", n, non), {31.47f, 1.44f/n, non});
    verifie_fenetre("hm", fenetre("hm", n, non), {42.68f, 1.30f/n, non});
    verifie_fenetre("tr", fenetre("tr", n, non), {26.50f, 1.27f/n, non});
    verifie_fenetre("bm", fenetre("bm", n, non), {-1, -1, non});
    verifie_fenetre("ch", fenetre("ch", n, non), {-1, -1, non});
  }

  verifie_fenetre("slepian", fenêtre_slepian(128, 0.1), {-1, -1});

  // Attention, ceci est faux : on aimerait plutôt spécifier une attén de filtre
  verifie_fenetre("cheby", fenêtre_chebychev(9,   60,  oui), {-1, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  40,  oui), {40, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  60,  oui), {60, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(32,  100, oui), {100, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(128, 40,  oui), {40, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(128, 60,  oui), {60, -1});


  verifie_fenetre("cheby", fenêtre_chebychev(1024,     60, oui), {60, -1});
  verifie_fenetre("cheby", fenêtre_chebychev(8*1024,   60, oui), {60, -1});

  // TODO
  // Ne fonctionne pas (1.6 dB d'erreur)
  //verifie_fenetre("cheby", fenêtre_chebychev(32*1024,  60, oui), {60, -1});
  //verifie_fenetre("cheby", fenêtre_chebychev(128*1024, 60, oui), {60, -1});
  //verifie_fenetre("cheby", fenêtre_chebychev(512*1024, 60, oui), {60, -1});

  si(tests_debug_actif)
  {
    soit N = 128*1024;
    soit f = fenêtre_chebychev(N, 60, oui);
    soit x1 = Veccf::ones(N) + randu(N);
    soit x2 = x1 * f;
    Figures fig;
    fig.subplot().plot_psd(x1, 1.0, "", "psd x1");
    fig.subplot().plot_psd(x2, 1.0, "", "psd x2");
    fig.afficher("Test Cheby / N=128K");
  }




  // Note : imprécision ici
  //verifie_fenetre("cheby", fenetre_chebychev(128, 100, oui), {100, -1});

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
  verifie_fenetre("kaiser", fenêtre_kaiser1(128, 1.25, non), {29.47f, 1.19f/128});


  // Pourquoi commenté ?
  //msg_avert("Vérification")
  //verifie_fenetre("kaiser", fenetre_kaiser1(n, ),  {30, -1});
  /*verifie_fenetre("kaiser", fenetre_kaiser(30, 0.05),  {30, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(60, 0.05),  {60, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(50, 0.05),  {50, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(80, 0.05),  {80, -1});
  verifie_fenetre("kaiser", fenetre_kaiser(100, 0.05), {100, -1});*/

  retourne 0;
}


