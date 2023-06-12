#include "tsd/tsd-all.hpp"
#include "tsd/fourier/tod.hpp"
#include "tsd/tests.hpp"

using namespace tsd::tf::tod;


void infos_ondelette(sptr<Ondelette<float>> ondelette)
{
  msg("Test de l'ondelette [{}]...", ondelette->nom);

  // construction de l'ondelette mère
  // interprétation comme un banc de filtre (et construction des deux filtres,
  // de leurs réponses fréquentielle et impulsionnelle)

  {
    entier n = 512;
    //ArrayXf x, x0 = ArrayXf::Zero(n);
    //x0(0) = 1;
    Vecf x, x0;
    x0 = randn(n);
    x = x0;
    dwt(ondelette, x, 9);
    iwt(ondelette, x, 9);

    soit errmax = abs(x-x0).valeur_max();
    soit err = sqrt(abs2(x-x0).moyenne());

    msg("Erreur RMS reconstruction : {}, erreur max = {}", err, errmax);
  }



  msg("  construction de l'ondelette mère...");


  entier n = 512;
  soit x = Vecf::zeros(n);
  x(0) = 1;

  iwt(ondelette, x, 9);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x);
    x.setZero();
    x(60) = 1;
    iwt(ondelette, x, 9);
    f.subplot().plot(x);
    f.afficher("Ondelette - " + ondelette->nom);
  }
}

void analyse_lift(const Lift &lift)
{
  msg_majeur("Analyse du schéma de lifting [{}]...", lift.nom);


  FormePolyphase fpol(lift);

  std::cout << fpol;


  QMF qmf(fpol);

  std::cout << qmf;

  soit h0 = qmf.H0.coefs, h1 = qmf.H1.coefs;

  //analyse_filtre(h0, 1, "h0");
  //analyse_filtre(h1, 1, "h1");

  si(tests_debug_actif)
  {
    Figure f;

    soit [fr0,xm0] = frmag(h0);
    soit [fr1,xm1] = frmag(h1);

    f.plot(fr0, xm0, "b-", "h0");
    f.plot(fr1, xm1, "g-", "h1");
    f.afficher(lift.nom + " - h0 vs h1");
  }

  infos_ondelette(ondelette_gen<float>(lift));
}


void test_tod()
{
  //infos_ondelette(ondelette_db4<float>());
  //infos_ondelette(ondelette_haar<float>());

  msg("Lift haar...");
  analyse_lift(lift_haar());

  msg("Lift db2...");
  analyse_lift(lift_db2());
}
