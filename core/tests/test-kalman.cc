#include "tsd/tsd-all.hpp"
#include "tsd/apps/kalman.hpp"
#include "tsd/tests.hpp"

using namespace tsd::kalman;


static void test_modele_ma()
{
  msg("test modele ma...");

  //! [ex_modele_ma]
  auto sys = modele_marche_aleatoire(2);

  sys->verifications();

  auto [x,y] = sys->steps(1000);

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x.row(0), x.row(1), "b.", "Vrais états");
    f.plot(y.row(0), y.row(1), "m.", "Observations");
  }

  //! [ex_modele_ma]
  //f.enregistrer("../doxy/images/ex-modele-ma.png");
}

static void test_modele_constante()
{
  msg("test modele cst...");

  //! [ex_modele_constante]
  auto sys = modele_constante();
  sys->verifications();


  auto [x,y] = sys->steps(100, 5 * VectorXf::Ones(1));

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x.row(0).transpose(), "b-", "Vrais états");
    f.plot(y.row(0).transpose(), "m*", "Observations");
  }

  //! [ex_modele_constante]
  //f.enregistrer("../doxy/images/ex-modele-constante.png");
}

/*ArrayXf logspace(float a, float b, int n)
{
  return (10 * ArrayXf::Ones(n)).pow(linspace(a, b, n));
}*/

static void test_kalman_ssg()
{
  msg("test kalman ssg...");

  //! [ex_kalman_ssg]
  auto ssm = modele_constante();
  ssm->verifications();
  int n = 50;
  ArrayXf Q = logspace(.5, -4, n);
  ArrayXf K = ArrayXf::Zero(n);
  for(auto i = 0; i < n; i++)
  {
    ssm->Q.resize(1,1);
    ssm->Q(0,0) = Q(i) * Q(i);
    K(i) = kalman_ssg(ssm)(0);
  }
  if(tests_debug_actif)
  {
    Figure f;
    f.plot(Q,K);
    f.titres("Steady-state Kalman gain vs process noise deviation", "Process noise (standard deviation relative to obs. noise)", "Steady-state gain");
  }

  //f.enregistrer("../doxy/images/ex-kalamn-ssg.png");

  //! [ex_kalman_ssg]
}

int test_kalman()
{
  msg_majeur("Tests Kalman...");
  test_modele_constante();
  test_modele_ma();
  test_kalman_ssg();
  msg_majeur("Fin.");
  return 0;
}


