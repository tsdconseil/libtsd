#include "tsd/tsd.hpp"
#include "tsd/stats.hpp"

using namespace tsd;
using namespace tsd::stats;

void test_levinson()
{
  msg_majeur("Test levinson - résolution R a = -r...");

  int p = 5;
  ArrayXf r = ArrayXf::Random(p);

  Eigen::MatrixXf R = r_vers_R(r.head(p-1));

  //msg("r =\n{}", r);
  //msg("R =\n{}", R);

  ArrayXf ac = levinson_reel(r);
  ArrayXf ar = R.lu().solve(-r.tail(p-1).matrix()).array();

  ArrayXf ac2 = levinson(r.head(p-1), r.head(p-1), -r.tail(p-1));

  tsd_assert(ar.rows() == p - 1);
  tsd_assert(ac.rows() == p);
  tsd_assert(ac2.rows() == p - 1);

  //msg("ac =\n{}", ac);
  //msg("ar =\n{}", ar);

  float err  = (ac.tail(p-1) - ar).abs().maxCoeff();
  float erry = ((R * ac.tail(p-1).matrix()).array() + r.tail(p-1)).abs().maxCoeff();
  float erryref = ((R * ar.matrix()).array() + r.tail(p-1)).abs().maxCoeff();
  float erry2 = ((R * ac2.matrix()).array() + r.tail(p-1)).abs().maxCoeff();

  msg("Erreur solution : {}, prédiction LU : {}, inversion levin : {}, inversion levin gen : {}", err, erryref, erry, erry2);

  if((err > 1e-5) || (erry > 1e-6) || (erryref > 1e-6) || (erry2 > 1e-6))
    echec("Echec levinson.");


  msg("ok.");
}

int test_stats()
{
  test_levinson();
  return 0;
}


