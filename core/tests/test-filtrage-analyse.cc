#include "tsd/tsd.hpp"
#include "tsd/tests.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"

using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::vue;


static void test_frmag()
{
  msg("Test frmag...");
  int n = 100;
  ArrayXf x = randn(n);
  auto [fr, mag] = frmag(x, 233);
  int m = fr.rows();
  tsd_assert(m == 233);
  tsd_assert(m == mag.rows());

  ArrayXf vmag(m);
  for(auto i = 0; i < m; i++)
  {
    // Théoriquement :
    // mag(Omega) = somme h_k e^{i Omega k}
    vmag(i) = abs((x * tsd::polar(fr(i) * 2 * π * linspace(0, n-1, n))).sum());
  }
  if(tests_debug_actif)
  {
    Figure f;
    f.plot(fr, mag,  "b-", "mag");
    f.plot(fr, vmag, "g-", "vmag");
    f.afficher();
    f.clear();
    f.plot(fr, mag - vmag, "r-", "erreur");
    f.afficher();
  }
  auto err = (vmag - mag).abs().maxCoeff();
  msg("Erreur max frmag = {}", err);
  tsd_assert(err < 1e-3);
}

int test_filtrage_analyse()
{
  test_frmag();
  return 0;
}

