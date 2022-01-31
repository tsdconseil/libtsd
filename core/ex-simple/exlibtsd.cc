#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"

using namespace tsd;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

int main()
{
  msg("Exemple simple libtsd.");

  int n     = 1000;
  ArrayXf t = linspace(0, 1, n);
  ArrayXf x = 2 * t - t.pow(3);
  ArrayXf y = x + 0.1 * randn(n);

  ArrayXf h = design_rif_fen(31, "lp", 0.05);
  auto filtre = filtre_rif<float>(h);
  ArrayXf z = filtre->step(y);


  Figure f;
  f.plot(t, x, "b-", "x (signal test)");
  f.plot(t, y, "g-", "y (bruit)");
  f.plot(t, z, "r-", "z (filtrage rif)");
  f.enregistrer("./exemple-figure.png");

  msg("Fin.");
}


