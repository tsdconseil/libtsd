#include "tsd/tsd-all.hpp"

entier main()
{
  msg("Exemple simple libtsd.");

  soit n = 1000;
  soit t = linspace(0, 1, n),
       x = 2 * t - t * t,
       y = x + 0.1 * randn(n);

  soit h = design_rif_fen(31, "lp", 0.05);
  soit filtre = filtre_rif<float>(h);
  soit z = filtre->step(y);


  Figure f;
  f.plot(t, x, "b-", "x (signal test)");
  f.plot(t, y, "g-", "y (bruit)");
  f.plot(t, z, "r-", "z (filtrage rif)");
  f.enregistrer("./exemple-figure.png");

  msg("Fin.");
}


