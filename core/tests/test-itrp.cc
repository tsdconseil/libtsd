#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/tests.hpp"


using namespace tsd;
using namespace tsd::telecom;
using namespace tsd::filtrage;
using namespace tsd::vue;



int test_itrp_irreg()
{
  ArrayXf x(2), y(2);
  x << 0, 1;
  y << 0, 1;

  ArrayXf x2 = linspace(0, 1, 100);

  ArrayXf yi = interp(x, y, x2);
  ArrayXf yj = interp(x, y, x2, InterpOption::CSPLINE);

  // TODO : automatiser

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, y, "ob", "Points connus");
    f.gcf().plot(x2, yi, "g-", "Interpolation");
    f.gcf().titre("Interpolation linéaire");
    f.subplot().plot(x, y, "ob", "Points connus");
    f.gcf().plot(x2, yj, "g-", "Interpolation");
    f.gcf().titre("Interpolation par splines naturelles");
    f.afficher();
  }
  return 0;

}



void test_itrp_retard()
{
  auto itrp = tsd::filtrage::itrp_sinc<cfloat>({15, 1024, 0.5, "hn"});
  ArrayXf delais = linspace(0, 1, 11);

  int M = 50;
  ArrayXf x = siggauss(M) * sigchirp2(0.001, 0.2, M);

  for(auto i = 0; i < 11; i++)
  {
    ArrayXf h = itrp->coefs(delais(i));
    auto filtre = filtre_rif<float,float>(h);

    ArrayXf y1 = filtre->step(x);
    y1 = tsd::fourier::delais(y1, -7);

    ArrayXf y2 = tsd::fourier::delais(x, delais(i));

    if(tests_debug_actif)
    {
      Figures fs;
      fs.subplot().plot(x, "", "x");
      fs.subplot().plot(y1, "g-", "délais RIF");
      fs.subplot().plot(y2, "m-", "délais FFT");
      fs.subplot().plot(y2-y1, "r-", "Différence");
      fs.subplot().plot(y2-x, "r-", "Différence sans corr");
      fs.afficher(fmt::format("Test retard = {}", delais(i)));
    }
  }
}


int test_itrp()
{
  msg_majeur("Test des interpolateurs...");

  auto itrps =
    {
        itrp_lineaire<cfloat>(),
        itrp_cspline<cfloat>(),
        itrp_lagrange<cfloat>(1),
        itrp_lagrange<cfloat>(2),
        itrp_lagrange<cfloat>(3),
        itrp_lagrange<cfloat>(5),
        itrp_lagrange<cfloat>(6),
        itrp_lagrange<cfloat>(7),
        itrp_sinc<cfloat>({15, 256, 0.5, "re"}),
        itrp_sinc<cfloat>({15, 256, 0.25, "re"}),
        itrp_sinc<cfloat>({15, 256, 0.5, "hn"}),
        itrp_sinc<cfloat>({15, 256, 0.25, "hn"}),
        itrp_sinc<cfloat>({15, 256, 0.4, "hn"}),
        itrp_sinc<cfloat>({15, 512, 0.5, "hn"}),

        itrp_sinc<cfloat>({31, 256, 0.5, "hn"}),
        itrp_sinc<cfloat>({63, 256, 0.5, "hn"}),
        itrp_sinc<cfloat>({127, 256, 0.5, "hn"}),


        itrp_sinc<cfloat>({31, 256, 0.48, "hn"}),
        itrp_sinc<cfloat>({63, 256, 0.48, "hn"}),
    };

  for(auto itrp: itrps)
  {
    msg("Test interpolateur [{:20s}] : nb points = {}, délais = {}.", itrp->nom, itrp->K, itrp->delais);

    Figure f, f2;

    for(auto δ = 0.0f; δ <= 1.0f; δ += 0.1f)
    {
      Couleur cl{255*δ,0,255*(1-δ)};
      ArrayXf h = itrp->coefs(δ);
      auto c = f.plot(h, "b-o", "Délais = {:.1f}", δ);
      c.def_couleur(cl);

      //msg("  frmag...");
      auto [fr,xm] = frmag(h);
      //msg("  ok.");
      c = f2.plot(fr, xm, "g-", "Délais = {:.1f}", δ);
      c.def_couleur(cl);
    }

    if(tests_debug_actif)
    {
      f.afficher(format("Interpolateur {} - réponse impulsionnelle", itrp->nom));
      f2.afficher(format("Interpolateur {} - réponse fréquentielle", itrp->nom));
    }
    //analyse_filtre(h, 1.0f, itrp->nom);
  }

  test_itrp_retard();
  return 0;
}
