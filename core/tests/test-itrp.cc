#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


entier test_itrp_irreg()
{
  soit x = Vecf::valeurs({0, 1});
  soit y = x.clone();
  soit x2 = linspace(0, 1, 100);

  soit yi = interp(x, y, x2);
  soit yj = interp(x, y, x2, InterpOption::CSPLINE);

  // TODO : automatiser

  si(tests_debug_actif)
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
  retourne 0;

}



void test_itrp_retard()
{
  soit itrp   = itrp_sinc<cfloat>({15, 1024, 0.5, "hn"});
  soit delais = linspace(0, 1, 11);
  soit M      = 50;
  soit x      = siggauss(M) * sigchirp(0.001, 0.2, M, 'q');

  pour(auto i = 0; i < 11; i++)
  {
    soit h = itrp->coefs(delais(i));
    soit filtre = filtre_rif<float,float>(h);

    soit y1 = filtre->step(x);
    y1 = délais(y1, -7);

    soit y2 = délais(x, delais(i));

    si(tests_debug_actif)
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


entier test_itrp()
{
  msg_majeur("Test des interpolateurs...");

  soit itrps =
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

  pour(auto itrp: itrps)
  {
    msg("Test interpolateur [{:20s}] : nb points = {}, délais = {}.", itrp->nom, itrp->K, itrp->delais);

    Figure f, f2;

    pour(auto δ = 0.0f; δ <= 1.0f; δ += 0.1f)
    {
      soit cl = Couleur{255*δ,0,255*(1-δ)};
      soit h  = itrp->coefs(δ);
      soit c  = f.plot(h, "b-o", "Délais = {:.1f}", δ);
      c.def_couleur(cl);

      //msg("  frmag...");
      soit [fr,xm] = frmag(h);
      //msg("  ok.");
      c = f2.plot(fr, xm, "g-", "Délais = {:.1f}", δ);
      c.def_couleur(cl);
    }

    si(tests_debug_actif)
    {
      f.afficher(format("Interpolateur {} - réponse impulsionnelle", itrp->nom));
      f2.afficher(format("Interpolateur {} - réponse fréquentielle", itrp->nom));
    }
    //analyse_filtre(h, 1.0f, itrp->nom);
  }

  test_itrp_retard();
  retourne 0;
}
