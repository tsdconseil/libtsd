#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/vue/axes.hpp"
#include "tsd/vue.hpp"


using namespace tsd;
using namespace tsd::vue;


void test_unites_unit(const std::vector<double> &tics, const std::string &unit = "")
{
  std::string s1 = "";
  for(auto t: tics)
    s1 += fmt::format("{},", t);



  auto [expo, ndigits] = tsd::vue::unites::calc_expo_nb_chiffres_commun(tics, unit);

  msg("test_unites : vec = [{}], unit={}, expo={}, ndigits={}", s1, unit, expo, ndigits);

  for(auto t : tics)
  {
    auto s = unites::valeur_vers_chaine(t, unit, expo, ndigits);
    msg(" label : {}", s);
  }
}

int test_unites()
{
  test_unites_unit({-0.02, 0.03});
  test_unites_unit({-0.02, 0.03, 0});

  test_unites_unit({-0.001, 0.001});
  test_unites_unit({-0.001, 0.001, 0});

  test_unites_unit({0.001, 0.02});
  test_unites_unit({0.001, 0.02, 0});

  test_unites_unit({-1, 2, 5});
  test_unites_unit({2, 5, 0});

  test_unites_unit({2, 50, 0});
  test_unites_unit({2, 50, 50000});

  test_unites_unit({2, 50, 0.12});


  test_unites_unit({20e-2, 15e-2, 10e-2, 5e-2, 0, -5e-2, -10e-2});


  test_unites_unit({-4e6, -2e6, 0, 2e6, 4e6});
  test_unites_unit({-4e6, -2e6, 0, 2e6, 4e6}, "Hz");


  test_unites_unit({1e6, 0, 0.75e6}, "Hz");


  test_unites_unit({-1e6, 1e6, 0.5e6, 0}, "Hz");

  test_unites_unit({0,1000000,-1000000,2000000,-2000000,3000000,-3000000,4000000,-4000000}, "Hz");

  return 0;
}

int test_figure()
{
  Figure f;

  ArrayXf x = linspace(0,1,10);
  ArrayXf y = ArrayXf::Ones(10);

  f.plot(x, y);
  f.enregistrer("./test-log/ones.png");

  // Vérifie pas de planté si vecteur vide
  {
    msg("Test pas de planté si vecteurs vides...");

    Figure f;
    ArrayXf x;
    f.plot_psd(x);
    f.afficher();
    f.plot(x);

    msg("ok.");
  }

  {
    msg("Test échelle logarithmique...");
    Figure f;
    int n = 30;
    ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));
    f.subplot(211);
    f.plot(x, "-ob");
    f.titre("Echelle linéaire");
    f.subplot(212);
    auto c = f.plot(x, "-ob");
    f.gca().def_echelle("lin", "log");
    f.titre("Echelle logarithmique");
    f.enregistrer("./test-log/log.png");
  }

  {
      msg("Test plot min / max...");
      Figure f;
      int n = 30;
      ArrayXf t = linspace(0, n-1, n);
      ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));
      x.tail(2).setZero();
      f.subplot(211);
      auto c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
      c.def_couleur({128,128,255});
      f.plot(t, x, "k-o", "x");

      c = f.plot_minmax(t, x - 0.3 * x + 0.4, x + 0.3 * x + 0.4);
      c.def_couleur({255,128,128});
      f.plot(t, x + 0.4, "k-o", "y");


      f.titre("Echelle linéaire");
      f.subplot(212);
      c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
      c.def_couleur({128,128,255});
      f.plot(t, x, "k-o");
      f.gca().def_echelle("lin", "log");
      f.titre("Echelle logarithmique");
      f.afficher();
    }

  {
      msg("Test plot min / max (2)...");
      Figure f;
      int n = 30;
      ArrayXf t = linspace(0, n-1, n);
      ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));



      auto c = f.plot(t, x, "b-o", "x");
      c.def_sigma(x * 0.1);
      c = f.plot(t, x + 0.4, "g-o", "y");
      c.def_sigma(x * 0.2);
      f.titre("Echelle linéaire");
      f.afficher();
    }

  {
    msg("test BER...");

    int n = 10;
    ArrayXf ber(n), sigma(n);
    ber   << 0.452283,0.323602,0.189487,0.0862274,0.0218109,0.00124748,1.00604e-05,0,0,0;
    sigma << 0.0388947,0.0587929,0.00844648,0.00233075,0.0018141,0.000457055,3.01811e-05,0,0,0;

    Figure f;
    ArrayXf t = linspace(0, n-1, n);
    auto c = f.plot(t, ber, "b-o", "x");
    f.plot(t, ber/2, "m-s", "x/2");
    c.def_sigma(sigma);
    f.titre("BER");
    f.gca().def_echelle("lin", "log");
    f.afficher();
  }

  {
    msg("Test signal constant...");
    ArrayXf x = Eigen::ArrayXf::Constant(100, -2.24e38f);
    x += 1e30f * randn(100);
    Figure f;
    f.plot(x);
    f.afficher();
  }

  return 0;
}
