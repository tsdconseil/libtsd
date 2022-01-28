#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/figure.hpp"


using namespace tsd;
using namespace tsd::vue;


static void test_unites_unit(const std::vector<double> &tics, const std::string &unit = "")
{
  std::string s1 = "";
  for(auto t: tics)
    s1 += fmt::format("{},", t);



  auto [expo, ndigits] = tsd::vue::unites::calc_expo_nb_chiffres_commun(tics, unit);

  msg("test_unites : vec = [{}], unit={}, expo={}, ndigits={}", s1, unit, expo, ndigits);

  for(auto t : tics)
  {
    auto s = tsd::vue::unites::valeur_vers_chaine(t, unit, expo, ndigits);
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
  // Test affichage long
  {
    Figure f;
    int n = 15000;
    ArrayXf x = 0.1 * randn(n);
    for(auto i = 1000; i < 15000; i += 1000)
      x(i) = 1;
    f.plot(x);
    f.afficher("Test affichage 15000 points");
  }



  // Test échelle
  {
    Figure f;
    ArrayXf vt(180);
    for(auto i = 0; i < 180; i++)
    {
      vt(i) = i / 60.0f;
    }
    f.plot(vt, vt, "", "Temps (minutes)");
    f.afficher("Test échelle");
  }



  {
    Figure f;

    ArrayXf x = linspace(0,1,100);
    //ArrayXf y = ArrayXf::Ones(10);

    f.plot(x, x, "", "valeurs de x");
    f.plot(x, x.square(), "", "valeurs de x^2");
    //f.enregistrer("./test-log/ones.png");

    f.titres("Titre principal", "Axe x", "Axe y");

    f.canva().texte(0.1, 0.2, "Texte @ 0,1x0,2", 0.2, 0.1);
    f.canva().set_couleur(Couleur::Rouge);
    f.canva().ligne({0.0f,0.5f}, {0.5f,0.0f});

    f.afficher();
  }


  {
    ArrayXf x = linspace(0,1,100);
    Figures figs;
    Figure f = figs.subplot();
    f.plot(x, x);
    f.titres("Titre principal", "Axe x", "Axe y");
    figs.subplot().plot(x, -x);
    figs.afficher();
  }



  // Vérifie pas de planté si vecteur vide
  {
    msg("Test pas de planté si vecteurs vides...");

    Figures f;
    ArrayXf x;
    f.subplot(2,1,1).plot(x);
    f.subplot(2,1,1).titre("Vecteur vide");
    f.subplot(2,1,2).plot_psd(x);
    f.subplot(2,1,2).titre("Vecteur vide (PSD)");
    f.afficher();
    msg("ok.");
  }

  {
    msg("Test signal constant...");
    ArrayXf x = Eigen::ArrayXf::Constant(100, -2.24e38f);
    x += 1e30f * randn(100);
    Figure f;
    f.plot(x);
    f.titre("Signal constant");
    f.afficher();
  }

  {
    msg("Test échelle logarithmique...");
    Figures figs;
    int n = 30;
    ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));


    x(n-1) = 0;
    auto f = figs.subplot();
    f.plot(x, "-ob");
    f.titre("Echelle linéaire");
    f = figs.subplot();
    f.axes().def_echelle("lin", "log");
    f.plot(x, "-ob");
    f.titre("Echelle logarithmique");
    figs.afficher();
  }

  {
    msg_majeur("Test échel log 2");
    ArrayXf x(4), y(4);
    x << 2.0206, 3.0206, 4.0206, 5.0206;
    //y << 0.0371617, 0.0226223, 6.84598e-37, 0;
    y << 0.5371617, 0.0226223, 0.01, 0.005;
    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 2");
    f.afficher("Test log 2");
  }

  {
    msg_majeur("Test échel log 3");
    ArrayXf x(4), y(4);
    x << 2.0206, 3.0206, 4.0206, 5.0206;
    y << 0.5371617, 0.0226223, 0.01, 0.0;
    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 3");
    f.afficher("Test log 3");
  }

  {
    msg_majeur("Test échel log 4");
    ArrayXf x(4), y(4);
    x << 2.0206, 3.0206, 4.0206, 5.0206;
    y << 0.5371617, 0.0226223, 6.84598e-37, 0;

    msg("y = {}", y);

    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 4");
    f.afficher("Test log 4");
  }


  {
    msg("Test plot min / max...");
    Figures figs;
    int n = 30;
    ArrayXf t = linspace(0, n-1, n);
    ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));
    x.tail(2).setZero();
    auto f = figs.subplot(211);
    auto c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
    c.def_couleur({128,128,255});
    f.plot(t, x, "k-o", "x");

    c = f.plot_minmax(t, x.sqrt() * 0.7, x.sqrt() * 1.3);
    c.def_couleur({255,128,128});
    f.plot(t, x.sqrt(), "k-o", "y");
    f.titre("Echelle linéaire");



    f = figs.subplot(212);
    f.axes().def_echelle("lin", "log");
    c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
    c.def_couleur({128,128,255});
    c = f.plot_minmax(t, x.sqrt() * 0.7, x.sqrt() * 1.3);
    c.def_couleur({255,128,128});
    f.plot(t, x, "k-o");
    f.plot(t, x.sqrt(), "k-o");
    f.titre("Echelle logarithmique");
    figs.afficher();
  }

  {
      msg("Test plot min / max (2)...");
      Figure f;
      int n = 30;
      ArrayXf t = linspace(0, n-1, n);
      ArrayXf x = (ArrayXf::Ones(n)*10).pow(linspace(0,-4,n));
      auto c = f.plot(t, x, "b-o", "x");
      c.def_σ(x * 0.1);
      c = f.plot(t, x + 0.4, "g-o", "y");
      c.def_σ(x * 0.2);
      f.titre("Plots min / max (linéaire)");
      f.afficher();
    }

  {
    msg("test BER...");
    int n = 10;
    ArrayXf ber(n), σ(n);
    ber   << 0.452283,0.323602,0.189487,0.0862274,0.0218109,0.00124748,1.00604e-05,0,0,0;
    σ << 0.0388947,0.0587929,0.00844648,0.00233075,0.0018141,0.000457055,3.01811e-05,0,0,0;
    Figure f;
    ArrayXf t = linspace(0, n-1, n);
    auto c = f.plot(t, ber, "b-o", "x");
    f.plot(t, ber/2, "m-s", "x/2");
    c.def_σ(σ);
    f.titre("BER");
    f.axes().def_echelle("lin", "log");
    f.afficher();
  }

  {
    msg("Test accu...");

  }

  {
    msg("Test cercle...");
  }



  return 0;
}
