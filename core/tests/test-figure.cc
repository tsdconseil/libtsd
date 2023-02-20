#include "tsd/tsd-all.hpp"
#include "Eigen/Core"

static void test_unites_unit(const std::vector<double> &tics, const std::string &unit = "")
{
  std::string s1 = "";
  pour(auto t: tics)
    s1 += fmt::format("{},", t);



  soit [expo, ndigits] = tsd::vue::unites::calc_expo_nb_chiffres_commun(tics, unit);

  msg("test_unites : vec = [{}], unit={}, expo={}, ndigits={}", s1, unit, expo, ndigits);

  pour(auto t : tics)
  {
    soit s = tsd::vue::unites::valeur_vers_chaine(t, unit, expo, ndigits);
    msg(" label : {}", s);
  }
}

entier test_unites()
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

  retourne 0;
}

entier test_figure()
{
  // Test affichage long
  {
    Figure f;
    soit n = 15000;
    soit x = 0.1 * randn(n);
    pour(auto i = 1000; i < 15000; i += 1000)
      x(i) = 1;
    f.plot(x);
    f.afficher("Test affichage 15000 points");
  }



  // Test échelle
  {
    Figure f;

    soit vt = linspace(0, 179.0/60, 180);
    f.plot(vt, vt, "", "Temps (minutes)");
    f.afficher("Test échelle");
  }



  {
    Figure f;

    soit x = linspace(0,1,100);
    //ArrayXf y = ArrayXf::Ones(10);

    f.plot(x, x, "", "valeurs de x");
    f.plot(x, square(x), "", "valeurs de x^2");
    //f.enregistrer("./test-log/ones.png");

    f.titres("Titre principal", "Axe x", "Axe y");

    f.canva().texte(0.1, 0.2, "Texte @ 0,1x0,2", 0.2, 0.1);
    f.canva().set_couleur(Couleur::Rouge);
    f.canva().ligne({0.0f,0.5f}, {0.5f,0.0f});

    f.afficher();
  }


  {
    soit x = linspace(0,1,100);
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
    Vecf x;
    f.subplot(2,1,1).plot(x);
    f.subplot(2,1,1).titre("Vecteur vide");
    f.subplot(2,1,2).plot_psd(x);
    f.subplot(2,1,2).titre("Vecteur vide (PSD)");
    f.afficher();
    msg("ok.");
  }

  {
    msg("Test signal constant...");
    soit x = Vecf::constant(100, -2.24e38f);
    x += 1e30f * randn(100);
    Figure f;
    f.plot(x);
    f.titre("Signal constant");
    f.afficher();
  }

  {
    msg("Test échelle logarithmique...");
    Figures figs;
    soit n = 30;
    soit x = pow(Vecf::ones(n)*10,linspace(0,-4,n));


    x(n-1) = 0;
    soit f = figs.subplot();
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
    ///*Eigen::ArrayXf*/Vecf x(4), y(4);
    //x << 2.0206, 3.0206, 4.0206, 5.0206;
    //y << 0.5371617, 0.0226223, 0.01, 0.005;

    soit x = Vecf::valeurs({2.0206, 3.0206, 4.0206, 5.0206});
    soit y = Vecf::valeurs({0.5371617, 0.0226223, 0.01, 0.005});

    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 2");
    f.afficher("Test log 2");
  }

  {
    msg_majeur("Test échel log 3");
    soit x = Vecf::valeurs({2.0206, 3.0206, 4.0206, 5.0206});
    soit y = Vecf::valeurs({0.5371617, 0.0226223, 0.01, 0.0});
    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 3");
    f.afficher("Test log 3");
  }

  {
    msg_majeur("Test échel log 4");
    soit x = Vecf::valeurs({2.0206, 3.0206, 4.0206, 5.0206});
    soit y = Vecf::valeurs({0.5371617, 0.0226223, 6.84598e-37, 0});

    msg("y = {}", y);

    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(x, y, "-ob", "Test log 4");
    f.afficher("Test log 4");
  }


  {
    msg("Test plot min / max...");
    Figures figs;
    soit n = 30;
    soit t = linspace(0, n-1, n);
    soit x = pow(Vecf::ones(n)*10,linspace(0,-4,n));
    x.tail(2).setConstant(0);
    soit f = figs.subplot(211);
    soit c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
    c.def_couleur({128,128,255});
    f.plot(t, x, "k-o", "x");

    c = f.plot_minmax(t, sqrt(x) * 0.7, sqrt(x) * 1.3);
    c.def_couleur({255,128,128});
    f.plot(t, sqrt(x), "k-o", "y");
    f.titre("Echelle linéaire");



    f = figs.subplot(212);
    f.axes().def_echelle("lin", "log");
    c = f.plot_minmax(t, x - 0.2 * x, x + 0.2 * x);
    c.def_couleur({128,128,255});
    c = f.plot_minmax(t, sqrt(x) * 0.7, sqrt(x) * 1.3);
    c.def_couleur({255,128,128});
    f.plot(t, x, "k-o");
    f.plot(t, sqrt(x), "k-o");
    f.titre("Echelle logarithmique");
    figs.afficher();
  }

  {
      msg("Test plot min / max (2)...");
      Figure f;
      soit n = 30;
      soit t = linspace(0, n-1, n);
      soit x = pow(Vecf::ones(n)*10,linspace(0,-4,n));
      soit c = f.plot(t, x, "b-o", "x");
      c.def_σ(x * 0.1);
      c = f.plot(t, x + 0.4, "g-o", "y");
      c.def_σ(x * 0.2);
      f.titre("Plots min / max (linéaire)");
      f.afficher();
    }



  {
    msg("test BER...");
    soit n = 10;
    soit ber = Vecf::valeurs({0.452283,0.323602,0.189487,0.0862274,0.0218109,0.00124748,1.00604e-05,0,0,0});
    soit σ   = Vecf::valeurs({0.0388947,0.0587929,0.00844648,0.00233075,0.0018141,0.000457055,3.01811e-05,0,0,0});
    Figure f;
    soit t = linspace(0, n-1, n);
    soit c = f.plot(t, ber, "b-o", "x");
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

  pour(auto ny : {10, 100, 1000})
  {
    msg("Test img, ny = {}", ny);
    soit nx = 100;
    Tabf M(nx,ny);
    soit x = linspace(0, 1, nx);
    soit y = linspace(0, 1, ny);
    pour(auto ix = 0; ix < nx; ix++)
      pour(auto iy = 0; iy < ny; iy++)
        M(ix,iy) = sin(x(ix)*10) * cos(y(iy)*10);


    Figure f;
    f.plot_img(M);
    f.afficher();
  }


  retourne 0;
}
