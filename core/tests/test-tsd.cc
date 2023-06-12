#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"

using namespace std;


void test_tsd()
{

  {
    msg("Test vec...");

    soit x = linspace(0, 4, 5);
    soit s = x.somme(), m = x.moyenne();

    msg("s = {}, m = {}", s, m);
    assertion((s == 10) && (m == 2));
  }

  {
    msg("Test linspace...");
    soit x = linspace(1, 9, 9);

    msg("x :       {}", x);
    msg("diff(x) : {}", diff(x));
    msg("ones    : {}", Vecf::ones(8));
    msg("err     : {}", diff(x) - Vecf::ones(8));

    assertion(x.rows() == 9);
    assertion(x(0) == 1);
    assertion(diff(x).est_approx(Vecf::ones(8)));
  }

  {
    msg("Test irange...");
    soit x = intervalle_entier(1, 9);
    assertion(
           (x.rows() == 9)
        && (x(0) == 1)
        && (diff(x).est_approx(Veci::ones(8))));
  }

  {
    msg("Test trange...");
    soit x = intervalle_temporel(10, 10);
    msg("x       = {}", x);
    msg("diff(x) = {}", diff(x));
    msg("cst     = {}", Vecf::constant(9, 0.1f));
    msg("err     = {}", Vecf::constant(9, 0.1f) - diff(x));
    assertion(x.rows() == 10);
    assertion(x(0) == 0);
    assertion(diff(x).est_approx(Vecf::constant(9, 0.1f)));
  }

  {
    msg("Test logspace...");
    soit x = logspace(-1, 2, 4);

    msg("x: {}", x);

    soit ε = 1e-8;
    assertion(x.rows() == 4);
    assertion(abs(x(0) - 0.1) < ε);
    assertion(abs(x(1) - 1) < ε);
    assertion(abs(x(2) - 10) < ε);
    assertion(abs(x(3) - 100) < ε);
  }

  {
    soit x = Vecf::ones(7);
    assertion(length(x) == 7);
  }


  {
    msg("Test vconcat...");
    soit a = Vecf::random(5),
         b = Vecf::random(7);
    soit c = a | b;

    msg("a = {}", a);
    msg("b = {}", b);
    msg("c = {}", c);

    assertion(c.rows() == 12);
    assertion(a.est_approx(c.head(5)));
    assertion(b.est_approx(c.tail(7)));
  }

  {
    msg("Test modulo...");
    soit check_modulo = [&](double a, double m, double b, double ε2 = 1e-13)
    {
      soit err = abs(modulo(a, m) - b);
      si(err >= ε2)
        échec("a = {}, b = {}, m = {}, modulo = {}, err = {}", a, b, m, modulo(a, m), err);
    };

    check_modulo(0, 1, 0);
    check_modulo(1.5, 1, 0.5);
    check_modulo(-1000.7, 1, 0.3);
  }

  {
    msg("Test modulo [0,2π[...");

    soit ε = 1e-10;

    soit check_wrap = [&](double a, double b, double ε2 = 1e-16)
    {
      soit err = abs(modulo_2π(a) - b);
      si(err >= ε2)
        échec("a = {}, b = {}, modulo = {}, err = {}", a, b, modulo_2π(a), err);
    };

    check_wrap(0, 0);
    check_wrap(π, π);
    check_wrap(2*π, 0);
    check_wrap(2*π-ε, 2*π-ε);
    check_wrap(2*π+ε, ε);
    check_wrap(5*π, π);
    check_wrap(-ε, 2*π-ε);
    check_wrap(-π, π);
    check_wrap(-2*π+ε, ε);
    check_wrap(-2*π-ε, 2*π-ε);
    check_wrap(1000*π+ε, ε, 1e-12);
  }

  {
    msg("Test modulo [-π,π[...");

    soit ε = 1e-10;

    soit check_wrap = [&](double a, double b, double ε2 = 1e-16)
    {
      soit err = abs(modulo_pm_π(a) - b);
      assertion(err < ε2);
    };

    check_wrap(0, 0);
    check_wrap(π, -π);
    check_wrap(2*π, 0);
    check_wrap(2*π-ε, -ε);
    check_wrap(2*π+ε, ε);
    check_wrap(5*π, -π);
    check_wrap(-ε, -ε);
    check_wrap(-π, -π);
    check_wrap(-2*π+ε, ε);
    check_wrap(-2*π-ε, -ε);
  }

  {
    msg("Test modulo...");
    soit n = 2000;
    soit x = linspace(-5,5,n);
    x = square(x);
    //TODO x = x.unaryExpr([](const float x) {retourne modulo_2π(x);}).eval();

    pour(auto i = 0; i < n; i++)
      x(i) = modulo_2π(x(i));

    soit y = déplie_phase(x);
    // Vérifie pas de saut de phase
    assertion(abs(diff(y)).valeur_max() < 0.1);
    // Vérifie valeur ok modulo 2pi
    pour(auto i = 0; i < n; i++)
    {
      soit err = abs(modulo_pm_π(y(i) - x(i)));
      //msg("z({}) = {}, wrap={}, x = {}, erreur = {}", i, y(i), wrap_pm_pi(y(i) - x(i)), x(i), err);
      assertion(err < 1e-5);
    }
  }

  {
    msg("Test 'trouve'...");
    soit x = Vecb::valeurs({0, 0, 1, 0, 1});
    soit idx = trouve(x);
    assertion((idx.size() == 2) && (idx[0] == 2) && (idx[1] == 4));
  }

  {
    msg("Test 'trouve premier'...");
    soit x = Vecb::valeurs({0, 0, 1, 0, 1});
    soit idx = trouve_premier(x);
    assertion(idx == 2);
  }

  {
    msg("Test pad zéros...");
    soit x1 = linspace(0, 4, 5),
         x2 = linspace(0, 3, 4);
    soit [x1p, x2p] = pad_zeros(x1, x2);

    assertion(x1p.est_approx(x1));
    Vecf x2p_ref(5);
    x2p_ref.head(4) = x2;
    x2p_ref(4) = 0;
    assertion(x2p.est_approx(x2p_ref));

    soit [x1b, x2b] = pad_zeros(x1, x2, oui);

    assertion(x1b.rows() == 8);
    assertion(x2b.rows() == 8);
    assertion(x1b.head(5).est_approx(x1));
    assertion(x2b.head(4).est_approx(x2));
  }

  {
    msg("Test rotation vec...");
    soit x = linspace(0, 5, 6);
    soit y = rotation_vec(x, 2);
    soit yref = Vecf::valeurs({2, 3, 4, 5, 0, 1});
    assertion(y.est_approx(yref));
  }

  {
    msg("Test diff...");
    soit x = randn(10);
    soit y = diff(x);
    soit yref = x.tail(9) - x.head(9);
    assertion(y.est_approx(yref));
  }

  {
    msg("Test cumsum...");
    soit x     = linspace(0, 99, 100);
    soit y     = cumsum(x);
    soit yref  = x * (x + 1) / 2;
    assertion(y.est_approx(yref));
  }


  /*{
    msg("Test subarray1d...");
    soit x  = linspace(0, 9, 10);
    soit s1 = subarray1d(x, 0, 5, 2);
    soit s2 = subarray1d(x, 1, 5, 2);

    soit pairs = linspace(0, 8, 5); // 0,2,4,6,8
    assertion(pairs.est_approx(s1, 1e-8f));
    soit impairs = linspace(1, 9, 5); // 1,3,5,7,9
    assertion(impairs.est_approx(s2, 1e-8f));

    s1.setZero();
    assertion(x(2) == 0);
  }*/

  {
    msg("Test pp2...");
    auto pp2 = prochaine_puissance_de_2;
    assertion(pp2(1) == 1);
    assertion(pp2(2) == 2);
    assertion(pp2(3) == 4);
    assertion(pp2(4) == 4);
    assertion(pp2(5) == 8);
    assertion(pp2(1 << 16) == (1 << 16));
    assertion(pp2((1 << 16) - 1) == (1 << 16));
  }

  {
    msg("Test sousech...");
    soit x      = linspace(0, 9, 10);
    soit xe     = sousech(x, 2),
         xe_ref = linspace(0, 8, 5); // 0,2,4,6,8
    assertion(xe_ref.est_approx(xe));
  }

  {
    msg("Test surech...");
    soit x    = linspace(0, 4, 5);
    soit y    = surech(x, 2),
         yref = Vecf::valeurs({0, 0, 1, 0, 2, 0, 3, 0, 4, 0});

    assertion(yref.est_approx(y));
  }

  {
    // db = 10 * log10(p)
    // p = 2  => 10 * log10(2) = 3 dB
    // p = 100 => 20 dB

    msg("pow2db(2) = {}", pow2db(2));
    msg("pow2db(100) = {}", pow2db(100));

    msg("db2pow(0) = {}",  db2pow(0));
    msg("db2pow(3) = {}",  db2pow(3));
    msg("db2pow(-3) = {}", db2pow(-3));

    assertion(abs(pow2db(2)   - 3)    < 0.02);
    assertion(abs(pow2db(100) - 20)   < 1e-6);
    assertion(abs(db2pow(0)   - 1)    < 1e-6);
    assertion(abs(db2pow(3)   - 2)    < 0.1);
    assertion(abs(db2pow(-3)   - 0.5)    < 0.1);
    assertion(abs(db2pow(20)  - 100)  < 1e-6);
    assertion(abs(db2pow(-20) - 0.01) < 1e-6);
  }


  {
    assertion(abs(rad2deg(π/4) - 45) < 1e-15);
    assertion(abs(rad2deg(-π/4) + 45) < 1e-15);
    msg("d2r(45) = {}, π/4={}, err={}", deg2rad(45), π/4, deg2rad(45) - π/4);
    assertion(abs(deg2rad(45) - π/4) < 1e-15);
    assertion(abs(deg2rad(-45) + π/4) < 1e-15);
  }

  {
    msg("Test randn...");
    soit n = 10000;
    soit x = randn(n);
    assertion(x.rows() == n);
    assertion(abs(x.somme()) < 4 * sqrt(n * 1.0));

    soit var = square(x).moyenne();
    msg("Variance mesurée : {}", var);
    assertion(abs(var - 1) < 1e-2);
    msg("Test randu...");
    msg("randu = {}", randu(10, 0, 1));
    x = randu(n);
    assertion(x.rows() == n);
    assertion(x.valeur_max() <= 1);
    assertion_msg(x.valeur_min() >= -1, "randu min = {}", x.valeur_min());
    assertion(abs(x.moyenne() - 0.0) < 1e-2);
    var = square(x).moyenne() - x.moyenne() * x.moyenne();
    msg("Variance mesurée : {} (théorique : {})", var, 4 * 1.0f/12);
    assertion(abs(var - 4 * 1.0f/12) < 1e-2);


    msg("Test randi...");
    n = 10000;
    soit xi = randi(10, n);

    //msg("xi = {}", xi);
    //msg("xi.as float = {}", xi.as<float>());


    assertion(xi.rows() == n);
    assertion(xi.valeur_min() == 0);
    assertion(xi.valeur_max() == 9);

    soit moy = xi.as<float>().moyenne();
    soit xif = xi.as<float>();
    soit moy2 = xif.moyenne();

    msg("moy = {}, moy2 = {}", moy, moy2);
    msg("moyd = {}", xi.as<double>().moyenne());

    /*{
      soit r  = randi(10, 10);
      msg("randi(10,10): {}", r);
      msg("randi(10,10).as float: {}", r.as<float>());
      msg("randi(10,10).as float moy: {}", r.moyenne());
    }*/

    assertion(abs(moy - 4.5) < 1e-1);
    {
      msg("Test randb");
      n = 10;
      x = randb(n).as<float>();
      //msg("randb = {} ({}) ({})", x, randb(n), randb(n).as<float>());
      assertion(x.rows() == n);
      assertion(x.valeur_max() == 1);
      assertion(x.valeur_min() == 0);
      assertion((x * (x - 1)).est_approx(0));
      assertion(abs(x.moyenne() - 0.5) < 1e-2);
    }
  }

  {
    msg("Test sigexp...");
    soit n = 1000;
    soit x    = sigexp(0.1, n);
    soit xref = tsd::polar(2 * π_f * 0.1f * linspace(0, n-1, n));

    si(tests_debug_actif)
    {
      Figure f;
      f.plot(real(x), "b-", "x");
      f.plot(real(xref), "g-", "xref");
      f.afficher();
      f.clear();
      f.plot(real(x) - real(xref), "r-", "erreur");
      f.afficher();
    }

    soit err = sqrt(abs2(x - xref).moyenne());
    soit mag = sqrt(min(abs2(xref).moyenne(), abs2(x).moyenne()));
    msg("Erreur RMS = {}, mag = {}, err relative = {}", err, mag, err / mag);

    assertion(x.est_approx(xref, 0.5e-4));
  }

  {
    msg("Test signyquist...");
    soit n = 100;
    soit x = signyquist(n);
    pour(auto i = 0; i < n; i += 2)
    {
      assertion(x(i) == -1);
      assertion(x(i+1) == 1);
    }
    x = signyquist(n+1);
    pour(auto i = 0; i < n; i += 2)
    {
      assertion(x(i) == -1);
      assertion(x(i+1) == 1);
    }
    assertion_msg(x(n) == -1, "x = {}", x);
  }

  {
    msg("Test sigimp...");
    soit n = 100;
    soit x = sigimp(n, 5);
    assertion(x.head(5).est_approx(0));
    assertion(x.tail(n-6).est_approx(0));
    assertion(x(5) == 1);
  }

  {
    msg("Test sigsin...");
    soit n    = 1000;
    soit x    = sigsin(0.1, n);
    soit xref = sin(2 * π_f * 0.1 * linspace(0, n-1, n));
    assertion(x.est_approx(xref, 0.5e-4));
  }

  {
    msg("Test sigcos...");
    soit n    = 1000;
    soit x    = sigcos(0.1, n);
    soit xref = cos(2 * π_f * 0.1 * linspace(0, n-1, n));
    assertion(x.est_approx(xref, 0.5e-4));
  }

  {
    msg("Test sigtri...");
    soit x = sigtri(10, 20);
    assertion(x.valeur_max() == 1);
    assertion(x.valeur_min() == -1);
    assertion(abs(x.moyenne()) < 1e-7);
    assertion(x.head(6).est_approx(linspace(-1,1,6)));
    assertion(x.segment(5,6).est_approx(linspace(1,-1,6)));
    assertion(x.segment(10,6).est_approx(linspace(-1,1,6)));
    assertion(x.tail(5).est_approx(linspace(1,-1+2.0f/5,5)));
  }

  {
    msg("Test sigcar...");
    soit x = sigcar(10, 20);
    assertion(x.head(5).est_approx(-Vecf::ones(5)));
    assertion(x.segment(5,5).est_approx(Vecf::ones(5)));
    assertion(x.segment(10,5).est_approx(-Vecf::ones(5)));
    assertion(x.tail(5).est_approx(Vecf::ones(5)));
  }

  {
    msg("Test source_ohc...");
    soit src = source_ohc(0.1);
    soit n   = 1000;
    Veccf x(n);
    x.head(n/2) = src->step(n/2);
    x.tail(n/2) = src->step(n/2);
    soit xref = tsd::polar(2 * π_f * 0.1 * linspace(0, n-1, n));
    assertion(x.est_approx(xref, 0.5e-4));
  }

  {
    msg("Test source_ohr...");
    soit src = source_ohr(0.1);
    soit n   = 1000;
    Vecf x(n);
    x.head(n/2) = src->step(n/2);
    x.tail(n/2) = src->step(n/2);
    Vecf xref = cos(2 * π_f * 0.1 * linspace(0, n-1, n));
    assertion(x.est_approx(xref, 0.5e-4));
  }

  {
    msg("Test tampon de données...");
    soit cnt = 0,
         n   = 16 * 512;
    soit X   = randn(n);
    soit t   = tampon_création<float>(512,
        // Définition de la callback
        [&](const Vecf &x)
        {
          assertion(x.rows() == 512);
          assertion(x.est_approx(X.segment(cnt, 512)));
          cnt += 512;
        });

    t->step(X.head(100));
    t->step(X.segment(100, 1000));
    t->step(X.tail(n - 100 - 1000));
    assertion(cnt == n);
  }
}
