#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"

using namespace std;


int test_tsd()
{

  {
    msg("Test linspace...");
    ArrayXf x = linspace(1, 9, 9);
    tsd_assert(
           (x.rows() == 9)
        && (x(0) == 1)
        && (diff(x).isApprox(ArrayXf::Ones(8))));
  }

  {
    msg("Test irange...");
    ArrayXi x = intervalle_entier(1, 9);
    tsd_assert(
           (x.rows() == 9)
        && (x(0) == 1)
        && (diff(x).isApprox(ArrayXi::Ones(8))));
  }

  {
    msg("Test trange...");
    ArrayXf x = intervalle_temporel(10, 10);
    msg("x = {}", x);
    tsd_assert((x.rows() == 10)
        && (x(0) == 0)
        && (diff(x).isApprox(ArrayXf::Constant(9, 0.1f))));
  }

  {
    msg("Test logspace...");
    ArrayXf x = logspace(-1, 2, 4);

    msg("x: {}", x);

    auto ε = 1e-8;
    tsd_assert(x.rows() == 4);
    tsd_assert(abs(x(0) - 0.1) < ε);
    tsd_assert(abs(x(1) - 1) < ε);
    tsd_assert(abs(x(2) - 10) < ε);
    tsd_assert(abs(x(3) - 100) < ε);
  }

  {
    ArrayXf x = ArrayXf::Ones(7);
    tsd_assert(length(x) == 7);
  }


  {
    msg("Test vconcat...");
    ArrayXf a = ArrayXf::Random(5);
    ArrayXf b = ArrayXf::Random(7);
    ArrayXf c = a | b;
    tsd_assert(c.rows() == 12);
    tsd_assert(a.isApprox(c.head(5)));
    tsd_assert(b.isApprox(c.tail(7)));
  }

  {
    msg("Test modulo...");
    auto check_modulo = [&](double a, double m, double b, double ε2 = 1e-13)
    {
      auto err = abs(modulo(a, m) - b);
      if(err >= ε2)
        echec("a = {}, b = {}, m = {}, modulo = {}, err = {}", a, b, m, modulo(a, m), err);
    };

    check_modulo(0, 1, 0);
    check_modulo(1.5, 1, 0.5);
    check_modulo(-1000.7, 1, 0.3);
  }

  {
    msg("Test modulo [0,2π[...");

    auto ε = 1e-10;

    auto check_wrap = [&](double a, double b, double ε2 = 1e-16)
    {
      auto err = abs(modulo_2π(a) - b);
      if(err >= ε2)
        echec("a = {}, b = {}, modulo = {}, err = {}", a, b, modulo_2π(a), err);
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

    auto ε = 1e-10;

    auto check_wrap = [&](double a, double b, double ε2 = 1e-16)
    {
      auto err = abs(modulo_pm_π(a) - b);
      tsd_assert(err < ε2);
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
    int n = 2000;
    ArrayXf x = linspace(-5,5,2000);
    x = x.square().eval();
    x = x.unaryExpr([](const float x) {return modulo_2π(x);}).eval();
    ArrayXf y = déplie_phase(x);
    // Vérifie pas de saut de phase
    tsd_assert(diff(y).abs().maxCoeff() < 0.1);
    // Vérifie valeur ok modulo 2pi
    for(auto i = 0; i < n; i++)
    {
      auto err = abs(modulo_pm_π(y(i) - x(i)));
      //msg("z({}) = {}, wrap={}, x = {}, erreur = {}", i, y(i), wrap_pm_pi(y(i) - x(i)), x(i), err);
      tsd_assert(err < 1e-5);
    }
  }

  {
    msg("Test 'trouve'...");
    ArrayXb x(5);
    x << 0, 0, 1, 0, 1;
    auto idx = trouve(x);
    tsd_assert((idx.size() == 2) && (idx[0] == 2) && (idx[1] == 4));
  }

  {
    msg("Test 'trouve premier'...");
    ArrayXb x(5);
    x << 0, 0, 1, 0, 1;
    auto idx = trouve_premier(x);
    tsd_assert(idx == 2);
  }

  {
    msg("Test pad zéros...");
    ArrayXf x1 = linspace(0, 4, 5),
            x2 = linspace(0, 3, 4);
    auto [x1p, x2p] = pad_zeros(x1, x2);

    tsd_assert(x1p.isApprox(x1));
    ArrayXf x2p_ref(5);
    x2p_ref.head(4) = x2;
    x2p_ref(4) = 0;
    tsd_assert(x2p.isApprox(x2p_ref));

    auto [x1b, x2b] = pad_zeros(x1, x2, true);

    tsd_assert(x1b.rows() == 8);
    tsd_assert(x2b.rows() == 8);
    tsd_assert(x1b.head(5).isApprox(x1));
    tsd_assert(x2b.head(4).isApprox(x2));
  }

  {
    msg("Test rotation vec...");
    ArrayXf x = linspace(0, 5, 6);
    ArrayXf y = rotation_vec(x, 2);
    ArrayXf yref(6);
    yref << 2, 3, 4, 5, 0, 1;
    tsd_assert(y.isApprox(yref));
  }

  {
    msg("Test diff...");
    ArrayXf x = randn(10);
    ArrayXf y = diff(x);
    ArrayXf yref = x.tail(9) - x.head(9);
    tsd_assert(y.isApprox(yref));
  }

  {
    msg("Test cumsum...");
    ArrayXf x     = linspace(0, 99, 100);
    ArrayXf y     = cumsum(x);
    ArrayXf yref  = x * (x + 1) / 2;
    tsd_assert(y.isApprox(yref));
  }


  {
    msg("Test subarray1d...");
    ArrayXf x = linspace(0, 9, 10);
    auto s1 = subarray1d(x, 0, 5, 2);
    auto s2 = subarray1d(x, 1, 5, 2);

    ArrayXf pairs = linspace(0, 8, 5); // 0,2,4,6,8
    tsd_assert(pairs.isApprox(s1, 1e-8f));
    ArrayXf impairs = linspace(1, 9, 5); // 1,3,5,7,9
    tsd_assert(impairs.isApprox(s2, 1e-8f));

    s1.setZero();
    tsd_assert(x(2) == 0);
  }

  {
    msg("Test pp2...");
    tsd_assert(prochaine_puissance_de_2(1) == 1);
    tsd_assert(prochaine_puissance_de_2(2) == 2);
    tsd_assert(prochaine_puissance_de_2(3) == 4);
    tsd_assert(prochaine_puissance_de_2(4) == 4);
    tsd_assert(prochaine_puissance_de_2(5) == 8);
    tsd_assert(prochaine_puissance_de_2(1 << 16) == (1 << 16));
    tsd_assert(prochaine_puissance_de_2((1 << 16) - 1) == (1 << 16));
  }

  {
    msg("Test sousech...");
    ArrayXf x  = linspace(0, 9, 10);
    ArrayXf xe = sousech(x, 2);

    ArrayXf xe_ref = linspace(0, 8, 5); // 0,2,4,6,8
    tsd_assert(xe_ref.isApprox(xe));
  }

  {
    msg("Test surech...");
    ArrayXf x = linspace(0, 4, 5);
    ArrayXf y = surech(x, 2);

    ArrayXf yref(10);
    yref << 0, 0, 1, 0, 2, 0, 3, 0, 4, 0;

    tsd_assert(yref.isApprox(y));
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

    tsd_assert(abs(pow2db(2)   - 3)    < 0.02);
    tsd_assert(abs(pow2db(100) - 20)   < 1e-6);
    tsd_assert(abs(db2pow(0)   - 1)    < 1e-6);
    tsd_assert(abs(db2pow(3)   - 2)    < 0.1);
    tsd_assert(abs(db2pow(-3)   - 0.5)    < 0.1);
    tsd_assert(abs(db2pow(20)  - 100)  < 1e-6);
    tsd_assert(abs(db2pow(-20) - 0.01) < 1e-6);
  }


  {
    tsd_assert(abs(rad2deg(π/4) - 45) < 1e-15);
    tsd_assert(abs(rad2deg(-π/4) + 45) < 1e-15);
    msg("d2r(45) = {}, π/4={}, err={}", deg2rad(45), π/4, deg2rad(45) - π/4);
    tsd_assert(abs(deg2rad(45) - π/4) < 1e-15);
    tsd_assert(abs(deg2rad(-45) + π/4) < 1e-15);
  }

  {
    msg("Test randn...");
    int n = 10000;
    ArrayXf x = randn(n);
    tsd_assert(x.rows() == n);
    tsd_assert(abs(x.sum()) < 4 * sqrt(n * 1.0));

    auto var = x.square().mean();
    msg("Variance mesurée : {}", var);
    tsd_assert(abs(var - 1) < 1e-2);
    msg("Test randu...");
    msg("randu = {}", randu(10));
    x = randu(n);
    tsd_assert(x.rows() == n);
    tsd_assert(x.maxCoeff() <= 1);
    tsd_assert(x.minCoeff() >= 0);
    tsd_assert(abs(x.mean() - 0.5) < 1e-2);
    var = x.square().mean() - x.mean() * x.mean();
    msg("Variance mesurée : {} (théorique : {})", var, 1.0f/12);
    tsd_assert(abs(var - 1.0f/12) < 1e-2);
    msg("Test randi...");
    ArrayXi xi = randi(10, n);
    tsd_assert(xi.rows() == n);
    tsd_assert(xi.minCoeff() == 0);
    tsd_assert(xi.maxCoeff() == 9);
    auto moy = xi.cast<float>().mean();
    msg("moy = {}", moy);
    tsd_assert(abs(moy - 4.5) < 1e-1);
    {
      msg("Test randb");
      x = randb(n).cast<float>();
      tsd_assert(x.rows() == n);
      tsd_assert(x.maxCoeff() == 1);
      tsd_assert(x.minCoeff() == 0);
      tsd_assert((x * (x - 1)).isApproxToConstant(0));
      tsd_assert(abs(x.mean() - 0.5) < 1e-2);
    }
  }

  {
    msg("Test sigexp...");
    int n = 1000;
    ArrayXcf x    = sigexp(0.1, n);
    ArrayXcf xref = tsd::polar(2 * π_f * 0.1f * linspace(0, n-1, n));

    if(tests_debug_actif)
    {
      Figure f;
      f.plot(x.real(), "b-", "x");
      f.plot(xref.real(), "g-", "xref");
      f.afficher();
      f.clear();
      f.plot(x.real() - xref.real(), "r-", "erreur");
      f.afficher();
    }

    auto err = sqrt((x - xref).abs2().mean());
    auto mag = sqrt(min(xref.abs2().mean(), x.abs2().mean()));
    msg("Erreur RMS = {}, mag = {}, err relative = {}", err, mag, err / mag);

    tsd_assert(x.isApprox(xref, 0.5e-4));
  }

  {
    msg("Test sigimp...");
    int n = 100;
    ArrayXf x = sigimp(n, 5);
    tsd_assert(x.head(5).isApproxToConstant(0));
    tsd_assert(x.tail(n-6).isApproxToConstant(0));
    tsd_assert(x(5) == 1);
  }

  {
    msg("Test sigsin...");
    int n = 1000;
    ArrayXf x    = sigsin(0.1, n);
    ArrayXf xref = (2 * π_f * 0.1 * linspace(0, n-1, n)).sin();
    tsd_assert(x.isApprox(xref, 0.5e-4));
  }

  {
    msg("Test sigcos...");
    int n = 1000;
    ArrayXf x    = sigcos(0.1, n);
    ArrayXf xref = (2 * π_f * 0.1 * linspace(0, n-1, n)).cos();
    tsd_assert(x.isApprox(xref, 0.5e-4));
  }

  {
    msg("Test sigtri...");
    ArrayXf x = sigtri(10, 20);
    tsd_assert(x.maxCoeff() == 1);
    tsd_assert(x.minCoeff() == -1);
    tsd_assert(abs(x.mean()) < 1e-7);
    tsd_assert(x.head(6).isApprox(linspace(-1,1,6)));
    tsd_assert(x.segment(5,6).isApprox(linspace(1,-1,6)));
    tsd_assert(x.segment(10,6).isApprox(linspace(-1,1,6)));
    tsd_assert(x.tail(5).isApprox(linspace(1,-1+2.0f/5,5)));
  }

  {
    msg("Test sigcar...");
    ArrayXf x = sigcar(10, 20);
    tsd_assert(x.head(5).isApprox(-ArrayXf::Ones(5)));
    tsd_assert(x.segment(5,5).isApprox(ArrayXf::Ones(5)));
    tsd_assert(x.segment(10,5).isApprox(-ArrayXf::Ones(5)));
    tsd_assert(x.tail(5).isApprox(ArrayXf::Ones(5)));
  }

  {
    msg("Test source_ohc...");
    auto src = source_ohc(0.1);
    int n = 1000;
    ArrayXcf x(n);
    x.head(n/2) = src->step(n/2);
    x.tail(n/2) = src->step(n/2);
    ArrayXcf xref = tsd::polar(2 * π_f * 0.1 * linspace(0, n-1, n));
    tsd_assert(x.isApprox(xref, 0.5e-4));
  }

  {
    msg("Test source_ohr...");
    auto src = source_ohr(0.1);
    int n = 1000;
    ArrayXf x(n);
    x.head(n/2) = src->step(n/2);
    x.tail(n/2) = src->step(n/2);
    ArrayXf xref = (2 * π_f * 0.1 * linspace(0, n-1, n)).cos();
    tsd_assert(x.isApprox(xref, 0.5e-4));
  }

  {
    msg("Test tampon de données...");
    int cnt = 0;
    int n = 16 * 512;
    ArrayXf X = randn(n);

    auto t = tampon_création<float>(512,
        // Définition de la callback
        [&](const ArrayXf &x)
        {
          tsd_assert(x.rows() == 512);
          tsd_assert(x.isApprox(X.segment(cnt, 512)));
          cnt += 512;
        });

    t->step(X.head(100));
    t->step(X.segment(100, 1000));
    t->step(X.tail(n - 100 - 1000));
    tsd_assert(cnt == n);
  }

  return 0;
}
