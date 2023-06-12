#include "tsd/tsd-all.hpp"

//#include "Eigen/Core"



// Pour vérif move sémantique
static int64_t xptr;
static Vecf bidon()
{
  soit x = Vecf::ones(5);
  xptr = (int64_t) x.data();
  retourne x;
}

void test_tab()
{
  //soit t = Tab::ones(ℝ, 32, 3,3);
  soit t = Tabf::ones(3,3);

  msg("t = {}", t);

  // Vérificatioon deep copy
  soit t2 = t.clone();

  assertion((t2.rows() == 3) && (t2.cols() == 3));

  t2(0,1) = 7;

  msg("t2 = {}", t2);

  pour(auto i = 0; i < 3; i++)
    pour(auto j = 0; j < 3; j++)
      assertion(t(i,j) == 1);

  {
    msg("Vérif affectation colonne d'un tableau...");
    soit M = Tabf::zeros(5,10);
    soit v1 = Vecf::ones(5),
         v2 = 2 * v1;
    M.col(0) = v1;
    M.col(1) = v2;

    msg("M = \n{}", M);
    msg("v1 = {}", v1);
    msg("v2 = {}", v2);

    assertion(M(0,0) == 1);
    assertion(M(0,1) == 2);
  }


  {
    msg("Vérif sémantique move...");
    soit y = bidon();
    soit yptr = (int64_t) y.data();
    msg("xptr = {}, yptr = {}", xptr, yptr);

    soit z = y.clone();
    soit zptr = (int64_t) z.data();
    msg("zptr = {}", zptr);


    assertion(xptr == yptr);
    assertion(zptr != xptr);
    assertion(y.est_approx(Vecf::ones(5)));

    // Maintenant, vérification que le move est désactivé si la destination
    // est une référence.
    soit X  = Vecf::zeros(10),
         Xh = X.head(5);

    assertion(X.data() == Xh.data());

    Xh = y;
    assertion(X.data() == Xh.data());

    msg("Xh = {}", Xh);
    msg("X = {}", X);

    assertion(X.head(5).est_approx(Vecf::ones(5)));
    assertion(X.tail(5).est_approx(Vecf::zeros(5)));
  }

  t *= 2;

  msg("t * 2 = {}", t);

  pour(auto i = 0; i < 3; i++)
    pour(auto j = 0; j < 3; j++)
      assertion(t(i,j) == 2);

  t += 1;

  msg("t * 2 + 1 = {}", t);

  pour(auto i = 0; i < 3; i++)
    pour(auto j = 0; j < 3; j++)
      assertion(t(i,j) == 3);



  {
    soit x  = Vecf::ones(10),
         x2 = cumsum(x);
    msg("cumsum({}) = {}", x, x2);
    assertion(x2.est_approx(linspace(1,10,10)));
  }

  {
    msg("Test head / tail...");
    soit x = Vecf::ones(10);
    x.head(5).setConstant(0);

    msg("x = {}", x);

    pour(auto i = 0; i < 5; i++)
    {
      msg("x({}) = {}", i, x(i));
      msg("x({}) = {}", i+5, x(i+5));
      assertion(x(i) == 0);
      assertion(x(i+5) == 1);
    }

    x.segment(1, 2).setConstant(5);
    assertion(x(0) == 0);
    assertion(x(1) == 5);
    assertion(x(2) == 5);
    assertion(x(3) == 0);



    x = linspace(0, 19, 20);
    x.head(19) = x.tail(19);
    assertion(x.head(19).est_approx(linspace(1,19,19)));
    assertion(x(19) == 19);
  }


  {
    msg("Test redim...");
    Vecf x;
    x.resize(10);
    assertion(x.dim() == 10);
    pour(auto i = 0; i < 10; i++)
      x(i) = i;
    pour(auto i = 0; i < 10; i++)
      assertion(x(i) == i);

  }

  {
    msg("Test cplx");
    soit x = Veccf::zeros(5);
    x.setConstant(cfloat{0,1});
    soit y = x.conjugate();

    msg("x = {}", x);
    msg("y = {}", y);
    msg("abs(x) = {}", abs(x));
    msg("abs(y) = {}", abs(y));
    msg("abs(y)-1 = {}", abs(y) - 1);
    assertion(abs(y).est_approx(Vecf::ones(5)));
    assertion((x + y).est_approx(Veccf::zeros(5)));
  }

  {
    msg("Test reverse");
    soit x = linspace(0, 9, 10);
    soit y = x.reverse();
    msg("rlin = {}", y);
    assertion((x+y).est_approx(9*Vecf::ones(10)));
  }

  {
    msg("Test segment...");
    soit x = Veccf::zeros(10);

    x.segment(2, 4) = Vecf::ones(4);


    msg("x = {}", x);
    assertion((x(0) == cfloat{0,0}) && (x(2) == cfloat{1,0}));

    x.tail(2) += x.segment(2, 2);

    msg("x2 = {}", x);
    assertion((x(8) == cfloat{1,0}));

  }

  {
    msg("Test abs...");
    soit n = 100;
    soit x = tsd::polar(linspace(0, 2 * π, n));
    soit a = abs(x);

    assertion(a.est_approx(Vecf::ones(n)));
  }

  {
    msg("Test * (Vecf, Veccf)");
    soit x = 10 * Vecf::ones(1);
    soit y = Veccf::ones(1);
    y(0) = cfloat{0.5,2};
    soit z = x * y;

    msg("x={}, y={}, x*y={}", x, y, z);

    assertion((z(0) == cfloat{5,20}));
  }

  {
    msg("Test prod cplx...");
    soit x = Veccf::ones(5);
    soit i = cfloat{0.0f, 1.0f};
    soit y = x * i;
    msg("x = {}, y = x * {} = {}", x, i, y);
    soit z = y * i;
    assertion(z.est_approx(-x));
  }

  {
    msg("Test *(Veccf, float)...");
    soit x = Veccf::ones(2);
    x(0) = cfloat{0.5,2};
    x *= 1.5f;
    assertion((x(0) == cfloat{0.75,3}));
    assertion((x(1) == cfloat{1.5,0}));
  }

  {
    msg("Prod / somme colonnes...");
    soit lsp = linspace(0,2,3);

    soit x = Tabf::eye(3);

    soit y = x.somme_colonnes();
    soit z = x.produit_colonnes(lsp);

    msg("x ={}", x);
    msg("y = {}", y);

    assertion(y.est_approx(Vecf::ones(3)));

    msg("z = {}", z);

    msg("x.moy _ligne = {}", x.moyenne_lignes());

    soit y2 = x.moyenne_colonnes();
    msg("y2 = {}", y2);
    soit err = (y2 - Vecf::ones(3)/3).rms();
    msg("Erreur RMS = {}", err);
    assertion(y2.est_approx(Vecf::ones(3)/3, 1e-7));


    msg("x.ajoute_col(0,1,2) = {}", x.ajoute_colonnes(lsp));

  }

  {
    msg("Manips colonnes...");
    soit x = Tabf::eye(3);
    x.col(1) += linspace(0,2,3);
    assertion(x.col(0).est_approx(Vecf::valeurs({1,0,0})));
    assertion(x.col(1).est_approx(Vecf::valeurs({0,2,2})));
    assertion(x.col(2).est_approx(Vecf::valeurs({0,0,1})));
  }

  {
    msg("dot...");

    soit x1 = Vecf::ones(3),
         x2 = linspace(0,2,3);
    x1(2) = -1;
    soit y  =  x1 + (x1.dot(x2)) * x2;

    msg("y={}", y);
    assertion(x1.dot(x2) == -1);
    assertion(y.est_approx(Vecf::valeurs({1,0,-3})));
  }

  {
    //soit noyau2 = noyau.rowwise().reverse().colwise().reverse();
    soit x = linspace(0,8,9);
    soit X = Tabf::map(x.data(), 3, 3);
    soit xr = x.reverse();
    soit Xr = Tabf::map(xr.data(), 3, 3);

    msg("X = {}, Xr = {}", X, Xr);
  }

  {
    msg("affectation constante...");
    soit x = linspace(0, 9, 10);
    x.setConstant(5);
    assertion(x.est_approx(5*Vecf::ones(10)));
  }

  {
    msg("affectation constante (2)...");
    soit x = linspace(0, 9, 10);
    x = 5.0f;
    msg("Fait.");
    assertion(x.est_approx(5*Vecf::ones(10)));
  }


  {
    msg("Map...");
    float f[9] = {0,1,2,3,4,5,6,7,8};
    soit A = Tabf::map(f, 3, 3);
    msg("map = {}", A);
    assertion((A(0,0) == 0) && (A(1,0) == 1));


    /*soit eB = Eigen::Map<Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic>>(A.data(), 3, 3);
    soit eA = Eigen::Map<Eigen::ArrayXXf>(A.data(), 3, 3);
    msg("eA = \n{}", eA);
    msg("eB = \n{}", eB);*/
  }


  {
    msg("Init par valeurs...");
    soit A = Tabf::valeurs(2, 3,
        {0, 1, 2,
         3, 4, 5});
    msg("A = {}", A);
    assertion((A.rows() == 2) && (A.cols() == 3)
               && (A(0,0) == 0) && (A(1,0) == 3));
  }
}




