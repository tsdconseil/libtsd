#include "tsd/tsd-all.hpp"
#include "tsd/eig-util.hpp"
#include "Eigen/LU"


/* Matrice de covariance par méthode directe (lente)
 * m : dimension de la matrice (nombre de délais examinés) */
template<typename T>
Tabf cov(const Vecteur<T> &x, entier m)
{
  soit n = x.rows();
  soit R = Tabf::zeros(m, m);
  pour(auto i = m; i < n; i++)
  {
    soit t = x.segment(i - m, m).reverse();
    R += t.matprod(t.adjoint());
  }
  retourne R / n;
}


static void test_cov()
{

}



void test_levinson()
{
  msg_majeur("Test levinson - résolution R a = -r...");

  soit p = 5;
  soit r = Vecf::random(p);
  soit R = r_vers_R(r.head(p-1));

  soit ac = levinson_reel(r);
  Eigen::VectorXf tmp = tab2etab(R).lu().solve(-vec2evec(r).tail(p-1));
  soit ar  = evec2vec(tmp);
  soit ac2 = levinson(r.head(p-1), r.head(p-1), -r.tail(p-1));

  tsd_assert(ar.rows() == p - 1);
  tsd_assert(ac.rows() == p);
  tsd_assert(ac2.rows() == p - 1);


  soit err     = abs(ac.tail(p-1) - ar).valeur_max(),
       erry    = abs((R.matprod(ac.tail(p-1))) + r.tail(p-1)).valeur_max(),
       erryref = abs((R.matprod(ar)) + r.tail(p-1)).valeur_max(),
       erry2   = abs((R.matprod(ac2)) + r.tail(p-1)).valeur_max();

  msg("Erreur solution : {}, prédiction LU : {}, inversion levin : {}, inversion levin gen : {}",
      err, erryref, erry, erry2);

  si((err > 5e-5) || (erry > 1e-5) || (erryref > 1e-5) || (erry2 > 1e-5))
    echec("Echec levinson.");


  msg("ok.");
}

entier test_stats()
{
  test_cov();
  test_levinson();
  retourne 0;
}


