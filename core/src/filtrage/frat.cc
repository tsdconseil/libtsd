#include "tsd/filtrage/frat.hpp"
#include "tsd/eig-util.hpp"

#include <iostream>

#define ENABLE_ROOTS 1

#if ENABLE_ROOTS
#include <unsupported/Eigen/Polynomials>
#endif

using namespace std;

namespace tsd
{


ostream& operator<<(ostream& os, const Poly<float> &p)
{
  p.display_content(os);
  retourne os;
}
ostream& operator<<(ostream& os, const FRat<float> &f)
{
  f.afficher(os);
  retourne os;
}
ostream& operator<<(ostream& os, const Poly<cfloat> &p)
{
  p.display_content(os);
  retourne os;
}
ostream& operator<<(ostream& os, const FRat<cfloat> &f)
{
  f.afficher(os);
  retourne os;
}

template<typename T>
Veccf Poly<T>::roots() const
{
  si(mode_racines)
  {
    retourne coefs;
  }
  si constexpr(est_complexe<T>())
  {
    échec("Pas possible de calculer les racines d'un polynôme complexe avec cette version de Eigen.");
    retourne Veccf();
  }
  sinon
  {
#   if ENABLE_ROOTS
    Eigen::PolynomialSolver<T, Eigen::Dynamic> solver;

    si(coefs.dim() <= 1)
      retourne Veccf();

    soit v = coefs / coefs.valeur_max();

    // Obligé, sinon la fonction Eigen se bloque (exemple : filtre Gaussien)
    tantque((v.rows() > 1) && (abs(v(v.rows()-1)) < 1e-8)) //== 0.0f))
      v = v.head(v.rows()-1).eval();

    // Polynome nul
    si((v.rows() == 1) && (v(0) == 0.0f))
      retourne Veccf();

    //msg("Calcul racines : v = {}", v);

    Eigen::VectorXf ev = vec2evec(v);
    solver.compute(ev);
    Eigen::VectorXcf r = solver.roots();

    retourne evec2vec(r);

#   else
    msg_erreur("TODO : FRat::roots()");
    retourne ArrayXcf(0);
#   endif
  }
}


template
Veccf Poly<float>::roots() const;

template
Veccf Poly<cfloat>::roots() const;


void un_essai()
{
  Poly<cfloat> pc;
  soit r = pc.roots();
  msg("r = {}", real(r));
}

}
