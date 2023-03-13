// Simulation de signaux à temps continu

#ifndef CT_HPP
#define CT_HPP

#include "tsd/tsd.hpp"

namespace tsd {

template<typename T>
  using FonctionAbstraite = std::function<T(float)>;

using FonctionRéelle = FonctionAbstraite<float>;


template<typename T>
struct FonctionEchantillonnée
{
  Vecteur<T> data;
  float tmin = -1, tmax = 1;

  static FonctionEchantillonnée<T> def(float tmin, float tmax, int n, FonctionAbstraite<T> fct);


  float dT() const{retourne tmax - tmin;}
  int n()    const{retourne data.rows();}
  float fe() const{retourne (n() - 1.0f) / dT();}
  Vecf t()   const{retourne linspace(tmin, tmax, n());}
};


template<typename T>
FonctionEchantillonnée<T> FonctionEchantillonnée<T>::def(
    float tmin, float tmax, int n, FonctionAbstraite<T> fct)
{
  FonctionEchantillonnée<float> res;
  res.tmin = tmin;
  res.tmax = tmax;
  res.data.resize(n);
  soit vt = res.t();
  pour(auto i = 0; i < n; i++)
    res.data(i) = fct(vt(i));
  retourne res;
}


extern FonctionRéelle fct_impulsion,
                      fct_échelon,
                      fct_0,
                      fct_1,
                      fct_sin;




/** @brief Intégration (approximation trapézoidale) */
template<typename T>
T intégrale_trap(const FonctionAbstraite<T> &f, float tmin, float tmax, int N);


// To remove?
extern FonctionEchantillonnée<cfloat> TF(const FonctionEchantillonnée<float> &f);
// To remove?
extern FonctionEchantillonnée<cfloat> TF(const FonctionAbstraite<float> &f, float fe, int N = 2048);


/** @brief Approximation de la transformée de Fourier à temps continu.
 *  @param δT L'intervalle temporel est [-δT/2, +δT/2].
 *  @param δF L'intervalle fréquentielle est [-δF/2, +δF/2].
 *  @param nt Nombre de points dans le domaine temporel.
 *  @param nf Nombre de points dans le domaine fréquentiel. */
template<typename T>
FonctionEchantillonnée<cfloat> tfc(const FonctionAbstraite<T> &fct, float δT, float δF, int nt, int nf);


}

#endif



