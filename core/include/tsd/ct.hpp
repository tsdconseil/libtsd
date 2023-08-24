// Simulation de signaux à temps continu

#ifndef CT_HPP
#define CT_HPP

#include "tsd/tsd.hpp"

namespace tsd {


/**  @addtogroup ct
  *  @{ */


/** @brief Fonction à temps continu (type de retour générique). */
template<typename T>
  using FonctionAbstraite = fonction<T(float)>;

/** @brief Fonction à temps continu (type de retour réel). */
using FonctionRéelle = FonctionAbstraite<float>;

/** @brief Approximation d'une fonction à temps continu par échantillonnage. */
template<typename T>
struct FonctionEchantillonnée
{
  /** @brief Valeurs échantillonnées */
  Vecteur<T> data;

  /** @brief Premier point d'échantillonnage */
  float tmin = -1,

  /** @brief Dernier point d'échantillonnage */
        tmax = 1;

  /** @brief Echantillonnage d'une fonction à temps continu */
  static FonctionEchantillonnée<T> def(float tmin, float tmax, entier n, FonctionAbstraite<T> f);

  /** @brief Intervalle temporel (@f$=t_{\textrm{max}} - t_{\textrm{min}}@f$) */
  float dT() const{retourne tmax - tmin;}

  /** @brief Nombre d'échantillons */
  int n()    const{retourne data.rows();}

  /** @brief Fréquence d'échantillonnage */
  float fe() const{retourne (n() - 1.0f) / dT();}

  /** @brief Vecteur des pas de temps = `linspace(tmin, tmax, n())` */
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


// Quelques fonctions classiques

extern FonctionRéelle
  /** @brief Dirac */
  fct_impulsion,
  /** @brief Echelon unité */
  fct_échelon,
  /** @brief Fonction nulle partout */
  fct_0,
  /** @brief Fonction égale à un partout */
  fct_1,
  /** @brief Sinusoide */
  fct_sin;


/** @brief Intégration (approximation trapézoidale)
 *  @param f Fonction à intégrer
 *  @param tmin Borne inférieure de l'intégration
 *  @param tmax Borne supérieure
 *  @param N    Nombre de points à utiliser pour l'approximation
 *
 *  Calcule :
 *  @f[
 *  \int_{t_{\textrm{min}}}^{t_{\textrm{max}}} f(t)\ dt \sim \sum_{k=1}^N ...
 *  @f]
 *
 *   */
template<typename T>
T intégrale_trap(const FonctionAbstraite<T> &f,
                 float tmin, float tmax, entier N = 1000);


/** @cond undoc */
// To remove?
extern FonctionEchantillonnée<cfloat> TF(const FonctionEchantillonnée<float> &f);
// To remove?
extern FonctionEchantillonnée<cfloat> TF(const FonctionAbstraite<float> &f, float fe, int entier = 2048);
/** @endcond */

/** @brief Approximation de la transformée de Fourier à temps continu.
 *  @param δT L'intervalle temporel est [-δT/2, +δT/2].
 *  @param δF L'intervalle fréquentielle est [-δF/2, +δF/2].
 *  @param nt Nombre de points dans le domaine temporel.
 *  @param nf Nombre de points dans le domaine fréquentiel. */
template<typename T>
FonctionEchantillonnée<cfloat> tfc(const FonctionAbstraite<T> &fct, float δT, float δF, entier nt, entier nf);


/**  @} */

}

#endif



