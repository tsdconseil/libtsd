#pragma once

/** @file tod.hpp
 *  @brief Transformée en Ondelettes Discrète
 */

#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"
#include <iostream>



namespace tsd::tf::tod {


/** @addtogroup temps-frequence-tod
 *  @{
 */


/** @brief Polynôme de %Laurent (puissance positives et négatives de z)
 *
 * @f[
 * P(z) = \sum_{n=0}^{D} p_n z^{n_0+n}
 * @f]
 * @f$D@f$ étant le degré du polynôme.
 * */
struct Laurent
{
  /** @brief Polynôme */
  tsd::Poly<float> polynome;
  /** @brief Index du premier terme */
  entier n0 = 0;
};


// Chaque étape est un filtre, mais pas forcément causal !!!!

/** @brief Décrit un étage de lifting */
struct LiftElem
{
  /** @brief Polynôme de %Laurent */
  Laurent polynome;
  /** @brief Si vrai, élément de type @f$\left(\begin{array}{cc}1& 0 \\ T(z)& 1\end{array}\right)@f$, sinon
   *  @f$\left(\begin{array}{cc}1 & S(z)\\ 0 & 1\end{array}\right)@f$ */
  bouléen predict = oui;
};

/** @brief Spécification d'une ondelette sous forme d'étapes de lifting.
 *
 *  <h3>Spécification d'une ondelette sous forme d'étapes de lifting</h3>
 *
 *  Série de matrices 2x2 constituée de termes
 *  de type @f$\left(\begin{array}{cc}1 & S(z)\\ 0 & 1\end{array}\right)@f$
 *  ou @f$\left(\begin{array}{cc}1& 0 \\ T(z)& 1\end{array}\right)@f$,
 *  plus un terme de normalisation de la forme :
 *  @f$\left(\begin{array}{cc}K& 0 \\ 0& K^{-1}\end{array}\right)@f$
 */
struct Lift
{
  string nom;
  /** @brief Terme de normalisation */
  float K;
  /** @brief Les différents étages */
  vector<LiftElem> etapes;
};

/** @brief Spécification d'une ondelette sous la forme d'un filtre polyphase (matrice 2x2 de polynômes de Laurent)
 *
 * <h3>Spécification d'une ondelette sous la forme d'un filtre polyphase</h3>
 *
 * @f[
 * H(z) = \left(\begin{array}{cc}
 * H_{00}(z) & H_{01}(z)\\
 * H_{10}(z) & H_{11}(z)
 * \end{array}
 * \right)
 * @f]
 *
 *  */
struct FormePolyphase
{
  /** @brief Développement des matrices implicite de lifting pour obtenir la forme polyphase */
  FormePolyphase(const Lift &lift);
  Laurent H00, H01, H10, H11;
};

/** @brief Spécification d'une ondelette sous forme de deux filtres fonctionnnant en quadrature
 *
 *  QMF = Quadrature Mirror Filter.
 */
struct QMF
{
  /** @brief Conversion forme polyphase @f$\to@f$ filtres QMF */
  QMF(const FormePolyphase &fp);

  /** @brief %Filtre d'analyse (basse fréquence) */
  Poly<float> H0;
  /** @brief %Filtre d'analyse (haute fréquence) */
  Poly<float> H1;
  /** @brief %Filtre de synthèse (basse fréquence) */
  Poly<float> G0;
  /** @brief %Filtre de synthèse (haute fréquence) */
  Poly<float> G1;
};


/** @cond */
/** @brief Affichage d'un polynôme de Laurent */
std::ostream& operator<<(std::ostream& os, const Laurent &p);

/** @brief Affichage d'une forme polyphase */
std::ostream& operator<<(std::ostream& os, const FormePolyphase &p);

/** @brief Affichage d'une forme QMF */
std::ostream& operator<<(std::ostream& os, const QMF &p);
/** @endcond */


/** @brief Spécification de l'ondelette de Haar (sous forme de schéma de lifting) */
extern Lift lift_haar();

/** @brief Spécification de l'ondelette Daubechie d'ordre 2 (sous forme de schéma de lifting) */
extern Lift lift_db2();



// TODO : à cacher ?
/** @brief %Ondelette générique */
template<typename T = float>
struct Ondelette
{
  string nom;
  virtual void lift_step(Vecteur<T> &x, entier n)  = 0;
  virtual void ilift_step(Vecteur<T> &x, entier n) = 0;
};


/** @brief Implémentation d'une ondelette à partir du schèma de lifting
 *
 *  <h3>Implémentation d'une ondelette à partir du schèma de lifting</h3>
 *
 *  @param lift Schèma de lifting
 *
 */
template<typename T>
  sptr<Ondelette<T>> ondelette_gen(const Lift &lift);




/** @brief Transformée en ondelette
 *
 * <h3>Transformée en ondelette</h3>
 *
 * @param ondelette Ondelette abstraite
 * @param[inout] x Données à traiter
 * @param profondeur Nombre d'étages de transformée
 *
 */
template<typename T>
  void dwt(sptr<Ondelette<T>> ondelette, Vecteur<T> &x, entier profondeur);

/** @brief Transformée en ondelette inverse
 *
 * <h3>Transformée en ondelette inverse</h3>
 *
 * @param ondelette Ondelette abstraite
 * @param[inout] x Données à traiter
 * @param profondeur Nombre d'étages de transformée
 */
template<typename T>
  void iwt(sptr<Ondelette<T>> ondelette, Vecteur<T> &x, entier profondeur);



/** @} */



}


ostream_formater(tsd::tf::tod::Laurent)
ostream_formater(tsd::tf::tod::FormePolyphase)
ostream_formater(tsd::tf::tod::QMF)


