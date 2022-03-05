#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"

namespace tsd {


/** @addtogroup divers
 *  @{
 */


/** @brief TF d'une fonction porte entre @f$-T/2@f$ et @f$T/2@f$ (fonction sinus cardinal).
 *
 *  <h3>TF d'une fonction porte</h3>
 *
 *  @param T Largeur de la porte dans le domaine temporel
 *  @param f Fréquence souhaitée pour l'évaluation de la TF
 *  @returns Valeur de la TF de la fonction porte en @f$f@f$ :
 *  @f[
 *   y = \mathcal{F}\left(\Pi_T\right)(f)
 *  @f]
 *  @f$\Pi_T@f$ étant une porte de largeur @f$T@f$ :
 *  @f[
 *  \Pi_T(t) = \begin{cases} 1 & \textrm{ si } -T/2 \leq t \leq T/2,\\ 0 & \textrm{sinon.} \end{cases}
 *  @f]
 *
 *  Au facteur d'échelle @f$T@f$ près, cette fonction n'est autre que le sinus cardinal :
 *  @f[
 *  \textrm{sinc}_T(f) = \frac{\sin \pi T f}{\pi f}
 *  @f]
 *
 *  @note Notez que cette fonction peut aussi être interprétée comme la TF inverse d'une porte fréquentielle.
 *  Dans ce cas, @f$T@f$ doit être interprété comme <b>deux fois</b> la fréquence de coupure @f$f_c@f$,
 *  et @f$f@f$ comme le temps @f$t@f$.
 *
 *  @par Exemple 1 : TF d'une fonction porte de largeur 1 (+/- 1/2)
 *  @code
 *  y = sinc(1, f);
 *  @endcode
 *
 *  @par exemple 2 : TFI d'une fonction porte de largeur 2fc
 *  @code
 *  y = sinc(2 * fc, t);
 *  @endcode
 *
 *  @sa sinc()
 *
 */
extern float sinc(float T, float f);

/** @brief Sinus cardinal normalisé avec fréquence de coupure 0,5.
 *
 *  <h3>Sinus cardinal normalisé</h3>
 *  Calcul d'un sinus cardinal (fréquence de coupure 0,5) :
 *  @f[
 *  y(t) = \frac{\sin(\pi t)}{\pi  t}
 *  @f]
 *
 *  @param t Point d'échantillonnage (en nombre d'échantillons)
 *  @returns Valeur de @f$y@f$
 *
 *  @sa sinc()
 *
 *  */
extern float sinc(float t);


/** @brief Noyau de Dirichlet (à documenter !) */
extern float Dirichlet(int N, float Ω);

/** @brief Polynôme de Chebychev du premier type : @f$T_n\left(\cos(\theta)\right) = \cos(n\theta)@f$ */
extern Poly<float> Chebychev_T(int n);

/** @brief Polynôme de Chebychev du deuxième type : @f$U_n\left(\sin(\theta)\right) = \sin(n\theta)@f$ */
extern Poly<float> Chebychev_U(int n);

/** @} */

}
