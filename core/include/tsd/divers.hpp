#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"

namespace tsd {


/** @addtogroup divers
 *  @{
 */


/** @brief Transformée de Fourier d'une fonction porte entre @f$-T/2@f$ et @f$T/2@f$ (fonction sinus cardinal).
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
 *  @note Notez que cette fonction peut aussi être interprétée comme la Transformée de Fourier inverse d'une porte fréquentielle.
 *  Dans ce cas, @f$T@f$ doit être interprété comme <b>deux fois</b> la fréquence de coupure @f$f_c@f$,
 *  et @f$f@f$ comme le temps @f$t@f$.
 *
 *  @par Exemple 1 : TF d'une fonction porte temporelle
 *  Ici, on calcule la TF d'une porte temporelle d'extension @f$\pm 1/2@f$,
 *  ce qui correspond par exemple à la réponse fréquentielle d'une moyenne glissante.
 *  @code
 *  y = sinc(1, f);
 *  @endcode
 *
 *  @par Exemple 2 : TFI d'une porte fréquentielle
 *  Ici, on calcule la TFI d'une porte fréquentielle de largeur @f$2f_c@f$ (extension @f$\pm f_c@f$),
 *  c'est-à-dire un filtre idéal coupant à la fréquence @f$f_c@f$.
 *  @code
 *  y = sinc(2 * fc, t);
 *  @endcode
 *
 *  @sa sinc(float)
 *
 */
extern float sinc(float T, float f);

/** @brief Sinus cardinal normalisé avec fréquence de coupure 0,5.
 *
 *  Cette fonction est un cas particulier de sinus cardinal (@ref sinc(float, float)), pour
 *  lequel la porte est d'extension 1 (@f$\pm 1/2@f$)~:
 *  @f[
 *  \textrm{sinc}(t) = \textrm{sinc}(1,t) = \frac{\sin(\pi t)}{\pi  t}
 *  @f]
 *
 *  @param t Point d'échantillonnage (en nombre d'échantillons)
 *  @returns Valeur de @f$y@f$
 *
 *  @sa sinc(float, float)
 *
 *  */
extern float sinc(float t);


/** @brief Noyau de Dirichlet (sinus cardinal périodique) (à documenter !)
 *
 *  @f[
 *  D_N(\Omega) = \sum_{k=-N}^N e^{-\mathbf{i}\Omega k}
 *  @f]
 *
 * */
extern double Dirichlet(entier N, double Ω);

/** @brief Polynôme de Chebychev du premier type : @f$T_n\left(\cos(\theta)\right) = \cos(n\theta)@f$
 *
 * Calculé d'après la récursion :
 * @f[
 *   T_n = 2 z \cdot T_{n-1} - T_{n-2}
 * @f]
 * et :
 * @f[
 *   T_0 = 1,\ T_1 = z
 * @f]
 *
 * @sa Chebychev_U
 *
 */
extern Poly<float> Chebychev_T(entier n);

/** @brief Polynôme de Chebychev du deuxième type : @f$U_n\left(\sin(\theta)\right) = \sin(n\theta)@f$
 *
 * Calculé d'après la récursion :
 * @f[
 *   T_n = 2 z \cdot T_{n-1} - T_{n-2}
 * @f]
 * et :
 * @f[
 *   T_0 = 1,\ T_1 = 2z
 * @f]
 *
 * @sa Chebychev_T
 *
 */
extern Poly<float> Chebychev_U(entier n);


//extern Poly<float> BesselRev(entier n);

/** @} */

}
