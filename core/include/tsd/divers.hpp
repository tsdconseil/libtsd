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
 *  @param ω Pulsation souhaitée pour l'évaluation de la TF
 *  @returns Valeur de la TF de la fonction porte en ω :
 *  @f[
 *   y = \mathcal{F}\left(\Pi_T\right)(\omega)
 *  @f]
 *  @f$\Pi_T@f$ étant une porte de largeur @f$T@f$ :
 *  @f[
 *  \Pi_T(t) = \begin{cases} 1 & \textrm{ si } -T/2 \leq t \leq T/2,\\ 0 & \textrm{sinon.} \end{cases}
 *  @f]
 *
 *  A un facteur d'échelle près, cette fonction n'est autre que le sinus cardinal :
 *  @f[
 *  \textrm{sinc}(T, \omega) = \frac{\sqrt{2}}{\pi}\cdot \frac{\sin T\omega/2}{\omega}
 *  @f]
 *
 */
extern float sinc(float T, float ω);

/** @brief Noyau de Dirichlet (à documenter !) */
extern float Dirichlet(int N, float Ω);

/** @brief Polynôme de Chebychev du premier type : @f$T_n\left(\cos(\theta)\right) = \cos(n\theta)@f$ */
extern Poly<float> Chebychev_T(int n);

/** @brief Polynôme de Chebychev du deuxième type : @f$U_n\left(\sin(\theta)\right) = \sin(n\theta)@f$ */
extern Poly<float> Chebychev_U(int n);

/** @} */

}
