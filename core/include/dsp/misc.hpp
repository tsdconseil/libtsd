#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "dsp/dsp.hpp"
#include "tsd/filtrage/frat.hpp"

namespace dsp {


/** @addtogroup divers
 *  @{
 */


/** @brief FT of door function between @f$-T/2@f$ and @f$T/2@f$ (sinc function).
 *
 *  <h3>FT of door function</h3>
 *
 *  @param T Door width in time domain.
 *  @param f Desired frequency for FT evaluation.
 *  @returns Value of the FT at frequency @f$f@f$:
 *  @f[
 *   y = \mathcal{F}\left(\Pi_T\right)(f)
 *  @f]
 *  @f$\Pi_T@f$ being a door of width @f$T@f$:
 *  @f[
 *  \Pi_T(t) = \begin{cases} 1 & \textrm{ if } -T/2 \leq t \leq T/2,\\ 0 & \textrm{otherwise.} \end{cases}
 *  @f]
 *
 *  Up to the scale factor @f$T@f$, this function is nothing else as the classical sinc function:
 *  @f[
 *  \textrm{sinc}_T(f) = \frac{\sin \pi T f}{\pi f}
 *  @f]
 *
 *  @note Note that this function can also be interpreted as the inverse FT of a frequential door function.
 *  In this case, @f$T@f$ must be interpreted as <b>twice</b> the cutt-off frequency @f$f_c@f$,
 *  and @f$f@f$ as the time @f$t@f$.
 *
 *  @par Example 1: FT of temporal door of width 1 (+/- 1/2)
 *  @code
 *  y = sinc(1, f);
 *  @endcode
 *
 *  @par Example 2: IFT of frequency door of width 2fc
 *  @code
 *  y = sinc(2 * fc, t);
 *  @endcode
 *
 *  @sa sinc()
 *
 */
extern float sinc(float T, float f);

/** @brief Normalized cardinal sine with cut-off frequency at 0.5.
 *
 *  <h3>Normalized cardinal sine</h3>
 *  Compute a cardinal sine with cut-off frequency at 0.5:
 *  @f[
 *  y(t) = \frac{\sin(\pi t)}{\pi  t}
 *  @f]
 *
 *  @param t Temporal sampling point (as fractionnal number of samples).
 *  @returns Normalized cut-off frequency, between 0 and 0.5.
 *
 *  @sa sinc2()
 *
 */
extern float sinc(float t);


/** @brief Noyau de Dirichlet (à documenter !) */
extern float Dirichlet(int N, float Ω);

/** @brief Polynôme de Chebychev du premier type : @f$T_n\left(\cos(\theta)\right) = \cos(n\theta)@f$ */
extern Poly<float> Chebychev_T(int n);

/** @brief Polynôme de Chebychev du deuxième type : @f$U_n\left(\sin(\theta)\right) = \sin(n\theta)@f$ */
extern Poly<float> Chebychev_U(int n);

/** @} */

}
