#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "dsp/dsp.hpp"
#include "tsd/divers.hpp"
#include "dsp/filter.hpp"

namespace dsp {


/** @addtogroup divers
 *  @{
 */


/** @brief FT of door function between @f$-T/2@f$ and @f$T/2@f$ (sinc function).
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
static inline float sinc(float T, float f)
{
  return tsd::sinc(T, f);
}

/** @brief Normalized cardinal sine with cut-off frequency at 0.5.
 *
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
static inline float sinc(float t)
{
  return tsd::sinc(t);
}


/** @brief Dirichlet kernel (periodic cardinal sinus)
 *
 *  @f[
 *  D_N(\Omega) = \sum_{k=-N}^N e^{-\mathbf{i}\Omega k}
 *  @f]
 */
static inline double Dirichlet(int N, double Ω)
{
  return tsd::Dirichlet(N, Ω);
}

/** @brief Chebychev polynomial of the first type: @f$T_n\left(\cos(\theta)\right) = \cos(n\theta)@f$
 *
 * Computed from the recursion:
 * @f[
 *   T_n = 2 z \cdot T_{n-1} - T_{n-2}
 * @f]
 * and:
 * @f[
 *   T_0 = 1,\ T_1 = z
 * @f]
 *
 * @sa Chebychev_U()
 */
static inline Poly<float> Chebychev_T(int n)
{
  return tsd::Chebychev_T(n);
}

/** @brief Chebychev polynomial of the second type: @f$U_n\left(\sin(\theta)\right) = \sin(n\theta)@f$
 *
 * Computed from the recursion:
 * @f[
 *   T_n = 2 z \cdot T_{n-1} - T_{n-2}
 * @f]
 * and:
 * @f[
 *   T_0 = 1,\ T_1 = 2z
 * @f]
 *
 * @sa Chebychev_T()
 */
static inline Poly<float> Chebychev_U(int n)
{
  return tsd::Chebychev_U(n);
}

/** @} */

}
