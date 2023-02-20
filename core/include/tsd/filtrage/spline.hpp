#pragma once

#include "tsd/tsd.hpp"

#include <vector>
#include <assert.h>

namespace tsd::filtrage {






/** @brief Compute the cardinal spline filter coefficients.
 *  @param c Tension parameter (0 : Catmull-Rom spline)
 *  @param t Delay between 0 and 1.0
 *  @param[out] g The four filter coefficents
 *
 *  ==> Now to interpolate @ t, knowing the 4 neighboring points :
 *
 * --> We know p(0), p(1)
 * --> We have to choose p'(0), p'(1)
 * --> Cardinal spline => We choose p'(i) = (1 - c) * 0.5 * (p(i+1) - p(i-1)).
 * With c = tension parameter.
 * c = 1 ==> all tangents = 0
 * c = 0 ==> Catmull-Rom spline
 * Thus, to interpolate @t, between 1 and 2:
 *       p'(1) = (1 - c) * 0.5 * (p(2) - p(0))
 *       p'(2) = (1 - c) * 0.5 * (p(3) - p(1))
 *
 * ==> p(t) = h(t,0) * p(1) + h(t,1) * p'(1) + h(t,2) * p(2) + h(t,3) * p'(2)
 * Wich gives, factoring by p(0), p(1), p(2), p(2):
 *
 *     p(t) = p(0) * -(1 - c) * 0.5 * h(t,1)
 *          + p(1) * (h(t,0) - h(t,3) * (1 - c) * 0.5)
 *          + p(2) * (h(t,2) + h(t,1) * (1 - c) * 0.5)
 *          + p(3) * (h(t,3) * (1 - c) * 0.5)
 */
extern Vecf cspline_filtre(float t, float c);


/** @brief Compute a LUT for cardinal spline interpolation
 *  This LUT can be used as a polyphase interpolation filter.
 *  @param n Number of delayed version of the filter (time resolution).
 *  @param c Tension parameter (0 : Catmull-Rom spline)
 *  @param[out] lut T-type matrix, with 4 rows, and n+1 columns.
 *  Each column contains the 4 coefficients of the filter for a given phase.
 *  @tparam T Type for the LUT elements.
 *
 *  Example: for n = 256, compute a 257 entries LUT,
 *  with entry 0   = coefficients for t = 0,
 *             1       ...          t = 1 / 256,
 *             ...
 *             256     ...          t = 1.0
 */
// n = rapport de sur-échantillonage
extern Tabf cspline_calc_lut(entier n, float c = 0);




} // namespace tsd::filtrage

