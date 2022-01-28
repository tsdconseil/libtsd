#pragma once

#include "tsd/tsd.hpp"

namespace tsd::tf::cqt {

  /** @addtogroup temps-frequence-cqt
   *  @{
   */


/** @brief Structure de configuration */
struct CQTConfig
{
  /** @brief Fréquence d'échantillonnage (Hz) */
  float fs = 1;
  /** @brief Fréquence d'analyse mini (Hz). */
  float fmin = 0.01;
  /** @brief Fréquence d'analyse maxi (Hz). */
  float fmax = 0.5;
  /** @brief Rapport entre deux fréquences d'analyse successives */
  float γ = std::pow(2.0f, 1.0f/12);
  /** @brief Facteur de qualité */
  float Q = 34;
  /** @brief Précision dans la représentation du noyau */
  float kernel_precision = 0.99;
};


/** @brief Classe pour calculer la %CQT.
 *
 * todo : exemple d'utilisation ici
 *
 *  */
class CQT
{
public:

  CQT();
  void configure(const CQTConfig &config);

  /** @brief Calcul */
  void step(const ArrayXf &x);

  /** @brief Interpolation so as to have a uniformly sampled time / frequency analysis
   *
   * @param ofs Desired output sampling frequency (in Hz)
   * @returns 3 elemetns tuple : time vector, frequency vector, amplitude versus time / frequency (matrix)

   * This function intepolate the result of the CQT so as to have the same number of points,
   * and aligned at the same time, for all frequency bins.
   */
  std::tuple<ArrayXf, ArrayXf, ArrayXXf> interpolation(float ofs);

  void affiche_noyaux();

private:
  struct Impl;
  sptr<Impl> impl;
};





/** @} */



}
