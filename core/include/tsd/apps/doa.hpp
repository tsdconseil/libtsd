#pragma once

#include "tsd/tsd.hpp"
#include "tsd/stats.hpp"
#include "tsd/fourier.hpp"

namespace tsd::apps::doa {



  /** @addtogroup apps-doa
   *  @{
   */



  /** @brief DOA linéaire 1d, avec micro / antennes équidistantes.
   *
   *  Calcul des angles d'incidences (DOA), pour un réseau de micros / antennes linéaire 1d et équidistantes.
   *
   *  @param R  Matrice d'auto-corrélation (dimension = nombre d'antennes * nombre d'antennes)
   *  @param d  Distance entre deux antennes (en unité de longueur d'onde de la fréquence porteuse)
   *  @param Ns Nombre de sources (ou -1 si détermination automatique)
   *  @returns  Vecteur des angles d'indicence de chaque source */
  extern ArrayXf musicdoa_1d(const MatrixXcf &R, float d, int Ns);


  /** @brief Calcul des steering vectors théoriques.
   *
   *  Calcul des steering vectors théoriques pour un réseau de micros / antennes linéaire 1d et équidistantes.
   *  @param pos Position des antennes (unité = longueur d'onde de la porteuse)
   *  @param angles Angles d'incidences des différentes sources (en radians)
   *  @returns Matrice des steering vectors (chaque colonne de la matrice correspond à une source, chaque ligne à une antenne)
   */
  extern MatrixXcf steervec_1d(const ArrayXf &pos, const ArrayXf &angles);


  /** @brief Calcul de la matrice de covariance théorique.
   *
   *  Calcul de la matrice de covariance théorique pour un réseau de micros / antennes linéaire 1d et équidistantes.
   *
   *  @param pos Position des antennes (unité = longueur d'onde de la porteuse)
   *  @param angles Angles d'incidences des différentes sources (en radians)
   *  @param SNR_db Rapport signal à bruit de chaque source (on suppose ici qu'elles ont toutes le même SNR)
   *  @returns Matrice de covariance (Nr par Nr, Nr étant le nombre d'antennes).
   */
  extern MatrixXcf sensorcov_1d(const ArrayXf &pos, const ArrayXf &angles, float SNR_db);


  /** @} */






}




