#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/temps.hpp"


#include <Eigen/Core>


namespace tsd::geo {


  using Eigen::Vector3d;

/** @addtogroup geometrie
 *  @{
 */

/** @brief %Quaternion unitaire pour la représentation d'une rotation 3d.
 *
 * @f[
 * q = \left(\begin{array}{c}\cos\theta\\ \sin\theta\cdot \hat e\end{array}\right)
 * @f]
 *
 * @f$\hat e@f$ étant l'axe d'Euler (normalisé).
 *
 *  */
struct Quaternion
{
  /** @brief Les 4 éléments du quaternion. */
  Eigen::Array4f q;

  /** @brief Constructeur (d'après les 4 éléments) */
  Quaternion(float q0 = 1, float q1 = 0, float q2 = 0, float q3 = 0);

  /** @brief Constructeur (d'après les 4 éléments) */
  Quaternion(const Eigen::Array4f &q_);

  /** @brief Constructeur (d'après une matrice de rotation) */
  Quaternion(const Eigen::Matrix3f &R);


  /** @brief Rotation inverse
   *
   * @f[
   * q^{-1} = \left(\begin{array}{c}q_0\\ -q_1\\ -q_2 \\-q_3\end{array}\right)
   * @f]
   *
   *  */
  Quaternion inv() const;

  /** @brief Element neutre.
   *
 * @f[
 * q = \left(\begin{array}{c}1\\ 0\\ 0 \\0\end{array}\right)
 * @f]
   */
  static Quaternion identite();

  /** @brief Matrice de rotation en coordonnées homogènes. */
  Eigen::Matrix4f mat() const;

  /** @brief Applique une rotation */
  Eigen::Vector3f rotate(const Eigen::Vector3f &x) const;



  /** @brief Matrice de rotation. */
  Eigen::Matrix3f rot_mat() const;
};


//extern Eigen::Matrix3f rotmat_3d(float α, entier axe);


/** @brief Angles de %Cardan pour la  représentation d'une rotation 3d. */
struct Cardan
{
  /** @brief Roulis (roll) */
  float φ;

  /** @brief Tangage (pitch) */
  float θ;

  /** @brief Dérive (yaw / heading) */
  float ψ;

  /** @brief Constructeur : à partir d'une matrice de rotation */
  Cardan(const Eigen::Matrix3f &R);

  /** @brief Constructeur : à partir des angles */
  Cardan(float φ, float θ, float ψ);

  /** @brief Constructeur : à partir d'un quaternion */
  Cardan(const Quaternion &q);

  /** @brief Calcule la matrice de rotation correpondante. */
  Tabf mat_rotation() const;
};

/** @cond undoc */

extern std::ostream& operator<<(std::ostream& ss, const Cardan &t);
extern std::ostream& operator<<(std::ostream& ss, const Quaternion &t);

/** @endcond */


/** @brief Matrice de rotation 3d autour de l'axe X */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R1(T α)
{
  auto ca = cos(α), sa = sin(α);
  Eigen::Matrix<T, 3, 3> R;
  R << 1, 0, 0,
       0, ca, sa,
       0, -sa, ca;
  return R;
}


/** @brief Matrice de rotation 3d autour de l'axe Y */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R2(T α)
{
  auto ca = cos(α), sa = sin(α);
  Eigen::Matrix<T, 3, 3> R;
  R << ca, 0, -sa,
       0, 1, 0,
       sa, 0, ca;
  return R;
}

/** @brief Matrice de rotation 3d autour de l'axe Z */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R3(T α)
{
  auto ca = cos(α), sa = sin(α);
  Eigen::Matrix<T, 3, 3> R;
  R << ca, sa, 0,
       -sa, ca, 0,
       0, 0, 1;
  return R;
}

/** @brief Matrice de rotation 3d autour d'un des axes canoniques
 *  @param axe Numéro d'axe (0, 1 ou 2). */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d(T α, entier axe)
{
  if(axe == 0)
    return rotmat_3d_R1<T>(α);
  else if(axe == 1)
    return rotmat_3d_R2<T>(α);
  return rotmat_3d_R3<T>(α);
}


/** @} */

}


ostream_formater(tsd::geo::Quaternion)
ostream_formater(tsd::geo::Cardan)


