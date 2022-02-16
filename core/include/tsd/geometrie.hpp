#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/temps.hpp"



namespace tsd::geo {


  using Eigen::Vector3d;

/** @addtogroup geometrie
 *  @{
 */

/** @brief %Quaternion unitaire pour la représentation d'une rotation 3d.
 *
 * <h3>%Quaternion unitaire</h3>
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

  Quaternion(float q0 = 1, float q1 = 0, float q2 = 0, float q3 = 0);
  Quaternion(const Eigen::Array4f &q_);
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

  /** @brief ? */
  Eigen::Matrix4f mat() const;

  Eigen::Vector3f rotate(const Eigen::Vector3f &x) const;



  /** @brief Matrice de rotation. */
  Eigen::Matrix3f rot_mat() const;
};


extern Eigen::Matrix3f rotmat_3d(float α, int axe);


/** @brief Angles de %Cardan pour la  représentation d'une rotation 3d. */
struct Cardan
{
  /** @brief Roulis (roll) */
  float φ;

  /** @brief Tangage (pitch) */
  float θ;

  /** @brief Dérive (yaw / heading) */
  float ψ;

  /** @brief A partir d'une matrice de rotation */
  Cardan(const Eigen::Matrix3f &R);

  /** @brief A partir des angles */
  Cardan(float φ, float θ, float ψ);

  /** @brief Matrice de rotation. */
  Eigen::Matrix3f mat_rotation() const;

  /** @brief A partir d'un quaternion */
  Cardan(const Quaternion &q);
};


extern std::ostream& operator<<(std::ostream& ss, const Cardan &t);
extern std::ostream& operator<<(std::ostream& ss, const Quaternion &t);


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


/** @} */

}




