#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "dsp/dsp.hpp"
#include "dsp/time.hpp"
#include "tsd/geometrie.hpp"


namespace dsp::geo {


  namespace nfr = tsd::geo;

  using Eigen::Vector3d;

/** @addtogroup geometrie
 *  @{
 */

/** @brief Unitary quaternion for the representation of a 3d rotation.
 *
 * @f[
 * q = \left(\begin{array}{c}\cos\theta\\ \sin\theta\cdot \hat e\end{array}\right)
 * @f]
 *
 * @f$\hat e@f$ being the normalized Euler axis.
 *
 */
struct Quaternion
{
  /** @brief The 4 elements of the quaternion. */
  Eigen::Array4f q;

  /** @cond undoc */
  nfr::Quaternion fr() const
  {
    return nfr::Quaternion{q};
  }


  Quaternion(const nfr::Quaternion &fr)
  {
    q = fr.q;
  }
  /** @endcond */

  /** @brief Constructor (from the 4 elements). */
  Quaternion(float q0 = 1, float q1 = 0, float q2 = 0, float q3 = 0)
  {
    *this = nfr::Quaternion(q0, q1, q2, q3);
  }

  /** @brief Constructor (from the 4 elements). */
  Quaternion(const Eigen::Array4f &q_)
  {
    *this = nfr::Quaternion(q);
  }

  /** @brief Constructor (from a rotation matrix). */
  Quaternion(const Eigen::Matrix3f &R)
  {
    *this = nfr::Quaternion(R);
  }


  /** @brief Inverse rotation
   *
   * @f[
   * q^{-1} = \left(\begin{array}{c}q_0\\ -q_1\\ -q_2 \\-q_3\end{array}\right)
   * @f]
   *
   */
  Quaternion inv() const
  {
    return fr().inv();
  }

  /** @brief Neutral element.
   *
   * @f[
   * q = \left(\begin{array}{c}1\\ 0\\ 0 \\0\end{array}\right)
   * @f]
   */
  static Quaternion identite()
  {
    return nfr::Quaternion::identite();
  }

  /** @brief Rotation matrix in homogeneous coordinates. */
  Eigen::Matrix4f mat() const
  {
    return fr().mat();
  }

  /** @brief Apply a rotation */
  Eigen::Vector3f rotate(const Eigen::Vector3f &x) const
  {
    return fr().rotate(x);
  }

  /** @brief Computes the rotation matrix. */
  Eigen::Matrix3f rot_mat() const
  {
    return fr().rot_mat();
  }
};

/** @brief %Cardan angles for the representation of a 3d rotation. */
struct Cardan
{
  /** @brief Roll */
  float φ;

  /** @brief Pitch */
  float θ;

  /** @brief Yaw / heading */
  float ψ;

  /** @cond undoc */
  nfr::Cardan fr() const
  {
    return nfr::Cardan(φ, θ, ψ);
  }

  Cardan(const nfr::Cardan &f)
  {
    φ = f.φ;
    θ = f.θ;
    ψ = f.ψ;
  }
  /** @endcond */

  /** @brief From a rotationmatrix. */
  Cardan(const Eigen::Matrix3f &R)
  {
    *this = nfr::Cardan(R);
  }

  /** @brief From Cardan angles. */
  Cardan(float φ, float θ, float ψ)
  {
    *this = nfr::Cardan(φ, θ, ψ);
  }

  /** @brief From a quaternion. */
  Cardan(const Quaternion &q)
  {
    *this = nfr::Cardan(q.fr());
  }

  /** @brief Compute the rotation matrix. */
  Tabf mat_rotation() const
  {
    return fr().mat_rotation();
  }


};

/** @cond undoc */
inline std::ostream& operator<<(std::ostream& ss, const Cardan &t)
{
  ss << t.fr();
  return ss;
}

inline std::ostream& operator<<(std::ostream& ss, const Quaternion &t)
{
  ss << t.fr();
  return ss;
}
/** @endcond */


/** @brief 3D rotation matrix, around X axis. */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R1(T α)
{
  return nfr::rotmat_3d_R1(α);
}

/** @brief 3D rotation matrix, around Y axis. */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R2(T α)
{
  return nfr::rotmat_3d_R2(α);
}

/** @brief 3D rotation matrix, around Z axis. */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R3(T α)
{
  return nfr::rotmat_3d_R3(α);
}

/** @brief 3D rotation matrix, around one of the 3 canonical axis.
 *  @param axe  Axis number (0, 1 or 2). */
template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d(T α, int axe)
{
  return nfr::rotmat_3d<T>(α, axe);
}



/** @} */

}

ostream_formater(dsp::geo::Quaternion)
ostream_formater(dsp::geo::Cardan)

