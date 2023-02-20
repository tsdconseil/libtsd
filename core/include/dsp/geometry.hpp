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

  nfr::Quaternion fr() const
  {
    return nfr::Quaternion{q};
  }

  Quaternion(const nfr::Quaternion &fr)
  {
    q = fr.q;
  }

  Quaternion(float q0 = 1, float q1 = 0, float q2 = 0, float q3 = 0)
  {
    *this = nfr::Quaternion(q0, q1, q2, q3);
  }
  Quaternion(const Eigen::Array4f &q_)
  {
    *this = nfr::Quaternion(q);
  }
  Quaternion(const Eigen::Matrix3f &R)
  {
    *this = nfr::Quaternion(R);
  }


  /** @brief Rotation inverse
   *
   * @f[
   * q^{-1} = \left(\begin{array}{c}q_0\\ -q_1\\ -q_2 \\-q_3\end{array}\right)
   * @f]
   *
   *  */
  Quaternion inv() const
  {
    return fr().inv();
  }

  /** @brief Element neutre.
   *
 * @f[
 * q = \left(\begin{array}{c}1\\ 0\\ 0 \\0\end{array}\right)
 * @f]
   */
  static Quaternion identite()
  {
    return nfr::Quaternion::identite();
  }

  /** @brief ? */
  Eigen::Matrix4f mat() const
  {
    return fr().mat();
  }

  Eigen::Vector3f rotate(const Eigen::Vector3f &x) const
  {
    return fr().rotate(x);
  }

  /** @brief Matrice de rotation. */
  Eigen::Matrix3f rot_mat() const
  {
    return fr().rot_mat();
  }
};




/** @brief Angles de %Cardan pour la  représentation d'une rotation 3d. */
struct Cardan
{
  /** @brief Roulis (roll) */
  float φ;

  /** @brief Tangage (pitch) */
  float θ;

  /** @brief Dérive (yaw / heading) */
  float ψ;

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

  /** @brief A partir d'une matrice de rotation */
  Cardan(const Eigen::Matrix3f &R)
  {
    *this = nfr::Cardan(R);
  }

  /** @brief A partir des angles */
  Cardan(float φ, float θ, float ψ)
  {
    *this = nfr::Cardan(φ, θ, ψ);
  }

  /** @brief Matrice de rotation. */
  Tabf mat_rotation() const
  {
    return fr().mat_rotation();
  }

  /** @brief A partir d'un quaternion */
  Cardan(const Quaternion &q)
  {
    *this = nfr::Cardan(q.fr());
  }
};


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


template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R1(T α)
{
  return nfr::rotmat_3d_R1(α);
}


template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R2(T α)
{
  return nfr::rotmat_3d_R2(α);
}

template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d_R3(T α)
{
  return nfr::rotmat_3d_R3(α);
}

template<typename T>
  Eigen::Matrix<T, 3, 3> rotmat_3d(T α, int axe)
{
  return nfr::rotmat_3d<T>(α, axe);
}



/** @} */

}

ostream_formater(dsp::geo::Quaternion)
ostream_formater(dsp::geo::Cardan)

