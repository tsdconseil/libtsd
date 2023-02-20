/*
 * eig-util.hpp
 *
 *  Created on: 20 sept. 2022
 *      Author: julien
 */

#ifndef CORE2_INCLUDE_TSD_EIG_UTIL_HPP_
#define CORE2_INCLUDE_TSD_EIG_UTIL_HPP_


#include "tsd/tableau.hpp"
#include "Eigen/Core"

namespace tsd {


template<typename T, int dim>
  auto evec2vec(const Eigen::Matrix<T,dim,1> &x)
{
  Vecteur<T> y(x.rows());
  memcpy(y.data(), x.data(), x.rows() * sizeof(T));
  return y;
}


inline auto evec2vec(const Eigen::Ref<const Eigen::ArrayXf> &x)
{
  Vecteur<float> y(x.rows());
  memcpy(y.data(), x.data(), x.rows() * sizeof(float));
  return y;
}

/*template<typename T, int dim>
  auto evec2vec(Eigen::Ref<const Eigen::Array<T,dim,1>> &x)
{
  Vecteur<T> y(x.rows());
  memcpy(y.data(), x.data(), x.rows() * sizeof(T));
  return y;
}

template<typename T, int dim>
  auto evec2vec(const Eigen::Array<T,dim,1> &x)
{
  Vecteur<T> y(x.rows());
  memcpy(y.data(), x.data(), x.rows() * sizeof(T));
  return y;
}*/

template<typename T>
  auto vec2evec(const Vecteur<T> &x)
{
  return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1>>(x.data(), x.dim());
}

template<typename T>
  auto etab2tab(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &x)
{
  TabT<T,2> y(x.rows(), x.cols());
  memcpy(y.data(), x.data(), x.rows() * x.cols() * sizeof(T));
  return y;
}

template<typename T>
  auto tab2etab(const TabT<T,2> &x)
{
  return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>(x.data(), x.rows(), x.cols());
}


}



#endif /* CORE2_INCLUDE_TSD_EIG_UTIL_HPP_ */
