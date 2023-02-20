// Simulation de signaux à temps continu

#ifndef CT_HPP
#define CT_HPP

#include "tsd/tsd.hpp"

namespace tsd {





template<typename T>
  using FonctionAbstraite = std::function<T(float)>;

using FonctionRéelle = FonctionAbstraite<float>;


template<typename T>
struct FonctionEchantillonnée
{
  Vecteur<T> data;
  float tmin = -1, tmax = 1;

  static FonctionEchantillonnée<T> def(float tmin, float tmax, int n, FonctionAbstraite<T> fct);


  float dT() const{retourne tmax - tmin;}
  int n()    const{retourne data.rows();}
  float fe() const{retourne (n() - 1.0f) / dT();}
  Vecf t()   const{retourne linspace(tmin, tmax, n());}
};


template<typename T>
FonctionEchantillonnée<T> FonctionEchantillonnée<T>::def(
    float tmin, float tmax, int n, FonctionAbstraite<T> fct)
{
  FonctionEchantillonnée<float> res;
  res.tmin = tmin;
  res.tmax = tmax;
  res.data.resize(n);
  soit vt = res.t();
  pour(auto i = 0; i < n; i++)
    res.data(i) = fct(vt(i));
  retourne res;
}

//using FR = FonctionEchantillonnée<float>;
//using FC = FonctionEchantillonnée<cfloat>;



extern FonctionRéelle fct_impulsion;
extern FonctionRéelle fct_échelon;
extern FonctionRéelle fct_0;
extern FonctionRéelle fct_1;
extern FonctionRéelle fct_sin;


//extern FR fct_impulsion(float tmin, float tmax, int n);

extern FonctionEchantillonnée<cfloat> TF(const FonctionEchantillonnée<float> &f);
extern FonctionEchantillonnée<cfloat> TF(const FonctionAbstraite<float> &f, float fe, int N = 2048);


template<typename T>
T intégrale_trap(const FonctionAbstraite<T> &f, float tmin, float tmax, int N);

template<typename T>
FonctionEchantillonnée<cfloat> tfc(const FonctionAbstraite<T> &fct, float δT, float δF, int nt, int nf);


}

#endif



