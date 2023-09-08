#include "tsd/ct.hpp"
#include "tsd/tsd-all.hpp"


namespace tsd
{


FonctionRéelle
  fct_impulsion = Λ(x, abs(x) < 1e-6 ? 1 : 0),
  fct_échelon   = Λ(x, x < 0 ? -1 : ((x == 0) ? 0 : 1)),
  fct_0         = Λ(x, 0),
  fct_1         = Λ(x, 1),
  fct_sin       = Λ(x, sin(x));


template<typename TT>
  struct augmente_precisison {};

template<>
  struct augmente_precisison<float> {using type = double;};
template<>
  struct augmente_precisison<cfloat> {using type = cdouble;};
template<>
  struct augmente_precisison<double> {using type = double;};
template<>
  struct augmente_precisison<cdouble> {using type = cdouble;};


template<typename T>
T intégrale_trap(const FonctionAbstraite<T> &f, float tmin, float tmax, int n)
{
  using Ti = augmente_precisison<T>::type;
  Ti y = 0;
  soit δ = (tmax - tmin) / (n-1);
  soit t = linspace(tmin, tmax, n);

  soit f0 = f(t(0));

  si(std::isnan(real(f0)))
  {
    échec("intégrale trap: f0=nan, t(0) = {}.", t(0));
  }

  pour(auto i = 0; i + 1 < n; i++)
  {
    soit f1 = f(t(i+1));

    si(std::isnan(real(f1)))
    {
      échec("intégrale trap: f1={}:{}, t = {}, i ={}.", real(f1), imag(f1), t(i+1), i);
    }

    si(std::isnan(imag(f1)))
    {
      échec("intégrale trap: f1={}:{}, t = {}, i ={}.", real(f1), imag(f1), t(i+1), i);
    }


    //y += δ * (f0 + (f1 - f0) / 2);
    y += ((T) δ) * (f0 + f1) * ((T) 0.5);
    f0 = f1;
  }

  si(std::isnan(real(y)))
  {
    échec("intégrale trap: nan.");
  }

  retourne (T) y;
}




FonctionEchantillonnée<cfloat> TF(const FonctionAbstraite<float> &f, float fe, int N)
{
  soit T  = N / fe;
  soit fec = FonctionEchantillonnée<float>::def(-T/2, T/2, N, f);
  retourne TF(fec);
}

FonctionEchantillonnée<cfloat> TF(const FonctionEchantillonnée<float> &f)
{
  // Approximation de la TF
  // Déformations de fft(f.data) :
  //  (1) échantillonage de la TF (période = 1/T) -> on ne peut pas y faire grand chose...
  //  (2) durée finie                 -> filtrage sinc -> pré-fenêtrage ?
  //  (3) échantillonage de t         -> réplication du spectre
  //  (4) signal décalé dans le temps -> modulation dans le domaine fréquentiel
  //      (mais pas d'impact sur l'amplitude)

  soit n = f.n();
  soit X = fft(f.data);

  // Corrige le décalage dans le temps
  // Position de t = 0, en nombre d'échantillons
  float l = -(f.tmin * (n-1)) / (f.tmax - f.tmin);
  //msg("l = {}", l);
  X *= tsd::polar(l * 2 * π * linspace(0,n-1,n) / n);

  // Corrige le décalage fréquentiel
  X = fftshift(X);

  FonctionEchantillonnée<cfloat> F;
  F.data = X;
  F.tmin = -f.fe() / 2;
  F.tmax = f.fe() / 2;

  retourne F;
}

template<typename T>
FonctionEchantillonnée<cfloat> tfc(const FonctionAbstraite<T> &fct, float δT, float δF, int nt, int nf)
{
  FonctionEchantillonnée<cfloat> res;
  soit f = linspace(-δF/2, δF/2, nf);

  Veccf H(nf);
  pour(auto i = 0; i < nf; i++)
  {
    soit tf = [i,f,fct](float x) -> cfloat
    {
      soit v = fct(x);

      if constexpr(est_complexe<T>())
     {
      si(std::isnan(real(v)))
      {
        échec("tfc: nan @ x = {}", x);
      }
     }
      retourne v * std::polar(1.0f, -2*π_f*x*f(i));
    };
    H(i) = intégrale_trap<cfloat>(tf, -δT/2, δT/2, nt);
  }
  res.data = H;
  res.tmin = -δF/2;
  res.tmax = δF/2;
  retourne res;
}




namespace caché {
soit intégrale_trap1 = intégrale_trap<float>;
soit intégrale_trap2 = intégrale_trap<cfloat>;
soit tfc1 = tfc<float>;
soit tfc2 = tfc<cfloat>;
}

}
