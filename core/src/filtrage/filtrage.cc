#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <set>

using namespace std;
using namespace tsd::fourier;

namespace tsd::filtrage
{

  bouléen debug_design = non;

  Vecf design_rif_diff(entier n, bouléen centré)
  {
    si(n < 2)
      échec("design_rif_diff(): n doit être supérieur ou égal à 2 (ici {})", n);
    sinon si(n == 2)
      retourne Vecf::valeurs({-1, 1});


    Vecf h(n);

    si(!centré)
    {
      pour(auto i = 0; i < n; i++)
      {
        soit k = i - n/2;
        h(i) = (((k + 1) & 1) ? -1 : 1) / (π * (0.5 + k) * (0.5 + k));
      }
    }
    sinon
    {
      pour(auto i = 0; i < n; i++)
      {
        soit k = i - n/2;
        si(k == 0)
          h(i) = 0;
        sinon
          h(i) = (((k) & 1) ? -1.0 : 1.0) / k;
      }
    }

    retourne h;
  }

  Vecf design_rif_prod(const Vecf &h1, const Vecf &h2)
  {
    soit n1  = h1.rows();
    soit h2p = vconcat(h2, Vecf::zeros(n1-1));
    retourne filtrer(h1, h2p);
  }

  Vecf design_rif_pb2ph_rs(const Vecf &h)
  {
    // Attention, nécessite que h soit de type I
    assertion_msg(type_rif(h) == 1, "design_rif_pb2ph_rs: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = h.clone();
    //soit m = (n-1) / 2;

    pour(auto k = 0; k < n; k++)
      si(est_impair(k))
        g(k) = -g(k);

    /*for(auto k = 0; k < m; k++)
    {
      si(est_impair(k))
      {
        g(m+k) = -g(m+k);
        g(m-k) = -g(m-k);
      }
    }*/

    retourne g;
  }

  Vecf design_rif_pb2ph_is(const Vecf &h)
  {
    // Attention, nécessite que h soit de type I
    assertion_msg(type_rif(h) == 1, "design_rif_pb2ph_is: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;
  }

  Vecf design_rif_pb2pb(const Vecf &h)
  {
    /*assertion_msg(type_rif(h) == 1, "design_rif_pb2ph: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;*/
    retourne h;
  }

  Vecf design_rif_pb2pm(const Vecf &h)
  {

    /*assertion_msg(type_rif(h) == 1, "design_rif_pb2ph: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;*/
    retourne h;
  }


  void verifie_frequence_normalisee(float f, cstring msg)
  {
    string s = "F";
    si(!msg.empty())
      s = msg + " : f";
    si((f < 0) || (f > 0.5))
      échec("{}réquence normalisée attendue, f = {} (devrait être comprise entre 0 et 0.5).", s, f);
  }



float lexp_coef(Fréquence fc)
{
  retourne 1.0 - exp(-fc.value * 2 * π);
}

float lexp_tc_vers_coef(float τ)
{
  retourne lexp_coef(1.0 / (2 * π * τ));
}

Fréquence lexp_fcoupure(float γ)
{
  retourne -log(1.0 - γ) / (2 * π);
}

float lexp_coef_vers_tc(float γ)
{
  retourne 1.0 / (2 * π * lexp_fcoupure(γ).value);
}

FRat<float> design_notch(float f0, float fc)
{
  soit γ = lexp_coef(Fréquence(fc)),
       r = 1 - γ,
       Ω0 = 2 * π_f * f0;

  soit z  = FRat<float>::z();

  retourne r * ((z*z - 2 * cos(Ω0) * z + 1) / (z*z - 2 * r * cos(Ω0) * z + r*r));
}

FRat<float> design_bloqueur_dc(float fc)
{
  soit γ = lexp_coef(Fréquence(fc));
  soit r = 1 - γ;
  soit z = FRat<float>::z();
  retourne r * ((z - 1) / (z - r));
}

FRat<float> design_lexp(Fréquence fc)
{
  // yn = γ * xn + (1-γ) yn-1
  float γ = lexp_coef(fc);
  soit z = FRat<float>::z();
  retourne (z * γ) / (z  - (1 - γ));
}

FRat<float> design_lexp(float γ)
{
  retourne design_lexp(lexp_fcoupure(γ));
}



//using mapint = std::map<string, int>;

/*struct Exp
{
  template<typename T>
    operator T();

  template<typename T>
    eval()
  {
  }
};

Exp k;

Exp operator *(const Exp &x) const
{

}

template<typename T>
T somme(int imin, int imax, auto exp)
{
  T s = 0;
  for(auto i = imin; i <= imax; i++)
    s += exp.eval<T>({"k", i});
  retourne somme;
}*/


FRat<float> design_mg(int K)
{
  FRat<float> h;
  h.numer.coefs = Vecf::ones(K) / K;
  retourne h.horner(FRat<float>::zinv());

  /*retourne
    somme<FRat<float>>(0, K-1, pow(z, -Exp{"k"}));*/
}


}
