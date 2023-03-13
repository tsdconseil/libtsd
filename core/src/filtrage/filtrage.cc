#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <set>

using namespace std;
using namespace tsd::fourier;

namespace tsd::filtrage
{

  bool debug_design = false;


  Vecf design_rif_prod(const Vecf &h1, const Vecf &h2)
  {
    soit n1  = h1.rows();
    soit h2p = vconcat(h2, Vecf::zeros(n1-1));
    retourne filtrer(h1, h2p);
  }

  Vecf design_rif_pb2ph_rs(const Vecf &h)
  {
    // Attention, nécessite que h soit de type I
    tsd_assert_msg(type_rif(h) == 1, "design_rif_pb2ph_rs: h doit être un filtre de type I.");
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
    tsd_assert_msg(type_rif(h) == 1, "design_rif_pb2ph_is: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;
  }

  Vecf design_rif_pb2pb(const Vecf &h)
  {
    /*tsd_assert_msg(type_rif(h) == 1, "design_rif_pb2ph: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;*/
    retourne h;
  }

  Vecf design_rif_pb2pm(const Vecf &h)
  {

    /*tsd_assert_msg(type_rif(h) == 1, "design_rif_pb2ph: h doit être un filtre de type I.");
    soit n = h.rows();
    soit g = -h;
    g((n-1)/2)++;
    retourne g;*/
    retourne h;
  }


  void verifie_frequence_normalisee(float f, const string &msg)
  {
    string s = "F";
    si(!msg.empty())
      s = msg + " : f";
    si((f < 0) || (f > 0.5))
      echec("{}réquence normalisée attendue, f = {} (devrait être comprise entre 0 et 0.5).", s, f);
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

float bloqueur_dc_coef(float fc)
{
  soit c = cos(2 * π * fc);
  retourne (1 - sqrt(1 - c * c)) / c;
}

FRat<float> design_notch(float f0, float fc)
{
  soit a  = bloqueur_dc_coef(fc);
  soit z  = FRat<float>::z();
  soit Ω0 = 2 * π_f * f0;
  soit h = ((z*z - 2 * cos(Ω0) * z + 1) / (z*z - 2 * a * cos(Ω0) * z + a*a));
  retourne h * (1 / h.horner(1.0f)); // Equilibre le gain dans la bande passante (du moins est exact en DC)
}

FRat<float> design_bloqueur_dc(float fc)
{
  soit a = bloqueur_dc_coef(fc);
  soit z = FRat<float>::z();
  retourne ((z - 1) / (z - a)) * ((1 + a) / 2);
  //soit zi = (FRat<float>::z()).inv();
  //retourne ((1 - zi) / (1 - a * zi)) * ((1 + a) / 2);
}

FRat<float> design_lexp(Fréquence fc)
{
  // yn = γ * xn + (1-γ) yn-1
  // h = γ / (1-(1-γ)z^-1)
  //   = γz / (z -(1-γ)
  float γ = lexp_coef(fc);
  soit z = FRat<float>::z();
  retourne (z * γ) / (z  - (1 - γ));
}

FRat<float> design_lexp(float γ)
{
  retourne design_lexp(lexp_fcoupure(γ));
}



}
