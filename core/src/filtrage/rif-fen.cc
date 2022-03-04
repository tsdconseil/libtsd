#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include "tsd/divers.hpp"
#include <cmath>

namespace tsd::filtrage {


// TODO : à clarifier !!!
float sinc2(float t, float fc)
{
  if(std::abs(t) < 1e-7f)
    return 2 * fc;
  return sin(2 * π * t * fc) / (π * t);

  //   return std::sin(T * omega / 2) / (omega / 2);
  // T*o/2 = 2pi*t*fc
  // o/2 = pi*t

  // <=> o = 2*pi*t
  //     T = 2 * fc

  // TFI d'une fonction porte fréquentielle, de la largeur 2fc
  //return std::sqrt(2*π_f) * tsd::sinc(2 * fc, 2 * π * t);
}

// TODO : à clarifier !!!
float sinc(float t)
{
  // sin(pi t) / pi t

  if(std::abs(t) < 1e-7f)
    return 1;
  return std::sin(π_f*t) / (π_f*t);

  //return std::sqrt(2*π_f) * tsd::sinc(1, 2 * π * t);
  //return sinc2(t, 0.5f);
}

static ArrayXf coefs_filtre_sinc(int n, float fc)
{
  ArrayXf h(n);

  // n impair
  if((n & 1) == 1)
  {
    for(auto i = 0; i < n; i++)
      h(i) = sinc2(i - n/2, fc);
  }
  else
  {
    // n = 4
    // -1..2
    for(auto i = 0; i < n; i++)
      h(i) = sinc2(i - (n-1)/2, fc);
  }
  return h;
}


static ArrayXf rif_fen_hp(int n, float fc)
{
  int no2  = (n-1)/2;
  ArrayXf h = -coefs_filtre_sinc(n, fc);
  h(no2)++;
  return h;
}



static ArrayXf rif_fen_lp(int n, float fc)
{
  // Fenêtrage d'un filtre idéal de type porte fréquentielle
  return coefs_filtre_sinc(n, fc);
}

static ArrayXf rif_fen_bp(int n, float fl, float fh)
{
  int no2  = (n-1)/2;
  auto ωc = π * (fh + fl);
  auto δf = (fh - fl) / 2;
  ArrayXf h = coefs_filtre_sinc(n, δf);
  return 2 * h * cos(ωc * intervalle_entier(-no2,no2).cast<float>());
}

static ArrayXf rif_fen_sb(int n, float fl, float fh)
{
  auto no2  = (n-1)/2;
  auto ωc = π * (fh + fl);
  auto δf = (fh - fl) / 2;
  ArrayXf h = -coefs_filtre_sinc(n,δf);
  h *= 2 * cos(ωc*intervalle_entier(-no2,no2).cast<float>());
  h(no2)++;
  return h;
}


static ArrayXf design_rif_fen(const std::string &type, float fcut, const ArrayXf &f, float fcut2)
{
  ArrayXf h;
  int n = f.rows();
  if(type == "lp")
    h = rif_fen_lp(n, fcut);
  else if(type == "hp")
    h = rif_fen_hp(n, fcut);
  else if(type == "bp")
    h = rif_fen_bp(n, fcut, fcut2);
  else if(type == "sb")
    h = rif_fen_sb(n, fcut, fcut2);
  else
    echec("design_rif_fen : type invalide ({}). doit être parmi : \"lp\", \"hp\", \"bp\", \"sb\"", type);
  tsd_assert(h.rows() == n);
  ArrayXf h2 = h * f;
  if(type == "lp")
    h2 /= h2.sum();
  return h2;
  // TODO : normalisation gain pour hp, bp, sb
}


ArrayXf design_rif_fen(unsigned int n, const std::string &type, float fcut, const std::string &fen, float fcut2)
{
  ArrayXf f = fenetre(fen, n, true);
  return design_rif_fen(type, fcut, f, fcut2);
}

ArrayXf design_rif_fen_kaiser(const std::string &type, float fcut, float atten_db, float df, float fcut2)
{
  ArrayXf f = fenêtre_kaiser(atten_db, df);
  return design_rif_fen(type, fcut, f, fcut2);
}

ArrayXf design_rif_fen_chebychev(int n, const std::string &type, float fcut, float atten_db, float fcut2)
{
  ArrayXf f = fenêtre_chebychev(n, atten_db, true);
  return design_rif_fen(type, fcut, f, fcut2);
}


}

