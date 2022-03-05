#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include "tsd/divers.hpp"
#include <cmath>

namespace tsd::filtrage {






static ArrayXf coefs_filtre_sinc(int n, float fc)
{
  ArrayXf h(n);

  // n impair
  if((n & 1) == 1)
  {
    for(auto i = 0; i < n; i++)
      h(i) = sinc(2*fc, i - n/2);
  }
  else
  {
    // n = 4
    // -1..2
    for(auto i = 0; i < n; i++)
      h(i) = sinc(2*fc, i - (n-1)/2);
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

