#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include "tsd/divers.hpp"
#include <cmath>

namespace tsd::filtrage {





  // sin(pi * 2 * fc * (i-n/2))
  // --------------------------
  //      pi * (i-n/2)


  //  sin(2 * pi * fc * k)
  //  -------------------
  //       pi * k

  // fc = 0.25:
  //  sin(pi * k * 0.5)
  //  -------------------
  //       pi * k



static Vecf coefs_filtre_sinc(entier n, float fc)
{
  si(est_impair(n))
    retourne Vecf::int_expr(n, IMAP(
        sinc(2 * fc, i - n/2)));

  // n = 4
  // -1..2
  retourne Vecf::int_expr(n, IMAP(
      sinc(2 * fc, i - (n-1)/2)));
}


static Vecf rif_fen_hp(entier n, float fc)
{
  entier no2  = (n-1)/2;
  soit h = -coefs_filtre_sinc(n, fc);
  h(no2)++;
  retourne h;
}



static Vecf rif_fen_pb(entier n, float fc)
{
  // Fenêtrage d'un filtre idéal de type porte fréquentielle
  retourne coefs_filtre_sinc(n, fc);
}

static Vecf rif_fen_pm(entier n, float fl, float fh)
{
  soit no2  = (n-1)/2;
  soit ωc = π * (fh + fl);
  soit δf = (fh - fl) / 2;
  soit h = coefs_filtre_sinc(n, δf);
  retourne 2 * h * cos(ωc * intervalle_entier(-no2,no2).as<float>());
}

static Vecf rif_fen_sb(entier n, float fl, float fh)
{
  soit no2  = (n-1)/2;
  soit ωc = π * (fh + fl);
  soit δf = (fh - fl) / 2;
  soit h = -coefs_filtre_sinc(n,δf);
  h *= 2 * cos(ωc*intervalle_entier(-no2,no2).as<float>());
  h(no2)++;
  retourne h;
}


static Vecf design_rif_fen(cstring type, float fcut, const Vecf &f, float fcut2)
{
  Vecf h;
  soit n = f.rows();
  si((type == "lp") || (type == "pb"))
    h = rif_fen_pb(n, fcut);
  sinon si((type == "hp") || (type == "ph"))
    h = rif_fen_hp(n, fcut);
  sinon si((type == "bp") || (type == "pm"))
    h = rif_fen_pm(n, fcut, fcut2);
  sinon si(type == "sb")
    h = rif_fen_sb(n, fcut, fcut2);
  sinon
    échec("design_rif_fen : type invalide ({}). doit être parmi : \"lp\" ou \"pb\" (passe-bas), \"hp\" ou \"ph\" (passe-haut), \"bp\" ou \"pm\" (passe-bande), \"sb\" (stoppe-bande)", type);
  assertion(h.rows() == n);
  soit h2 = h * f;
  si(type == "lp")
    h2 /= h2.somme();
  retourne h2;
  // TODO : normalisation gain pour hp, bp, sb
}


Vecf design_rif_fen(entier n, const std::string &type, float fcut, cstring fen, float fcut2)
{
  soit f = fenetre(fen, n, oui);
  retourne design_rif_fen(type, fcut, f, fcut2);
}

Vecf design_rif_fen_kaiser(cstring type, float fcut, float atten_db, float df, float fcut2)
{
  soit f = fenêtre_kaiser(atten_db, df);
  retourne design_rif_fen(type, fcut, f, fcut2);
}

Vecf design_rif_fen_chebychev(entier n, cstring type, float fcut, float atten_db, float fcut2)
{
  soit f = fenêtre_chebychev(n, atten_db, oui);
  retourne design_rif_fen(type, fcut, f, fcut2);
}


}

