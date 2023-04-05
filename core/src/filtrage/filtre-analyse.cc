#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include "tsd/fourier.hpp"
#include <iostream>

using namespace tsd::vue;
using namespace tsd::fourier;
using namespace std;

namespace tsd::filtrage {


  std::ostream& operator<<(std::ostream& ss, const TypeFiltre &t)
  {
    switch(t)
    {
    case TypeFiltre::PASSE_BAS:
      ss << "passe-bas";
      break;
    case TypeFiltre::PASSE_HAUT:
      ss << "passe-haut";
      break;
    case TypeFiltre::PASSE_BANDE:
      ss << "passe-bande";
      break;
    case TypeFiltre::COUPE_BANDE:
      ss << "coupe-bande";
      break;
    case TypeFiltre::AUTRE:
      ss << "autre";
      break;
    default:
      ss << "?";
    }
    retourne ss;
  }


  int type_rif(const Vecf &h)
  {
    soit nc = h.rows();
    Vecf hrr = h.reverse();
    si(hrr.est_approx(h))
      retourne est_impair(nc) ? 1 : 2;
    sinon si(hrr.est_approx(-h))
      retourne est_impair(nc) ? 3 : 4;
    msg_avert("type_rif(): type inconnu (ni symétrique ni anti-symétrique).");
    retourne -1;
  }


  template<typename T>
  AnalyseLIT analyse_LIT(const FRat<T> &h, bool avec_plots)
  {
    AnalyseLIT res;
    soit npts = 4096;//2048;
    soit [fr,mag] = frmag(h, npts);
    Vecf hr       = repimp(h);
    soit nc       = hr.rows();

    // (1) Détermine le type de filtre (passe-bas, passe-haut, etc.)

    soit m0 = mag(0), m1 = mag(npts-1);

    soit SEUIL = 0.2;

    si((m0 > 1 - SEUIL) && (m1 < SEUIL))
      res.type = TypeFiltre::PASSE_BAS;
    sinon si((m1 > 1 - SEUIL) && (m0 < SEUIL))
      res.type = TypeFiltre::PASSE_HAUT;
    sinon si((m1 > 1 - SEUIL) && (m1 > 1 - SEUIL))
      res.type = TypeFiltre::COUPE_BANDE;
    sinon si((m1 < SEUIL) && (m1 < SEUIL))
      res.type = TypeFiltre::PASSE_BANDE;
    sinon si((m0 > 1 - SEUIL) && (m1 < m0))
      res.type = TypeFiltre::PASSE_BAS;
    sinon si((m1 > 1 - SEUIL) && (m0 < m1))
      res.type = TypeFiltre::PASSE_HAUT;
    sinon
    {
      msg_avert("analyse_LIT(): type de filtre non géré (ni passe-bas, ni passe-haut, ni passe-bande, ni coupe-bande).");
      msg("mag(0) = {}, mag($) = {}", m0, m1);
      retourne res;
    }

    msg("type de filtre : {}", res.type);

    si(h.est_rif())
    {
      res.type_rif = type_rif(hr);
      msg("Type RIF: {} ({} coefficients)", res.type_rif, nc);
    }

    //res.gain_dc      = hr.somme();
    //res.gain_Nyquist = (signyquist(nc) * hr).somme();
    //msg("Gain dc = {}, Nyquist = {}", res.gain_dc, res.gain_Nyquist);

    // Directement d'après la fonction de transfert
    res.gain_dc       = abs(h.horner((T) 1.0f));
    res.gain_Nyquist  = abs(h.horner((T) -1.0f));
    msg("Gain dc = {}, Nyquist = {}", res.gain_dc, res.gain_Nyquist);

    // magdB = énergie en log (carré de la magnitude)
    soit magdB = mag2db(mag);

    // 3 dB = 3 dB en énergie, soit énergie = 1/2, mag = 1/sqrt(2)
    // 6 dB = ...                             1/4, mag = 1/2

    entier index_3dB = -1, index_6dB = -1;

    si(res.type == TypeFiltre::PASSE_BAS)
    {
      index_3dB = trouve_premier(magdB < -3);
      index_6dB = trouve_premier(magdB < -6);
    }
    sinon si(res.type == TypeFiltre::PASSE_HAUT)
    {
      index_3dB = trouve_dernier(magdB >= -3);
      index_6dB = trouve_dernier(magdB >= -6);
    }


    soit index_premier_lobe = -1;

    // Recherche premier lobe secondaire
    si(index_3dB >= 0)
    {
      index_premier_lobe = trouve_premier_max_local(magdB.tail(npts - index_3dB));
      si(index_premier_lobe >= 0)
        index_premier_lobe += index_3dB;
      sinon
      {
        // Sinon prends la dernière valeur
        index_premier_lobe = npts-1;
      }
    }

   Figure f;

   si(avec_plots)
   {
     f = res.fig.subplot(121);
     f.plot(fr, mag, "b-");
     si(index_premier_lobe >= 0)
       f.plot(fr(index_premier_lobe), mag(index_premier_lobe), "rs");
     f.titre("Vue linéaire");
   }

   soit idx2freq = [&](entier idx) -> float
   {
     si(idx < 0)
       retourne -1;
     retourne idx / (2.0 * npts);
   };

   res.largeur_lp     = idx2freq(index_3dB);
   res.largeur_lp_6dB = idx2freq(index_6dB);

   soit index_pire_lobe = index_premier_lobe;

   si(index_premier_lobe >= 0)
   {
     res.pire_ls.atten  = -magdB.tail(npts - index_premier_lobe).maxCoeff(&index_pire_lobe);
     index_pire_lobe += index_premier_lobe;
     res.premier_ls.atten = -magdB(index_premier_lobe);
   }

   res.premier_ls.freq  = idx2freq(index_premier_lobe);
   res.pire_ls.freq     = idx2freq(index_pire_lobe);

   soit idx_ls = trouve_premier(magdB < -res.pire_ls.atten);
   res.freq_debut_ls = idx2freq(idx_ls);


   si(avec_plots)
   {
     si(index_3dB >= 0)
       f.plot(fr(index_3dB), mag(index_3dB), "gs");

     f.canva().set_couleur(Couleur::Bleu);
     f.canva().set_dim_fonte(0.6);
     f.canva().texte({0.1f, 0.9f},
             fmt::format("Largeur lobe principal (-3 dB) : {:.5f} (={:.2f}/N)\nDébut bande atténuée : {:.3f}\nAttén. premier lobe sec. : {:.1f} dB\nAttén. pire lob sec. : {:.1f} dB",
             res.largeur_lp, res.largeur_lp * nc, res.freq_debut_ls, res.premier_ls.atten, res.pire_ls.atten),
             {0.4f, 0.2f});

     f = res.fig.subplot(122);

     f.plot(fr, magdB, "b-");
     f.titre("Vue logarithmique");
     si(index_premier_lobe >= 0)
       f.plot(fr(index_premier_lobe), magdB(index_premier_lobe), "rs");
     si(index_pire_lobe >= 0)
       f.plot(fr(index_pire_lobe), magdB(index_pire_lobe), "ms");
   }

   msg("  Atténuation lobe principal :          \033[33m{:.1f}\033[0m dB.", -magdB(0));
   msg("  Largeur lobe principal (à -3 dB, demi-énergie)   :    \033[33m{:.4f}\033[0m (j={}, = {}/N)", res.largeur_lp, index_3dB, res.largeur_lp * nc);
   msg("  Largeur lobe principal (à -6 dB, demi-amplitude) :    \033[33m{:.4f}\033[0m (j={}, = {}/N)", res.largeur_lp_6dB, index_6dB, res.largeur_lp_6dB * nc);
   si(index_premier_lobe >= 0)
   {
     msg("  Début bande atténuée : \033[33m{:.3f}\033[0m", res.freq_debut_ls);
     msg("  Atténuation premier lobe secondaire : \033[33m{:.1f}\033[0m dB (@f={}).",
         res.premier_ls.atten, res.premier_ls.freq);
     msg("  Atténuation pire lobe secondaire :    \033[33m{:.1f}\033[0m dB (@f={}).",
         res.pire_ls.atten, res.pire_ls.freq);
   }
   sinon
     msg("  Pas de lobe secondaire trouvé.");

    retourne res;
  }




  template<typename T>
  tuple<Vecf, Vecf> frphase(const FRat<T> &h, entier n)
  {
    soit fr  = linspace(0, 0.5 - (0.5 / n), n);
    soit dfr = polar(2 * π * fr);

    soit xm = Vecf::int_expr(n,
      IMAP(arg(h.horner(dfr(i)))));

    retourne {fr, déplie_phase(xm, 2 * π)};
  }



  template<typename T> tuple<Vecf, Vecf> frgroup(const Vecteur<T> &h, entier npts)
  {
    retourne frgroup(FRat<T>::rif(h), npts);
  }

  template<typename T> tuple<Vecf, Vecf> frmag(const Vecteur<T> &h, entier npts)
  {
    npts = max(npts, h.rows() / 2);
    soit h2 = Vecteur<T>::zeros(2*npts);
    h2.head(h.rows()) = h;

    retourne {
       linspace(0, 0.5 - (0.5 / npts), npts),
      (abs(fft(h2))).head(npts) * sqrt(2.0f*npts)};
  }

  template<typename T> tuple<Vecf, Vecf> frphase(const Vecteur<T> &h, entier npts)
  {
    retourne frphase(FRat<T>::rif(h), npts);
  }

  template<typename T>
  tuple<Vecf, Vecf> frgroup(const FRat<T> &h, entier npts)
  {
    soit [fr, ϕ] = frphase(h, npts);

    // Pas de phase pour la dérivée
    float delta_ϕ = ϕ(1) - ϕ(0);

    Vecf ϕp(npts);
    ϕp.head(npts-1) = diff(ϕ) / delta_ϕ;
    ϕp(npts-1)      = ϕp(npts-2);

    retourne {fr, ϕp};
  }


  template<typename T>
  Veccf repfreq(const FRat<T> &h, const Vecf &fr)
  {
    retourne Veccf::int_expr(fr.rows(),
         IMAP(h.horner(std::polar(1.0f, 2 * π_f * fr(i)))));
  }

  Veccf repfreq(const Vecf &h, const Vecf &fr)
  {
    retourne repfreq(FRat<float>::rif(h), fr);
  }

  template<typename T>
  tuple<Vecf, Vecf> frmag(const FRat<T> &h, entier npts)
  {
    soit fr = linspace(0, 0.5 - (0.5 / npts), npts);
    retourne{fr, abs(repfreq(h, fr))};
  }



template<typename T>
Vecf repimp(const FRat<T> &h, entier npts)
{
  si(h.est_rif())
    retourne h.coefs_rif();

  // Réponse infinie
  si(npts == -1)
    npts = h.numer.coefs.rows() * 20;

  // si T = complexe, la réponse est forcément complexe !!!

# if 0
  si constexpr (est_complexe<T>())
  {
    soit x = sigimp(npts);
    soit y = filtrer(h, x);
    tantque((y.rows() > 1) && (abs(y(y.rows()-1)) < abs(y).moyenne() * 0.01))
      y = y.head(y.rows()-1).clone();
    retourne real(y);

  }
  sinon
# endif
  {
    soit x = sigimp(npts);
    soit y = filtrer(h, x);
    tantque((y.rows() > 1) && (abs(y(y.rows()-1)) < abs(y).moyenne() * 0.01))
      y = y.head(y.rows()-1).clone();
    retourne real(y);
  }
}

template<typename T> Vecf repimp(const Vecteur<T> &h, entier npts)
{
  retourne repimp(FRat<T>::rif(h), npts);
}





template Vecf repimp(const Vecteur<float> &h, entier npts);
template Vecf repimp(const Vecteur<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frgroup<float>(const FRat<float> &h, entier npts);
template tuple<Vecf, Vecf> frgroup<float>(const Vecteur<float> &h, entier npts);
template tuple<Vecf, Vecf> frmag<float>(const FRat<float> &h, entier npts);
template tuple<Vecf, Vecf> frmag<float>(const Vecteur<float> &h, entier npts);
template tuple<Vecf, Vecf> frphase<float>(const FRat<float> &h, entier npts);
template tuple<Vecf, Vecf> frphase<float>(const Vecteur<float> &h, entier npts);
template tuple<Vecf, Vecf> frgroup<cfloat>(const FRat<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frgroup<cfloat>(const Vecteur<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frmag<cfloat>(const FRat<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frmag<cfloat>(const Vecteur<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frphase<cfloat>(const FRat<cfloat> &h, entier npts);
template tuple<Vecf, Vecf> frphase<cfloat>(const Vecteur<cfloat> &h, entier npts);
template AnalyseLIT analyse_LIT(const FRat<float> &h, bool avec_plots);
template AnalyseLIT analyse_LIT(const FRat<cfloat> &h, bool avec_plots);

}
