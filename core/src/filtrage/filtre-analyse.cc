#include "tsd/tsd-all.hpp"

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


  entier type_rif(const Vecf &h)
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


  AnalyseFiltre analyse_filtre(const Design &d, bouléen avec_plots)
  {
    AnalyseFiltre res;
    soit npts     = 4096;
    soit [fr,mag] = frmag(d, npts);
    Vecf hr       = repimp(d);
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

      res.description += "Type de filtre inconnu (ni passe-bas, ni passe-haut, ni passe-bande, ni coupe-bande.\n";

      retourne res;
    }

    res.description += sformat("Type de filtre : {}.\n", res.type);

    si(d.est_rif)
    {
      res.type_rif = type_rif(hr);
      res.description += sformat("Type RIF: {} ({} coefficients).\n", res.type_rif, nc);
    }

    //res.gain_dc      = hr.somme();
    //res.gain_Nyquist = (signyquist(nc) * hr).somme();
    //msg("Gain dc = {}, Nyquist = {}", res.gain_dc, res.gain_Nyquist);

    soit h = d.complex_frat();

    // Directement d'après la fonction de transfert
    res.gain_dc       = abs(h.horner((cfloat) 1.0f));
    res.gain_Nyquist  = abs(h.horner((cfloat) -1.0f));
    res.description += sformat("Gain dc = {}, Nyquist = {}\n", res.gain_dc, res.gain_Nyquist);

    // magdB = énergie en log (carré de la magnitude)
    soit magdB = mag2db(mag);

    // 3 dB = 3 dB en énergie, soit énergie = 1/2, mag = 1/sqrt(2)
    // 6 dB = ...                             1/4, mag = 1/2

    soit index_3dB = -1, index_6dB = -1;

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


    si(avec_plots)
    {
      plot_plz(res.figures.plz, d, non);
      plot_plz(res.figures.plz_cmap, d, oui);
      plot_rimp(res.figures.repimp, d);
      plot_frmag(res.figures.repfreq_lin, d, non);
      plot_frmag(res.figures.repfreq_log, d, oui);
      plot_frphase(res.figures.repphase, d);
      plot_frdelay(res.figures.delay, d);
      plot_rech(res.figures.repech, d);
    }



   Figure f;

   si(avec_plots)
   {
     f = res.fig.subplot(121);
     f.plot(fr, mag, "b-");
     si(index_premier_lobe >= 0)
       f.plot(fr(index_premier_lobe), mag(index_premier_lobe), "rs");
     f.titre(est_fr() ? "Vue linéaire" : "Linear view");
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
         est_fr() ?
             sformat("Largeur lobe principal (-3 dB) : {:.5f} (={:.2f}/N)\nDébut bande atténuée : {:.3f}\nAttén. premier lobe sec. : {:.1f} dB\nAttén. pire lob sec. : {:.1f} dB",
             res.largeur_lp, res.largeur_lp * nc, res.freq_debut_ls, res.premier_ls.atten, res.pire_ls.atten) :
             sformat("Principal lob bandwidth (-3 dB): {:.5f} (={:.2f}/N)\nStart of stop-band: {:.3f}\nAtten. first secondary lobe : {:.1f} dB\nAtten. worst secondary lobe : {:.1f} dB",
                 res.largeur_lp, res.largeur_lp * nc, res.freq_debut_ls, res.premier_ls.atten, res.pire_ls.atten)
             ,
             {0.4f, 0.2f});

     f = res.fig.subplot(122);

     f.plot(fr, magdB, "b-");
     f.titre(est_fr() ? "Vue logarithmique" : "Logarithmic view");
     si(index_premier_lobe >= 0)
       f.plot(fr(index_premier_lobe), magdB(index_premier_lobe), "rs");
     si(index_pire_lobe >= 0)
       f.plot(fr(index_pire_lobe), magdB(index_pire_lobe), "ms");
   }

   res.description += sformat("  Atténuation lobe principal :          \033[33m{:.1f}\033[0m dB.\n", -magdB(0));
   res.description += sformat("  Largeur lobe principal (à -3 dB, demi-énergie)   :    \033[33m{:.4f}\033[0m (j={}, = {}/N)\n", res.largeur_lp, index_3dB, res.largeur_lp * nc);
   res.description += sformat("  Largeur lobe principal (à -6 dB, demi-amplitude) :    \033[33m{:.4f}\033[0m (j={}, = {}/N)\n", res.largeur_lp_6dB, index_6dB, res.largeur_lp_6dB * nc);
   si(index_premier_lobe >= 0)
   {
     res.description += sformat("  Début bande atténuée : \033[33m{:.3f}\033[0m\n", res.freq_debut_ls);
     res.description += sformat("  Atténuation premier lobe secondaire : \033[33m{:.1f}\033[0m dB (@f={}).\n",
         res.premier_ls.atten, res.premier_ls.freq);
     res.description += sformat("  Atténuation pire lobe secondaire :    \033[33m{:.1f}\033[0m dB (@f={}).\n",
         res.pire_ls.atten, res.pire_ls.freq);
   }
   sinon
     res.description += sformat("  Pas de lobe secondaire trouvé.\n");

    msg("{}", res.description);
    retourne res;
  }





  tuple<Vecf, Vecf> frphase(const Design &d, entier n)
  {
    soit fr  = linspace(0, 0.5 - (0.5 / n), n);
    soit dfr = polar(2 * π * fr);

    soit xm = Vecf::int_expr(n,
      IMAP(arg(d.complex_frat().horner(dfr(i)))));

    retourne {fr, déplie_phase(xm, 2 * π)};
  }


  tuple<Vecf, Vecf> frgroup(const Design &d, entier npts)
  {
    soit [fr, ϕ] = frphase(d, npts);
    soit [fr2, g] = frmag(d, npts);

    // Pas de phase pour la dérivée
    float delta_ϕ = 2 * π * (fr(1) - fr(0));

    soit  dϕ = diff(ϕ);

    Vecf ϕp(npts);
    ϕp.head(npts-1) = -dϕ / delta_ϕ;

    // En cas de phase discontinue, la dérivée n'est pas définie.
    pour(auto k = 0; k < npts-1; k++)
    {
      si(abs(dϕ(k)) > π_f/2)
        ϕp(k) = std::nanf("");
    }
    // Approximaton !
    ϕp(npts-1)      = ϕp(npts-2);

    retourne {fr, ϕp};
  }


  Veccf repfreq(const Design &d, const Vecf &fr)
  {
    retourne Veccf::int_expr(fr.rows(),
         IMAP(d.complex_frat().horner(std::polar(1.0f, 2 * π_f * fr(i)))));
  }


tuple<Vecf, Vecf> frmag_rif(const Vecf &h, entier npts)
{
  npts = max(npts, h.rows() / 2);
  soit h2 = Vecf::zeros(2*npts);
  h2.head(h.rows()) = h;

  retourne {
     linspace(0, 0.5 - (0.5 / npts), npts),
    (abs(fft(h2))).head(npts) * sqrt(2.0f*npts)};
}

tuple<Vecf, Vecf> frmag(const Design &d, entier npts)
{
  si(d.est_rif)
    retourne frmag_rif(d.coefs, npts);

  soit fr = linspace(0, 0.5 - (0.5 / npts), npts);
  retourne{fr, abs(repfreq(d, fr))};
}

Vecf repech(const Design &d, entier npts)
{
  si(npts == -1)
  {
    si(d.est_rif)
      npts = 2 * d.coefs.rows();
    sinon
      npts = 200;
  }

  retourne filtrer(d, Vecf::ones(npts));
}

Vecf repimp(const Design &d, entier npts)
{
  si(d.est_rif)
    retourne d.coefs;

  // Réponse infinie
  si(npts == -1)
    npts = (d.est_complexe ? d.frat_c.numer.coefs.rows() : d.frat.numer.coefs.rows()) * 20;

  soit x = sigimp(npts),
       r = filtrer(d, x);

  tantque((r.rows() > 1) && (abs(r(r.rows()-1)) < abs(r).moyenne() * 0.01))
    r = r.head(r.rows()-1).clone();

  retourne r;
}

}
