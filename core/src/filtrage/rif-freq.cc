#include "tsd/tsd-all.hpp"

namespace tsd::filtrage {


  Vecf design_rif_freq_freqs(entier n)
  {
    soit m = (n+1) / 2;
    retourne linspace(0,(m-1.0)/(2*m-1),m);
  }

  Figures design_rif_freq_analyse(entier n, const Vecf &d)
  {
    soit m = d.rows();
    assertion(2 * m - 1 == n);

    soit h = design_rif_freq(n, d);
    assertion(h.rows() == n);
    soit [fr,xm] = frmag(h, 2048);

    Figures fig;
    soit f = fig.subplot();
    f.plot(fr,abs(xm),"b-", "Réponse réelle");

    soit fr1 = design_rif_freq_freqs(n);
    assertion(fr1.rows() == m);
    f.plot(fr1, abs(d), "gs", "Points d'échantillonnage");
    f.titre("Réponse fréquentielle");

    retourne fig;
  }


  Vecf design_rif_freq(entier Nc, const Vecf &f)
  {
    soit F = f.clone();
    soit M = f.rows();

    si(Nc == 0)
      Nc = 2 * M - 1;

    msg("Design RIF FREQ: N (nombre de coefficients) = {}, M (nombre de points du gabarit) = {}", Nc, M);

    si(Nc != 2 * M - 1)
    {
      msg("design_rif_freq: interpolation gabarit {} points -> {} points.", M, (Nc+1)/2);
      si((Nc & 1) == 0)
        Nc++;
      // Adapte le nombre de coefficients
      M = (Nc + 1) / 2;
      // Réinterpolation
      F = interp(linspace(0,0.5*(1-1.0/f.rows()),f.rows()), f, linspace(0,0.5*(1-1.0/M),M));
    }

    assertion(Nc == 2 * M - 1);

    soit Hd = Veccf::zeros(Nc);
    // pour un filtre à phase linéaire
    //soit θ  = linspace(0, -(2 * π * M * (M-1)) / N, M);
    Hd.head(M).set_real(F);// * polar(θ);

    // M = (Nc+1)/2

    // 0               -> partie imaginaire = 0 (c déjà le cas)
    // 1..(Nc-1)/2     -> tel quel
    // (Nc+1)/2...Nc-1 -> conjugé sym
    csym_forçage(Hd);

    // Attention, problème de précision avec ifft !!!
    soit hc = fftshift(ifft(Hd));

    soit err = abs(imag(hc)).valeur_max();
    assertion_msg(err < 1e-3,
        "design_rif_freq : filtre réel attendu (err = {}): hc = \n{}\nHd={}, F={}, f={}.",
        err, hc, Hd, F, f);

    soit h = real(hc) / sqrt(Nc);

    si(debug_design)
    {
      msg("Design RIF freq : N={}, M={}", Nc, M);
      Figures fig;
      fig.subplot().plot(f, "", "Avant interpolation");
      fig.subplot().plot(F, "", "Après interpolation");
      fig.afficher("rif-freq-gabarit");


      soit [fr,xm] = frmag(h, 2048);

      soit f = fig.subplot();
      f.plot(fr,xm, "b-", "Réponse obtenue");

      soit fr1 = design_rif_freq_freqs(Nc);
      f.plot(fr1, F, "gs", "Points d'échantillonnage (gabarit)");
      f.titre("Réponse fréquentielle");

    }

    retourne h;
  }


}
