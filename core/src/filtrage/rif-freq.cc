#include "tsd/tsd-all.hpp"

namespace tsd::filtrage {


  Vecf design_rif_freq_freqs(entier n)
  {
    entier m = (n+1) / 2;
    retourne linspace(0,(m-1.0)/(2*m-1),m);
  }

  // TODEL
  Figures design_rif_freq_analyse(entier n, const Vecf &d)
  {
    soit h = design_rif_freq(n, d);
    soit [fr,xm] = frmag<float>(h, 2048);

    Figures fig;
    soit f = fig.subplot();
    f.plot(fr,abs(xm),"b-", "Réponse obtenue");

    soit fr1 = design_rif_freq_freqs(n);
    f.plot(fr1, abs(d), "gs", "Points d'échantillonnage (gabarit)");
    f.titre("Réponse fréquentielle");

    retourne fig;
  }


  Vecf design_rif_freq(entier N, const Vecf &f)
  {
    soit F = f.clone();
    soit M = f.rows();

    si(N == 0)
      N = 2 * M - 1;

    msg("Design RIF FREQ: N = {}, M = {}", N, M);

    si(N != 2 * M - 1)
    {
      msg("design_rif_freq: interpolation gabarit {} points -> {} points.", M, (N+1)/2);
      si((N & 1) == 0)
        N++;
      // Adapte le nombre de coefficients
      M = (N + 1) / 2;
      // Réinterpolation
      F = interp(linspace(0,0.5*(1-1.0/f.rows()),f.rows()), f, linspace(0,0.5*(1-1.0/M),M));
    }

    tsd_assert(N == 2 * M - 1);

    soit Hd = Veccf::zeros(N);
    // pour un filtre à phase linéaire
    //soit θ  = linspace(0, -(2 * π * M * (M-1)) / N, M);
    Hd.head(M).set_real(F);// * polar(θ);

    csym_forçage(Hd);

    // Attention, problème de précision avec ifft !!!
    soit hc = fftshift(ifft(Hd));

    soit err = abs(imag(hc)).valeur_max();
    tsd_assert_msg(err < 1e-3,
        "design_rif_freq : filtre réel attendu (err = {}): hc = \n{}\nHd={}, F={}, f={}.", err, hc, Hd, F, f);

    soit h = real(hc) / sqrt(N);

    si(debug_design)
    {
      msg("Design RIF freq : N={}, M={}", N, M);
      Figures fig;
      fig.subplot().plot(f, "", "Avant interpolation");
      fig.subplot().plot(F, "", "Après interpolation");
      fig.afficher("rif-freq-gabarit");


      soit [fr,xm] = frmag<float>(h, 2048);

      soit f = fig.subplot();
      f.plot(fr,xm, "b-", "Réponse obtenue");

      soit fr1 = design_rif_freq_freqs(N);
      f.plot(fr1, F, "gs", "Points d'échantillonnage (gabarit)");
      f.titre("Réponse fréquentielle");

    }

    retourne h;
  }


}
