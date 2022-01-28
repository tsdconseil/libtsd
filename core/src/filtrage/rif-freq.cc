#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/figure.hpp"

using namespace tsd::vue;

namespace tsd::filtrage {


  ArrayXf design_rif_freq_freqs(int n)
  {
    int m = (n+1) / 2;
    return linspace(0,(m-1.0)/(2*m-1),m);
  }

  Figures design_rif_freq_analyse(int n, const ArrayXf &d)
  {
    ArrayXf h = design_rif_freq(n, d);

    auto [fr,xm] = frmag<float>(h, 2048);
    Figures fig;
    auto f = fig.subplot();
    f.plot(fr,abs(xm),"b-", "Réponse obtenue");
    ArrayXf fr1 = design_rif_freq_freqs(n);
    f.plot(fr1, d.abs(), "gs", "Points d'échantillonnage (gabarit)");
    f.titre("Réponse fréquentielle");

    return fig;
  }


  ArrayXf design_rif_freq(int N, const ArrayXf &f)
  {
    ArrayXf F = f;
    int M = f.rows();
    if(N == 0)
      N = 2 * M - 1;

    msg("Design RIF FREQ: N = {}, M = {}", N, M);

    if(N != 2 * M - 1)
    {
      msg("design_rif_freq: interpolation gabarit {} points -> {} points.", M, (N+1)/2);
      if((N & 1) == 0)
        N++;
      // Adapte le nombre de coefficients
      M = (N + 1) / 2;
      // Réinterpolation
      F = interp(linspace(0,0.5*(1-1.0/f.rows()),f.rows()), f, linspace(0,0.5*(1-1.0/M),M));
    }

    tsd_assert(N == 2 * M - 1);

    /*if(debug)
    {
      msg("Design RIF freq : ntaps = {}", ncoefs);
      Figures fig;
      fig.subplot().plot(f.abs(), "", "Avant interpolation");
      fig.subplot().plot(F.abs(), "", "Après interpolation");
      fig.afficher("rif-freq-gabarit");
    }*/


    ArrayXcf Hd = ArrayXcf::Zero(N);
    // Pour un filtre à phase linéaire
    //ArrayXf θ  = linspace(0, -(2 * π * M * (M-1)) / N, M);
    Hd.head(M) = F;// * polar(θ);

    tsd::fourier::force_csym(Hd);

    // Attention, problème de précision avec ifft !!!
    ArrayXcf hc = tsd::fourier::fftshift(tsd::fourier::ifft(Hd));
    //ArrayXcf hc = tsd::fourier::ifft(Hd);

    auto err = hc.imag().abs().maxCoeff();
    tsd_assert_msg(err < 1e-3, "design_rif_freq : filtre réel attendu (err = {}).", err);

    ArrayXf h = hc.real();

    // Exemple : gabarit = constant = 1
    // => fft(1) = 0.racine(N).000000
    h /= std::sqrt(N);

    return h;
  }


}
