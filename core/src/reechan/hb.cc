#include "tsd/tsd-all.hpp"

namespace tsd::filtrage
{

  Vecf design_rif_demi_bande(int n, float fc)
  {
    ////////////////////////////////////////////////////////
    // D'après :                                          //
    // A “TRICK” for the Design of FIR Half-Band Filters  //
    // P. P. VAIDYANATHAN AND TRUONG Q. NGUYEN, 1987      //
    ////////////////////////////////////////////////////////

    // Attention, on a forcément n impair
    // Et idéalement (n-1)/2 impair, sinon ça ajoute deux zéros inutiles en bout de filtre.

    // Deux cas possibles :
    //  - m pair
    //    exemple : n = 5, m = 2
    //      h = 0 + h(-2) z^-1 + 1/2 + h(2) z^1 + 0 [coefs d'ordre pair sont nul, sauf h_m]
    //  - m impair
    //    exemple : n = 3, m = 1
    //      h = h(-2) z^-1 + 1/2 + h(2) z^1 [coefs d'ordre impair sont nul, sauf h_m]

    tsd_assert_msg(est_impair(n), "design_rif_hb(n={}) : filtre de type I, n doit être impair.", n);

    soit m = (n-1) / 2;

    soit m_impair = est_impair(m);

    si(!m_impair)
    {
      msg_avert("design_rif_hb(n={}): m=(n-1)/2={}. m n'est pas impair => deux zéros inutiles en bout de filtre.", n, m);
      // Se ramène à m impair en ajoutant deux zéros à la fin
      n -= 2;
      m--;
    }

    soit L = 16 * m;
    soit D = Vecf::ones(L),
         W = Vecf::ones(L);

    // fp=0.5 <=> L

    soit idp = (int) ceil((2 * fc * L) / 0.5);
    idp = clamp(idp, 0, L-2);

    msg("design_hb: n={}, m={}, L={}, index pass={}", n, m, L, idp);
    msg("fc={} => 2 * fc = {}", fc, 2 * fc);

    W.segment(idp, L-idp).setZero();

    // Non : D(L-1), ce n'est pas le Nyquist, c'est un peu avant...
    //D(L-1) = 0;
    //W(L-1) = 1;
    // Mais pas besoin, car le fait que le nombre de coefficients soit pair,
    // et que le filtre soit symétrique, implique nécessairement que la réponse
    // est nulle à la fréquence de Nyquist.


    soit g = design_rif_eq(m + 1, D, W);

    tsd_assert(g.rows() == m + 1);




    //msg("g: {}", g);

    Vecf h(n);
    pour(auto k = 0; k < n; k++)
    {
      si(k == m)
        h(k) = 0.5;
      sinon si(est_impair(k))
        h(k) = 0;
      sinon
        h(k) = 0.5 * g(k/2);
    }

    si(!m_impair)
      h = vconcat(Vecf::zeros(1), vconcat(h, Vecf::zeros(1)));


    si(debug_design)
    {
      plot_filtre(g).afficher("Filtre g");
      plot_filtre(h).afficher("Filtre h");
    }

    //msg("h: {}", h);
    retourne h;
  }





}
