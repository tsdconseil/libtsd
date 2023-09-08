#include "tsd/tsd-all.hpp"


namespace tsd::filtrage {



Vecf design_rif_hilbert(entier n, cstring fen)
{
  // Attention, filtre de type IV (n impair)
  assertion_msg((n & 1) == 1, "design_rif_hilbert: n doit être impair (ici, n = {})", n);

  soit h = Vecf::int_expr(n, [&](int i) -> float
  {
    soit k = i - n / 2;

    retourne est_pair(k) ? 0.0 : 2 / (k * π);
  });

  retourne h * fenêtre(fen, n);
}





struct HilbertTransformeur: Filtre<float, cfloat, HilbertTransformeurConfig>
{
  sptr<FiltreGen<float>> retard, hilbert;

  HilbertTransformeur(entier ntaps, cstring fenetre)
  {
    configure({ntaps, fenetre});
  }
  void configure_impl(const HilbertTransformeurConfig &config)
  {
    soit h = design_rif_hilbert(config.ntaps, config.fenetre);
    hilbert   = filtre_rif<float>(h);
    retard    = ligne_a_retard<float>(config.ntaps/2);
  }
  void step(const Vecf &x, Veccf &y)
  {
    soit x2 = hilbert->step(x),
         // Ligne à retard sur y, pour être synchrone avec le signal filtré
         y2 = retard->step(x);
    y.resize(x.rows());
    y.set_real(y2);
    y.set_imag(x2);
  }
};


sptr<Filtre<float, cfloat, HilbertTransformeurConfig>> hilbert_transformeur(entier ntaps, cstring fenetre)
{
  retourne std::make_shared<HilbertTransformeur>(ntaps, fenetre);
}

Veccf hilbert(const Vecf &x, entier ncoefs)
{
  retourne hilbert_transformeur(ncoefs)->step(x);
}

Veccf hilbert_tfd(const Vecf &x)
{
  // TODO : vérifier n pair et impair
  soit n = x.rows();
  soit X = fft(x);
  X.tail(n/2).setZero();
  retourne ifft(2 * X);
}

}





