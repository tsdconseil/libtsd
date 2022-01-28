#include "tsd/filtrage.hpp"
#include "tsd/filtrage/spline.hpp"

#define VERB(AA)
//AA
#define RA_SAFE_MODE 0
//1

namespace tsd::filtrage {

/** Maximum decimation ration = 1/2^16 */
const auto RATE_ADAPTATION_NB_HALF_MAX = 16;



template<typename T>
struct AdaptationRythmeSimple: FiltreGen<T>
{
  float phase;
  float ratio;
  /** Inverse du ratio precedent */
  float increment;
  Vecteur<T> fenetre;
  int nfen = 0;
  sptr<Interpolateur<T>> interpolateur;

  AdaptationRythmeSimple(float ratio, sptr<Interpolateur<T>> itrp)
  {
    this->ratio   = ratio;
    increment     = 1.0f / ratio;
    phase         = 0;
    interpolateur = itrp;
    nfen          = itrp->npts;
    fenetre.setZero(nfen);
    if(nfen <= 1)
      msg_avert("AdaptationRythmeSimple : npts interpolation = {} ({}).", nfen, itrp->nom);
  }

  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    int n = x.rows();

    if(n == 0)
    {
      y.resize(0);
      return;
    }

    int tmp_len = ceil(ratio * n) + 10;

    Vecteur<T> temp(tmp_len);

    auto j = 0;
    auto optr = temp.data();
    auto iptr = x.data();

    for(auto i = 0; i < n; i++)
    {
      // Rotation de la fenêtre
      // TODO : ça serait pas mal que les interpolateurs puissent fonctionner de manière circulaire
      // pour éviter cette rotation
      fenetre.head(nfen-1) = fenetre.tail(nfen-1).eval();
      fenetre(nfen-1) = *iptr++;

      while(phase < 1)
      {
        // phase = index entre deux échantillons
        auto si = interpolateur->calcule(fenetre, 0, phase);
        tsd_assert(j < tmp_len);
        *optr++ = si;
        j++;
        phase += increment;
      }
      phase--;
    } // end for i
    y = temp.head(j);
  }

};



/** Adapt the sampling rate, arbitrary ratio */
template<typename T>
struct AdaptationRythmeArbitraire: FiltreGen<T>
{
  sptr<FiltreGen<T>> interpolateur;
  std::vector<sptr<FiltreGen<T>>> décimateurs, suréchantilloneurs;

  /** Ratio d'interpolation / post-décimation */
  float facteur_post_interpolation;


  /** Decimation / interpolation ratio.
      e.g. d < 1.0 for decimation, d > 1.0 for interpolation */
  float ratio;


  AdaptationRythmeArbitraire(float ratio)
  {
    configure(ratio);
  }

  /** Setup the rate adaptation */
  int configure(const float &ratio_)
  {
    ratio = ratio_;

    if((ratio_ <= 0) || std::isinf(ratio_) || (ratio_ >= 1e9))
    {
      msg_erreur("Facteur de décimation invalide : {}.", ratio_);
      ratio = 1;
    }

    VERB(msg("ratio de décimation / interpolation demandé : {}.", ratio_);)

    float ratio2 = ratio;
    auto nb_suréchantilloneurs = 0u, nb_décimateurs = 0u;

    while(ratio2 < 0.5)
    {
      nb_décimateurs++;
      ratio2 *= 2.0;
    }

    nb_suréchantilloneurs = 0;
    while(ratio2 >= 2)
    {
      nb_suréchantilloneurs++;
      ratio2 /= 2;
    }

    // Se place dans l'intervalle [0,5 ; 2[

    facteur_post_interpolation = ratio2;

    // Structure polyphase avec la moitié des coefs nuls
    ArrayXf coefs = design_rif_fen(15, "lp", 0.25, "hn");

    décimateurs.resize(nb_décimateurs);
    for(auto &d: décimateurs)
      d = filtre_rif_demi_bande<float, T>(coefs);
    // TODO : idem décimateur, moitié des coefs nuls
    suréchantilloneurs.resize(nb_suréchantilloneurs);
    for(auto &d: suréchantilloneurs)
      d = filtre_rif_ups<float, T>(coefs, 2);

    //float fcut = 0.4;
    //if(ratio2 < 1)
    float fcut = std::min(0.4f, facteur_post_interpolation / 2);

    auto itrp = itrp_sinc<T>({15, 256, fcut, "hn"});
    interpolateur = filtre_itrp<T>(facteur_post_interpolation, itrp);

    VERB(msg(" nb décimateurs = {}, nb suréchantilloneurs = {}, post-interpolation = {:.4f}.",
              nb_décimateurs, nb_suréchantilloneurs, facteur_post_interpolation););

    return 0;
  }

  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    y = x;

    if(ratio == 1)
      return;

    // Normamelement, à peu près deux fois moins de samples à chaque passe
    for(auto &d: décimateurs)
      y = d->step(y);

    for(auto &s: suréchantilloneurs)
      y = s->step(y);

    // Interpolateur
    if((facteur_post_interpolation < 1.000001)
        && (facteur_post_interpolation > 0.9999))
      return;

    y = interpolateur->step(y);
  }


};

template<typename T> sptr<FiltreGen<T>> filtre_reechan(float ratio)
{
  return std::make_shared<AdaptationRythmeArbitraire<T>>(ratio);
}

template<typename T> sptr<FiltreGen<T>> filtre_itrp(float ratio, sptr<Interpolateur<T>> itrp)
{
  return std::make_shared<AdaptationRythmeSimple<T>>(ratio, itrp);
}

namespace hidden {
auto filtre_ra1  = filtre_reechan<float>;
auto filtre_ra2  = filtre_reechan<cfloat>;
auto filtre_ra3  = filtre_itrp<float>;
auto filtre_ra4  = filtre_itrp<cfloat>;
}

}
