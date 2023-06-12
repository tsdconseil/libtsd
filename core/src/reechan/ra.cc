#include "tsd/filtrage.hpp"
#include "tsd/filtrage/spline.hpp"

#define VERB(AA)
#define RA_SAFE_MODE 0


using namespace std;

namespace tsd::filtrage {


template<typename T>
struct AdaptationRythmeSimple: FiltreGen<T>
{
  float phase = 0,
        ratio = 1,
  /** Inverse du ratio precedent */
        increment = 1;

  Vecteur<T> fenetre;
  entier nfen = 0;
  sptr<Interpolateur<T>> interpolateur;

  AdaptationRythmeSimple(float ratio, sptr<Interpolateur<T>> itrp)
  {
    this->ratio   = ratio;
    increment     = 1 / ratio;
    phase         = 0;
    interpolateur = itrp;
    nfen          = itrp->K;
    fenetre.setZero(nfen);
    si(nfen <= 1)
      msg_avert("AdaptationRythmeSimple : npts interpolation = {} ({}).", nfen, itrp->nom);
  }



  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n       = x.rows(),
         j       = 0,
         tmp_len = (entier) (ceil(ratio * n) + 10);

    si(n == 0)
    {
      y.resize(0);
      retourne;
    }

    Vecteur<T> temp(tmp_len);

    soit optr = temp.data();
    soit iptr = x.data();

    pour(auto i = 0; i < n; i++)
    {
      // Rotation de la fenêtre
      // TODO : ça serait pas mal que les interpolateurs puissent fonctionner de manière circulaire
      // pour éviter cette rotation
      fenetre.head(nfen-1) = fenetre.tail(nfen-1).clone();
      fenetre(nfen-1) = *iptr++;

      tantque(phase < 1)
      {
        // phase = index entre deux échantillons
        soit valeur_interpolée = interpolateur->step(fenetre, 0, phase);
        assertion(j < tmp_len);
        *optr++ = valeur_interpolée;
        j++;
        phase += increment;
      }
      phase--;
    } // end pour i
    y.resize(j);
    y = temp.head(j);
  }

};



// Adaptation de rythme, ratio arbitraire
template<typename T>
struct AdaptationRythmeArbitraire: Filtre<T,T,float>
{
  sptr<FiltreGen<T>> interpolateur;
  vector<sptr<FiltreGen<T>>> décimateurs, suréchantilloneurs;

  /** Ratio d'interpolation / post-décimation */
  float facteur_post_interpolation;

  /** Decimation / interpolation ratio.
      e.g. d < 1.0 pour decimation, d > 1.0 pour interpolation */
  float ratio;


  AdaptationRythmeArbitraire(float ratio)
  {
    Configurable<float>::configure(ratio);
  }

  /** Setup the rate adaptation */
  void configure_impl(const float &ratio_) override
  {
    ratio = ratio_;

    si((ratio <= 0) || isinf(ratio) || (ratio >= 1e9))
    {
      msg_erreur("AdaptationRythmeArbitraire::configurer() : facteur de décimation invalide : {}.", ratio);
      ratio = 1;
    }

    VERB(msg("ratio de décimation / interpolation demandé : {}.", ratio_);)

    facteur_post_interpolation = ratio;
    soit nb_suréchantilloneurs = 0,
         nb_décimateurs = 0;


    // Place le ratio d'interpolation dans l'intervalle [0,5 ; 2[
    tantque(facteur_post_interpolation < 0.5)
    {
      nb_décimateurs++;
      facteur_post_interpolation *= 2;
    }

    tantque(facteur_post_interpolation >= 2)
    {
      nb_suréchantilloneurs++;
      facteur_post_interpolation /= 2;
    }

    // Structure polyphase avec la moitié des coefs nuls
    // TODO : design plus efficace (fir_eq, etc.)
    soit coefs = design_rif_fen(15, "lp", 0.25, "hn");

    décimateurs.resize(nb_décimateurs);
    pour(auto &d: décimateurs)
      d = filtre_rif_demi_bande<float, T>(coefs);
    // TODO : idem décimateur, moitié des coefs nuls
    suréchantilloneurs.resize(nb_suréchantilloneurs);
    pour(auto &d: suréchantilloneurs)
      d = filtre_rif_ups<float, T>(coefs, 2);


    // TODO: documenter
    // [0.5 ; 2[
    soit fcut = min(0.4f, facteur_post_interpolation / 2);
    soit itrp = itrp_sinc<T>({15, 256, fcut, "hn"});
    // TODO: voir meilleur interpolateur?
    interpolateur = filtre_itrp<T>(facteur_post_interpolation, itrp);

    VERB(msg(" nb décimateurs = {}, nb suréchantilloneurs = {}, post-interpolation = {:.4f}.",
              nb_décimateurs, nb_suréchantilloneurs, facteur_post_interpolation););
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y) override
  {
    y = x;

    si(ratio == 1)
      retourne;

    // Normamelement, à peu près deux fois moins de samples à chaque passe
    pour(auto &d: décimateurs)
      y = d->step(y);

    pour(auto &s: suréchantilloneurs)
      y = s->step(y);

    // Interpolateur
    si(abs(facteur_post_interpolation - 1) < 1e-6f)
      retourne;

    y = interpolateur->step(y);
  }
};

template<typename T> sptr<Filtre<T,T,float>> filtre_reechan(float ratio)
{
  retourne make_shared<AdaptationRythmeArbitraire<T>>(ratio);
}

template<typename T> sptr<FiltreGen<T>> filtre_itrp(float ratio, sptr<Interpolateur<T>> itrp)
{
  retourne make_shared<AdaptationRythmeSimple<T>>(ratio, itrp);
}

namespace hidden {
soit filtre_ra1  = filtre_reechan<float>;
soit filtre_ra2  = filtre_reechan<cfloat>;
soit filtre_ra3  = filtre_itrp<float>;
soit filtre_ra4  = filtre_itrp<cfloat>;
}

}
