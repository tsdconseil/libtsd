#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/moniteur-cpu.hpp"


/** @brief Transformées classiques : TFR, TQC, ...  */
namespace tsd::fourier
{

/** @addtogroup fourier
 *  @{
 */

  /** @brief Classe abstraite pour un algorithme de calcul de FFT. */
  struct FFTPlan
  {
    /** @brief Configuration
     *  @param n          Nombre de points
     *  @param avant      Si faux, FFT inverse
     *  @param normalize  Si vrai, facteur de normalisation @f$1/\sqrt{N}@f$
     */
    virtual void configure(entier n, bouléen avant, bouléen normalize = oui) = 0;

    /** @brief Calcul de la FFT ou IFFT. */
    virtual void step(const Veccf &x, Veccf &y, bouléen avant = oui) = 0;

    Veccf step(const Veccf &x, bouléen avant = oui){Veccf y(x.dim()); step(x,y,avant); retourne y;}
  };

  /** @cond undoc */
  extern fonction<sptr<FFTPlan>()> fftplan_defaut;
  /** @endcond */

/** @brief Création d'un plan de calcul FFT (pour calculer efficacement plusieurs FFT).
 *
 * Cette fonction permet de créer un bloc qui sera efficace pour calculer
 * plusieurs FFT de même dimension (les facteurs de rotation sont calculés une seule fois).
 *
 * @param n           Dimension des vecteurs (paramètre optionnel)
 * @param avant       Si vrai, transformée directe, sinon transformée inverse.
 * @param normaliser  Si vrai, alors la transformée inclut un facteur de normalisation en @f$1/\sqrt{N}@f$ de manière à préserver l'énergie.
 * @return            Un pointeur vers une structure abstraite de type @ref FFTPlan.
 *
 * Le paramètre n n'a pas besoin d'être précisé, et la dimension des vecteurs peut changer en cours de route, si c'est le cas, les facteurs de rotations sont automatiquement réévalués.
 *
 * @note Si les échantillons d'entrée sont réels,
 * alors un plan plus efficace existe : @ref rtfrplan_création().
 *
 * @par Exemple
 * @code
 * entier N = 1024; // Dimension des blocs
 * auto plan = fftplan_création();
 * for(i = 0; i < n; i++)
 * {
 *    soit x = Vecf::random(N);
 *    // X est la FFT de x
 *    soit X = plan->step(x);
 *    // (équivalent à X = fft(x), mais plus efficace)
 * }
 *
 * @endcode
 *
 * @sa rfftplan_création(), fft(), ifft()
 **/
  extern sptr<FFTPlan> tfrplan_création(entier n = -1, bouléen avant = oui, bouléen normaliser = oui);

/** @brief Création d'un plan de calcul FFT pour des signaux réels (pour calculer efficacement plusieurs FFT).
 *
 * Cette fonction vous permet de créer un filtre qui sera efficace pour calculer
 * plusieurs FFT <b>sur des vecteurs réels</b> de même dimension (les facteurs de rotation sont calculés une seule fois).
 *
 * @param n       Dimension des vecteurs (paramètre optionnel)
 * @return Un %Filtre float vers @ref cfloat.
 *
 * Le paramètre n n'a pas besoin d'être précisé, et la dimension des vecteurs peut changer en cours de route, si c'est le cas, les facteurs de rotations sont automatiquement réévalués.
 *
 *
 * @par Exemple
 * @code
 * entier N = 1024; // Dimension des blocs
 * auto plan = rfftplan_création();
 * for(i = 0; i < n; i++)
 * {
 *    // Le signal d'entrée doit être réel
 *    soit x = Vecf::random(N);
 *    // X est la FFT de x
 *    soit X = plan->step(x);
 *    // (équivalent à X = rfft(x), ou X = fft(x), mais plus efficace)
 * }
 *
 * @endcode
 *
 * @sa fftplan_création(), fft(), ifft()
 */
extern sptr<FiltreGen<float, cfloat>> rtfrplan_création(entier n = -1);


/** @brief TFD d'un vecteur réel.
 *
 *  Cette fonction calcule de manière efficace la TFD d'un vecteur réel <b>de dimension paire</b>
 *  (si le vecteur fourni comporte un nombre impair d'éléments, la routine standard de la TFR/FFT est appelée).
 *
 *  @param x Vecteur réel
 *  @returns TFD (vecteur complexe)
 *
 *  @note La fonction @ref fft() utilise automatiquement cette routine
 *  si le signal d'entrée est de type réel.
 *
 *  @sa fft(), ifft()
 *
 */
template<typename T>
Veccf rfft(const Vecteur<T> &x)
{
  static_assert(!est_complexe<T>(), "RFFT : le vecteur d'entrée ne peut pas être complexe.");
  auto plan = tsd::fourier::rtfrplan_création(x.rows());
  return plan->step(x);
}

/** @brief Ré-échantillonage zéro phase à partir de la TFD.
 *
 * Cette fonction change la fréquence d'échantillonnage d'un signal, sans délais, grâce à un calcul dans le domaine fréquentiel.
 *
 * @param x     Vecteur d'entrée.
 * @param ratio Rapport entre la fréquence d'échantillonnage de sortie et celle d'entrée
 * (un ratio supérieur à 1 indique une interpolation, tandis que dans le cas contraire, on a affaire à une décimation).
 *
 * @warning Notez que, du fait des hypothèses sous-jacentes à la TFD, des artefacts peuvent apparaitrent au bord du signal.
 * Pour atténuer ce phénomène, vous pouvez pré-appliquer une fenêtre ou bien insérer des zéros avant et après le signal
 * (ce ne sont cependendant pas des solutions parfaites).
 *
 * @par Exemple illustrant les artefacts aux bords du signal
 * @snippet exemples/src/fourier/ex-fourier.cc exemple_resample_freq
 * @image html fourier-resample.png
 *
 * @sa rééchan()
 */
template<typename T>
  Vecteur<T> rééchan_freq(const Vecteur<T> &x, float ratio);

/** @brief Transformée de Fourier Discrète (TFD) rapide
 *
 *  Cette fonction calcule la TFD <b>normalisée</b> d'un vecteur réel ou complexe :
 *  @f[
 *    X_n = \frac{1}{\sqrt{N}}\cdot\sum_{k=0}^{N-1} x_k \cdot e^{\frac{-2\pi\mathbf{i}kn}{N}}
 *  @f]
 *
 *  @param x  Vecteur d'entrée : domaine temporel, réel ou complexe.
 *  @returns  X = TFD(x) : domaine fréquentiel, complexe.
 *
 *  @par Exemple
 *  @code
 *  Vecf  x = randn(100); // Vecteur réel
 *  Veccf X = fft(x);     // Vecteur complexe
 *  @endcode
 *
 *  @sa rfft(), ifft()
 */
template<typename T>
auto fft(const Vecteur<T> &x)
{
  Si constexpr(est_complexe<T>())
    retourne (tfrplan_création())->step(x);
  sinon
    retourne rfft(x);
}

/** @brief Fréquences normalisées associées à chaque bin d'une TFD / FFT.
 *
 *  @param n          Dimension de la TFD / FFT
 *  @param avec_shift Si vrai, suppose que fftshift a été appelé (fréquence négatives au début).
 *
 */
extern Vecf tfd_freqs(entier n, bouléen avec_shift = non);


/** @brief Transformée de Fourier Discrète Inverse
 *
 *  Cette fonction calcule la TFD inverse <b>normalisée</b> d'un signal réel ou complexe :
 *  @f[
 *    X_n = \frac{1}{\sqrt{N}}\cdot\sum_{k=0}^{N-1} x_k \cdot e^{\frac{+2\pi\mathbf{i}kn}{N}}
 *  @f]
 *
 *  @par Exemple
 *  @code
 *  soit x1 = randcn(100),
 *       X  = fft(x1),
 *       x2 = ifft(X);
 *
 *  assert((x2-x1).norm2() < 1e-7);
 *  @endcode
 *
 *  @sa fft(), rfft()
 *
 */
template<typename T>
auto ifft(const Vecteur<T> &X)
{
  soit plan = tsd::fourier::tfrplan_création();
  return plan->step(X.as_complex(), non);
}

/** @brief Décalage du spectre de manière à centrer les basses fréquences au milieu.
 *
 *  Cette fonction décale le spectre passé en entrée de manière à centrer le 0 Hz au milieu.
 *  Plus précisémment, si le nombre d'échantillons @f$N@f$ est pair :
 *  @f[
 *   X'_k = X_{k+N/2\ [N]}
 *  @f]
 *  soit :
 *  @f[
 *   X' = \left[X_{N/2}, X_{N/2+1}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2-1}\right]
 *  @f]
 *
 *  Et si @f$N@f$ est impair :
 *  @f[
 *   X' = \left[X_{N/2+1}, X_{N/2+2}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2}\right]
 *  @f]
 *
 *
 *  @param X spectre d'entrée, non centré (typiquement la sortie de la fonction @ref fft)
 *  @returns Spectre centré autour de 0 Hz.
 *
 *  @par Exemple
 *  @snippet exemples/src/fourier/ex-fourier.cc exemple_fftshift
 *  @image html fftshift.png
 **/
template<typename T>
Vecteur<T> fftshift(const Vecteur<T> &X)
{
  auto n = X.rows();
  Vecteur<T> res = Vecteur<T>::zeros(n);
  if((n & 1) == 0)
  {
    res.tail(n/2) = X.head(n/2);
    res.head(n/2) = X.tail(n/2);
  }
  else
  {
    res.tail(1+n/2) = X.head(1+n/2);
    res.head(n/2)   = X.tail(n/2);
  }
  return res;
}

/** @brief Modifie un vecteur de manière à ce qu'il soit conjugé symétrique.
 *
 *  @param X Vecteur complexe (modifié par la fonction)
 *
 *  Ecrase la deuxième partie du vecteur (qui correspond aux fréquences négatives),
 *  ainsi que la partie imaginaire du premier élément,
 *  et de l'élément médian si @f$n@f$ est pair.
 *  de manière à forcer la symétrie conjugée :
 *  @f[
 *  X_{N-k} = X_k^\star
 *  @f]
 *
 *  @note La TFD inverse d'un tel signal sera forcément réelle.
 */
template<typename T>
  void csym_forçage(Vecteur<T> &X)
{
  // Rien à faire dans ce cas
  if constexpr (!est_complexe<T>())
    return;
  else
  {
    soit n = X.rows();
    X(0).imag(0);

    si((n & 1) == 0)
      X(n/2).imag(0);
    sinon
      X(n/2+1) = conj(X(n/2));

    X.tail(n/2-1) = X.segment(1,n/2-1).reverse().conjugate();
  }
}

/** @cond undoc */

  // A SUPPRIMER

/** @brief Calcul de la cohérence spectrale entre deux vecteurs
 *
 *  Calcul de la cohérence spectrale entre deux vecteurs :
 *  @f[
 *  C_k = \frac{\left|X_k \cdot Y_k^\star\right|}{\left|X_k\right|\cdot \left|Y_k\right|} = \cos \left(\widehat{X_k, Y_k}\right)
 *  @f]
 *
 *  @f$X@f$ et @f$Y@f$ étant les TFD de @f$x@f$ et @f$y@f$.
 *
 *  Les @f$C_k@f$ sont donc compris entre -1 et 1, et indiquent à quel point les deux transformées de Fourier sont alignées.
 *
 *  @todo à tester et trouver une appli, un exemple */
extern Vecf coherence(const Veccf &x, const Veccf &y);

/** @endcond */

/** @brief Structure de configuration pour un filtre FFT (voir @ref filtre_fft()). */
struct FiltreFFTConfig
{
  /** @brief Nombre d'échantillons par bloc temporel (ou zéro si déterminé automatiquement) */
  entier dim_blocs_temporel = 0;

  /** @brief Nombre minimum de zéros insérés en début de bloc (<i>"zéro-padding"</i>) */
  entier nb_zeros_min = 0;

  /** @brief Si vrai, une fenêtre de Hann est utilisée, avec un recouvrement de 1/2
   *  pour un reconstruction parfaite. */
  bouléen avec_fenetrage = non;

  /** @brief Callback appelée pour le traitement dans le domaine fréquentiel
   *  (à réaliser par l'utilisateur). */
  fonction<void (Veccf &)> traitement_freq;
};

/** @brief Création d'un filtre dans le domaine fréquentiel (technique OLA / OverLap-and-Add).
 *
 *  Pour différentes raisons, il peut être plus intéressant de filtrer dans le domaine fréquentiel
 *  que temporel, par exemple pour alléger la charge de calcul
 *  (les convolutions dans le domaine temporel deviennent de simples
 *  produits termes à termes en fréquentiel).
 *  Par filtrage fréquentiel, on entends le fait de pouvoir modifier terme par terme les composantes
 *  de la transformée de Fourier d'un signal, puis de pouvoir restituer un signal temporel.
 *
 *  Afin de travailler dans le domaine fréquentiel, les données sont regroupées par blocs
 *  de @f$N_e@f$ échantillons, et pour chacun de ces blocs les opérations suivantes sont effectuées :
 *   - <b>Fenêtrage et recouvrement</b> (optionnel) : Afin de lisser les traitements d'un blocs à
 *      l'autre, une fenêtre d'apodisation peut être appliquée (par défaut,
 *      une fenêtre de Hann avec un recouvrement de 50 %, de manière à avoir une reconstruction parfaite).
 *   - <b>Zéro-padding</b> : De manière à avoir pouvoir calculer efficacement les FFT,
 *   des zéros sont ajoutés afin d'avoir des blocs dont la dimension est une puissance de 2. Après zéro-padding, la dimension des blocs vaut @f$N=2^k\geq N_e@f$.
 *   - <b>Transformée de Fourier :</b> via une FFT.
 *   - <b>Traitement dans le domaine fréquentiel</b> : à faire par l'utilisateur (callback passée en paramètre, voir ci-dessous).
 *   - <b>Transformée de Fourier inverse : </b> via une FFT inverse.
 *   - <b>Recomposition des blocs temporels</b> : par OLA (<i>OverLap-and-Add</i>).
 *
 *
 *
 *
 *  @param config La structure de configuration (@ref FiltreFFTConfig) doit préciser la dimension @f$N_e@f$
 *  des blocs temporels,
 *  le nombre minimum de zéros à insérer avant chaque bloc et avant de passer dans le domaine fréquentiel,
 *  ainsi qu'une callback permettant à l'appelant de modifier le signal dans le domaine fréquentiel.
 *
 *  @return Un tuple de deux éléments :
 *    - le filtre OLA,
 *    - un entier indiquant la dimension @f$N@f$ des FFT (c'est-à-dire
 *      la dimension des vecteurs qui seront passés à la callback)
 *
 *
 *  @parblock
 *  @note Pour le problème plus spécifique du calcul efficace d'un filtre RIF (à partir de la FFT), utilisez plutôt
 *  @ref filtre_rif_fft(), qui est est basée sur cette fonction, et qui plus simple à utiliser.
 *  @endparblock
 *
 *  @par Exemple : Filtrage passe-bas
 *  Dans cet exemple, on filtre passe-bas un signal en mettant tout simplement à zéro la moitié des composantes fréquentielles.
 *  @snippet exemples/src/fourier/ex-ola.cc ex_ola_pb
 *  @image html ex-ola-pb.png
 *
 *
 *  @sa filtre_rif_fft(), filtre_rfft(), ola_complexité(), ola_complexité_optimise()
 */
extern tuple<sptr<Filtre<cfloat, cfloat, FiltreFFTConfig>>, entier> filtre_fft(const FiltreFFTConfig &config);



/** @brief Calcul de la complexité d'un filtre OLA, par échantillon d'entrée, en FLOPS.
 *
 *  @param M  Taille de filtre (ou motif à détecter)
 *  @param Ne Taille de bloc d'entrée
 *  @param[out] C  Complexité, en FLOPS
 *  @param[out] Nf Dimension de la FFT
 *  @param[out] Nz Nombre de zéros insérés pour chaque bloc de FFT (@f$N_f = N_e + N_z@f$)
 *
 *  @sa ola_complexité_optimise()
 */
extern void ola_complexité(entier M, entier Ne, float &C, entier &Nf, entier &Nz);

/** @brief Calcul des paramètres optimaux pour un filtre par OLA.
 *
 * Cette fonction calcul la taille de bloc optimale pour le calcul d'un filtrage par OLA (voir @ref filtre_fft()).
 * L'optimum est trouvé pour @f$N_e = 2^k - (M-1)@f$, @f$k@f$ étant déterminé par recherche exhaustive.
 *
 *  @param M  Taille de filtre (ou motif à détecter)
 *  @param[out] C  Complexité, en FLOPS
 *  @param[out] Nf Dimension de la FFT
 *  @param[out] Nz Nombre de zéros insérés pour chaque bloc de FFT (@f$N_f = N_e + N_z@f$)
 *  @param[out] Ne Taille de bloc d'entrée
 *
 *  @sa ola_complexité()
 *
 */
extern void ola_complexité_optimise(entier M, float &C, entier &Nf, entier &Nz, entier &Ne);

// Création d'un filtre dans le domaine fréquentiel (technique OLA / OverLap-and-Add) - signaux réels.
// *
// *  <h3>Création d'un filtre fréquentiel</h3>
// *
// *  Ce filtre est identique à @ref filtre_fft(), à l'exception du fait que les types d'entrées / sorties sont réels (implémentation plus efficace).
// *
// *  @sa filtre_fft()
//extern tuple<sptr<Filtre<float, float, FiltreFFTConfig>>, entier> filtre_rfft(const FiltreFFTConfig &config);

/** @brief Transformée en z-chirp
 *
 *  Calcule la transformée en z sur les points suivants :
 *    @f$z_0 \cdot W^n@f$, pour @f$n=0,1,...,m-1@f$
 *
 *  @param x  Séquence d'entrée
 *  @param m  Nombre de points d'évaluation
 *  @param W  Raison de la série géométrique
 *  @param z0 Terme initial
 *  @returns Les valeurs de la transformée en @f$z@f$ de @f$x@f$ :
 *  @f[
 *    X\left(z_0 * W^n\right) = \sum_{k=0}^{N-1} x_k \cdot z_0^k \cdot W^{nk},\quad \textrm{pour }n=0,1,...,m-1
 *  @f] */
extern Veccf czt(const Veccf &x, entier m, cfloat W, cfloat z0 = 1.0f);

/** @} */

/** @addtogroup fourier-corr
 *  @{
 */

/** @brief Corrélation circulaire (normalisée) entre deux vecteurs complexes.
 *
 *  Calcul de la corrélation circulaire (normalisée) entre deux signaux complexes (calcul efficace via des FFT) :
 * @f[
 * c_n = \frac{1}{N} \cdot \sum_{k=0}^{N-1} x_k y^{\star}_{k+n[N]},\ \ n=0\dots N-1
 * @f]
 *
 *  @param x Premier vecteur
 *  @param y Deuxième vecteur (éventuellement égal au premier pour avoir l'auto-corrélation)
 *  @returns Premier vecteur : index des retards @f$n@f$ (soit @f$0, 1, ..., N-1@f$), deuxième vecteur : @f$c_n@f$.
 *
 *
 * @sa xcorr(), xcorrb()
 */
extern tuple<Vecf, Veccf> ccorr(const Veccf &x, const Veccf &y = Veccf());


/** @brief Corrélation (non biaisée) entre deux vecteurs complexes
 *
 *  Calcule la corrélation entre deux signaux complexes (calcul efficace via des FFT),
 *  soit pour un délais négatif (@f$\tau=-(m-1)\dots -1@f$) :
 *  @f[
 *  c(\tau) = \frac{1}{N-\tau} \cdot \sum_{k=0}^{N-1-\tau} x_{k-\tau} y^{\star}_k
 *  @f]
 *  Et, pour les délais positifs (@f$\tau=0\dots m-1@f$):
 *  @f[
 *  c(\tau) = \frac{1}{N-\tau} \cdot \sum_{k=0}^{N-1-\tau} x_k y^{\star}_{k+\tau}
 *  @f]
 *
 *  Les délais étant échantillonnés de manière régulière entre @f$-(m-1)@f$ et @f$m-1@f$ :
 *  @f[
 *   \tau = \textrm{linspace}(-(m-1), m-1, 2m-1) = -(m-1), -(m-2), \dots, -1, 0, 1, \dots, (m-2), (m-1)
 *  @f]
 *
 *
 *
 *  @param x Premier vecteur
 *  @param y Deuxième vecteur (éventuellement égal au premier pour avoir l'auto-corrélation)
 *  @param m Nombre de délais à examiner (si négatif, alors fait comme si @f$m=N@f$).
 *  @returns Un tuple de deux vecteurs de dimension @f$2m-1@f$ : délais @f$\tau_k@f$ et corrélation @f$c(\tau_k)@f$.
 *
 *
 * @note Les signaux sont complétés en interne avec des zéros pour pouvoir effectuer efficacement la corrélation dans le domaine fréquentiel.
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_xcorr
 * @image html xcorr.png
 *
 * @sa xcorrb(), ccorr(), détecteur_création()
 */
extern tuple<Vecf, Veccf> xcorr(const Veccf &x, const Veccf &y = Veccf(), entier m = -1);


/** @brief Corrélation (biaisée) entre deux vecteurs complexes.
 *
 *  Calcule la corrélation entre deux signaux complexes (calcul efficace via des FFT) :
 *  @f[
 *  c_n = \frac{1}{N} \cdot \sum_{k=0}^{N-1-k} x_k y^{\star}_{k+n},\ \ n=-m\dots m-1
 *  @f]
 *
 *  @param x Premier vecteur
 *  @param y Deuxième vecteur (éventuellement égal au premier pour avoir l'auto-corrélation)
 *  @param m Nombre de délais à examiner (si négatif, alors fait comme si @f$m=N@f$).
 *  @returns Un tuple de deux vecteurs de dimension @f$2m@f$ : délais et corrélation.
 *
 *  @sa xcorr(), ccorr(), détecteur_création()
 */
extern tuple<Vecf, Veccf> xcorrb(const Veccf &x, const Veccf &y = Veccf(), entier m = -1);

/** @brief Délais entier (décalage du signal) ou fractionnaire (basé sur la FFT).
 *
 *  Réalise un délais entier (décalage du signal) ou fractionnaire (basé sur la FFT) :
 *  @f[
 *  y_k = x_{k-\tau}
 *  @f]
 *
 *
 *  Le délais peut être positif (pour retarder le signal), ou négatif (pour avancer le signal).
 *
 *  @note Si le délais passé en paramètre est un entier, un simple décalage est effectué.
 *  Sinon, le délais est appliqué dans le domaine fréquentiel (après avoir complété le signal d'entrée avec des zéros afin d'atténuer les effets de bords dus à la périodicité de la FFT) :
 *  @f[
 *  Y_n = X_n \cdot e^{\frac{-2\pi\mathbf{i}\tau \cdot n}{N}}
 *  @f]
 *
 *  @param x Signal à retarder ou avancer dans le temps.
 *  @param τ Délais à appliquer, en nombre d'échantillons (peut être positif ou négatif, et n'est pas forcément un nombre entier).
 *  @returns Signal décalé dans le temps.
 *
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_delais
 * @image html ex-delais.png
 *
 **/
template<typename T = float>
  Vecteur<T> délais(const Vecteur<T> &x, float τ);


// Fonction supprimée
//extern entier estimation_délais_entier(IArrayXcf x, IArrayXcf y, float &score);


/** @brief Estimation du délais le
 *  plus probable entre deux signaux (via une corrélation) */
extern tuple<float, float> estimation_délais(const Veccf &x, const Veccf &y);


/** @brief Alignement de deux signaux */
template<typename T>
  tuple<Vecteur<T>, Vecteur<T>, entier, float> aligne_entier(const Vecteur<T> &x, const Vecteur<T> &y);


/** @brief Informations calculées à partir du motif détecté. */
struct Detection
{
  /** @brief Position en nombre d'échantillons (compris entre @f$-N_e@f$ et @f$N_e-1@f$),
  *    depuis le début du bloc de données en cours. Par exemple :
  *    - @f$0\  \Leftrightarrow@f$    Début du bloc en cours
  *    - @f$1\  \Leftrightarrow@f$    Deuxième échantillon du bloc en cours
  *    - @f$-1\  \Leftrightarrow@f$   Dernier échantillon du bloc précédent
  *    - ... */
  entier position;

  /** @brief Idem position, avec interpolation quadratique pour plus de précision */
  float position_prec;

  /** @brief Valeur absolue de la corrélation normalisée (valeur positive, entre 0 et 1) */
  float score;

  /** @brief Gain du signal (par rapport au motif passé en paramètre) */
  float gain;

  /** @brief Déphasage du signal (entre @f$-\pi/2@f$ et @f$\pi/2@f$) */
  float θ;

  /** @brief SNR estimé */
  float SNR_dB;

  /** @brief Ecart-type du bruit */
  float σ_noise;
};


/** @brief Structure de configuration pour un corrélateur par FFT */
struct DetecteurConfig
{
  /* @brief Nombre minimum d'échantillons / par bloc temporel (si 0 : déterminé automatiquement). */
  uint32_t Ne = 0;

  /* @brief Motif à détecter */
  Veccf motif;

  /* @brief Seuil de détection, entre 0 et 1. */
  float seuil = 0.5;

  bouléen debug_actif = non;

  /** @brief Mode de calcul du produit de corrélation : FFT avec OLA ou simple filtre RIF */
  enum Mode
  {
    MODE_OLA = 0,
    MODE_RIF = 1
  } mode = MODE_OLA;

  /** @brief Callback utilisateur appellée à chaque fois que le motif est détecté. */
  fonction<void (const Detection &det)> gere_detection;

  bouléen calculer_signal_correlation = non;
};

/** @cond undoc */
extern std::ostream &operator <<(std::ostream &os, const Detection &det);
/** @endcond */


/** @brief Structure abstraite pour un corrélateur à base de FFT ou de filtre RIF */
struct Detecteur: Filtre<cfloat, float, DetecteurConfig>
{
  virtual MoniteursStats moniteurs() = 0;
};

/** @brief Détecteur par corrélation.
 *
 * Cette fonction calcule, au fil de l'eau, la corrélation normalisée entre un signal @f$(x_k)@f$
 * reçu en continu et
 * un motif fixe @f$(h_k)@f$ de dimension @f$M@f$ :
 * @f[
 * c_n = \frac{\left|\displaystyle{\sum_{k=0}^{M-1} x_{n+k} \cdot h_{k}^\star}\right|}{\displaystyle{ \sqrt{\left(\sum_{k=0}^{M-1}\left|x_{n+k}\right|^2\right) \left(\sum_{k=0}^{M-1}\left|h_k\right|^2\right)}}}
 * @f]
 *
 * Les @f$c_n@f$ sont donc compris entre 0 et 1, et seront égaux à 1 uniquement si le signal est
 * exactement égal au motif (à un facteur d'échelle constant près) sur l'intervalle @f$[n\dots n + M[@f$.
 *
 *
 * Le calcul pouvant être fait soit grâce à un filtre RIF
 * (avec une complexité de l'ordre de @f$M@f$ opérations / échantillon),
 * soit de manière efficace dans le domaine fréquentiel
 * via la technique <i>Overlap-And-Add</i> (bloc @ref filtre_fft()),
 * la complexité étant alors de l'ordre de @f$\log_2(M)@f$ opérations
 * par échantillon (si @f$M@f$ est une puissance de 2).

 *
 * Ensuite, les pics isolés correspondant à des positions possibles de motif sont détectés.
 * Les deux critères suivants sont vérifiés~:
 * - Score @f$c_n@f$ supérieur à un seuil paramétrable,
 * - Absence de pic de score supérieur dans un voisinage de @f$\pm@f$ la moitié de la durée
     du motif de synchronisation.
 *
 * Pour chaque position possible, les paramètres suivants sont alors estimés :
 *   <b>Gain complexe du signal</b> (rapport complexe entre le signal reçu et le motif théorique attendu) :
 *     @f[
 *      g \sim \frac{\displaystyle{\sum_{k=0}^{M-1} x_{n+k} \cdot h_k^\star}}{\displaystyle{\sum_{k=0}^{M-1} |h_k|^2}}
 *     @f]
 *     De cette valeur complexe, on en déduit le gain en amplitude et le déphasage :
 *     @f[
 *     A = \left|g\right|,\quad \phi = \arg g
 *     @f]
 *  <b>Ecart-type du bruit et SNR :</b>
 *       Pour estimer le niveau de bruit, on note simplement que
 *       @f$x_n - g^\star h_n \sim x_n - g h_n = b_n@f$, ainsi il suffit de
 *       mesurer l'énergie du signal @f$x_n - g h_n@f$ pour avoir le niveau du bruit :
 *       @f[
 *       \epsilon_b \sim \mathbb{E}\left[\left|x_n - g h_n\right|^2\right]
 *       @f]
 *       Cependant, le signal reçu sera en général affecté d'un retard non entier qui,
 *       si il ne posera pas de problème pour l'estimation de gain (puisque le
 *       signal est suréchantillonné), pourrait conduire sur-estimer le bruit. C'est pourquoi avant estimation de celui-ci,
 *       le signal est décalé d'un pas fractionnaire afin de minimiser l'erreur.
 *
 * Pour chaque motif détecté, une callback utilisateur est appelée,
 * avec l'ensemble des paramètres estimés ainsi que le score de corrélation normalisé.
 *
 *
 * @par Schéma-bloc
 * <img src="detecteur.png" align="left" width="300px"/>
 * <div style="clear: both"></div>
 *
 * @param config Structure de configuration
 * @return Filtre cfloat @f$\to@f$ float, qui prend en entrée des échantillons, et sort, au fil de l'eau,
 * la corrélation normalisée avec le motif.
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_fft_correlateur
 * @image html fft_correlateur.png
 */
extern sptr<Detecteur>
  détecteur_création(const DetecteurConfig &config = DetecteurConfig());



// Alignement de deux signaux suivant un délais variable
//extern sptr<ProcesseurConfigurable<cfloat, cfloat, entier>> creation_aligneur();

// ?
struct AlignementSignal
{
  AlignementSignal();

  void configure(entier N);
  void step(const Tabcf &x, const Tabcf &y, entier delais);

  _PIMPL_
};




/** @} */



/** @addtogroup fourier-freq
 *  @{
 */

/** @brief Calcul des fréquences normalisées associée à une psd
 *
 *
 *  @param n Nombre de points dans le domaine fréquentiel
 *  @param est_complexe Spectre avec fréquences négatives et positives (<pre>est_complexe=oui</pre>) ou fréquences positives uniquement (<pre>est_complexe=non</pre>).
 *  @return Vecteur de fréquences normalisées, entre -0,5 et 0,5 ou 0 à 0,5.
 *
 *
 *  @sa psd() */
extern Vecf psd_freqs(entier n, bouléen est_complexe = oui);

/** @brief PSD (corrélogramme)
 *
 * Calcul du corrélogramme avec une fenêtre de Hann :
 *
 * @f[
 * S(k) = \left|TFD(x \cdot w)(k)\right|
 * @f]
 *
 * Le spectre est centré autour de 0 Hz (vecteur de fréquences : linspace(-0.5,0.5-1.0/N,N)).
 *
 * @param x   Signal à analyser.
 * @return Un tuple de deux éléments : vecteur de fréquences (normalisées, entre -0.5 et 0.5), et le spectre.
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_psd2
 *
 * PSD signal triangulaire :
 * @image html ex-psd2.png
 *
 * @sa psd_welch(), psd_sousesp()
 */
template<typename T>
tuple<Vecf, Vecf> psd(const Vecteur<T> &x)
{
  soit fen = tsd::filtrage::fenêtre("hn", x.rows(), non);
  soit xf = x * fen;

  if constexpr(est_complexe<T>())
  {
    soit Y = pow2db(abs2(fft(xf)));
    return {psd_freqs(x.rows(), oui), fftshift(Y)};
  }
  else
  {
    soit Y = pow2db(abs2(rfft(xf)));
    return {psd_freqs(x.rows(), non), Y.head(Y.rows() / 2).eval()};
  }
}

/** @brief PSD (méthode de Welch - par moyennage).
 *
 * Cette fonction est une estimation de la PSD par la technique du moyennage :
 * plusieur PSD sont calculées sur des fenêtres (pondérées et avec du recouvrement), puis moyennées.
 *
 * @param x   Signal à analyser.
 * @param N   Dimension des sous-fenêtres (autrement dit, la résolution de la PSD résultante).
 * @param fen Choix du fenêtrage.
 * @return Un tuple de deux éléments : vecteur de fréquences (normalisées, entre -0.5 et 0.5), et le spectre.
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_psd3
 *
 * PSD et moyennage :
 * @image html ex-psd3.png
 *
 * @sa psd(), psd_sousesp()
 *
 */
extern tuple<Vecf, Vecf> psd_welch(const Veccf &x, entier N, cstring fen = "hn");


/** @brief Calcul d'un spectre par la méthode des sous-espaces
 *
 *  Cette technique d'estimation spectrale est une méthode paramétrique, qui suppose que le signal est constitué
 *  de la somme d'un nombre fini et <b>connu</b> d'exponentielles pures.
 *
 *  La formule utilisée ici est celle de l'algorithme MUSIC (MUltiple SIgnal Classification).
 *
 *  @param x  Signal à analyser
 *  @param Ns Nombre d'exponentielles pures attendues
 *  @param Nf Résolution fréquentielle (en nombre de points)
 *  @param m  Dimension de la matrice d'auto-corrélation (si non précisé = dimension du signal)
 *  @returns Un tuple de deux éléments : vecteur de fréquences (normalisées, entre -0.5 et 0.5), et le spectre.
 *
 *  @note Cette fonction est une spécialisation de la méthode des sous-espaces pour la détection d'exponentielles pures.
 *  Cependant, la méthode des sous-espaces elle-même est bien plus générale, et peut
 *  permettre de détecter d'autres types de signaux
 *  (voir la fonction @ref tsd::stats::subspace_spectrum()).
 *
 *  @warning Cette technique implique la décomposition en valeurs propres d'une matrice @f$m\times m@f$
 *  (la matrice d'auto-corrélation du signal), soit de l'ordre de @f$m^3@f$ opérations.
 *  Par conséquent, si la dimension du signal est importante,
 *  il est conseillé de choisir un @f$m@f$ plus petit.
 *
 *  @par Exemple : détection de 3 exponentielles pures
 *  @snippet exemples/src/fourier/ex-fourier.cc ex_psd_subspace
 *
 *  Spectre PSD vs MUSIC :
 *  @image html ex-subspace-freq-spectrum.png
 *
 *  Notez comme le spectre MUSIC est plus propre ; il ne faut cependant pas oublier
 *  que cette technique suppose connu le nombre de signaux (ici 3).
 *
 *  @sa tsd::stats::subspace_spectrum()
 *
 *  @note Pour la détection de sinusoïdes pures (signaux réels), il faudra choisir @f$N_s@f$ comme étant <b>deux fois</b>
 *  le nombre de signaux à détecter (en effet, un cosinus ou sinus peut se représenter comme la somme de deux exponentielles complexes).
 *
 */
extern tuple<Vecf, Vecf> psd_sousesp(const Veccf &x, entier Ns, entier Nf = 1024, entier m = 0);//, const SubSpaceSpectrumConfig &config);


/** @brief Choix d'un algorithme pour l'estimation de fréquence */
enum FreqEstimMethode
{
  /** @brief Maximum de la valeur absolue de la FFT */
  FFT,
  /** @brief Méthode de Candan cf TODO */
  CANDAN2,
  /** @brief Interpolation quadratique */
  QUADRATIC
  // TODO
};

/** @brief Estimation de fréquence d'un signal périodique
 *
 *  Cette fonction essaye d'estimer la fréquence d'un signal périodique (exponotielle pure).
 *
 *  @note Si le signal d'entrée est réel (sinusoide pure de fréquence @f$f@f$),
 *  la fonction renvoie une des deux valeurs possibles de fréquence (@f$f@f$ ou @f$-f@f$).
 *
 *  @warning Si le signal d'entrée est composé de plusieurs fréquences distinctes, cette fonction
 *  n'est pas adaptée. Utilisez plutôt un algorithme de type sous-espace (voir @ref psd_sousesp()).
 *
 *
 *  @param x Signal dont on veut estimer la fréquence
 *  @param m Méthode d'estimation (voir l'énumération @ref FreqEstimMethode)
 *  @return Une fréquence normalisée (entre -0,5 et 0,5) */
extern float freqestim(const Veccf &x, FreqEstimMethode m = CANDAN2);




/** @brief %Filtre de Goertzel (calcul sur un buffer).
 *
 *  Calcul de la densité spectrale, normalisée par rapport à l'énergie totale du signal, à la fréquence @f$f@f$ :
 *  @f[
 *  y = \frac{2\cdot\left|X(f)\right|^2}{\sum{x_k^2}}
 *  @f]
 *
 *  La normalisation permet d'obtenir un score entre 0 (énergie nulle pour la fréquence considérée) et 1
 *  (sinusoïde pure de fréquence @f$f@f$).
 *
 *  Le calcul de @f$X(f)@f$ est effectué grâce à un filtre récursif réel d'ordre 2, moins coûteux que le calcul d'une TFD complète.
 *
 *  @note Le facteur 2 dans la formule ci-dessous est nécessaire pour compenser le fait que pour un signal réel, l'énergie est partagée entre les fréquences positives et négatives.
 *  On pourrait d'ailleurd écrire de manière équivalente :
 *  @f[
 *  y = \frac{\left|X(f)\right|^2 + \left|X(-f)\right|^2}{\sum{x_k^2}}
 *  @f]
 *
 *  @param x         Signal à analyser
 *  @param frequence Fréquence à  détecter (normalisée à  la fréquence d'échantillonnage, entre -0,5 et 0,5).
 *  @returns         Score normalisé entre 0 et 1.
 *  @sa filtre_goertzel() */
extern float goertzel(const Vecf &x, float frequence);

/** @brief %Filtre de Goertzel (calcul au fil de l'eau)
 *
 *  Ce filtre calcule au fil de l'eau la densité de puissance (normalisée à la puissance totale du signal)
 *  à une fréquence donnée (voir @ref goertzel()).
 *
 *  Le calcul est fait par blocs de @f$N@f$ échantillons, et la fréquence d'échantillonnage
 *  en sortie est donc divisée par @f$N@f$ par rapport
 *  à la fréquence d'entrée (1 échantillon produit pour @f$N@f$ échantillons consommés) :
 *  @f[
 *  y_i = \frac{\left|2\cdot\mathcal{F}(x_{Ni\dots Ni + N-1})(f)\right|^2}{\displaystyle{\sum_{k=Ni}^{Ni+N-1}{x_k^2}}}
 *  @f]
 *
 *  @f$\mathcal{F}(.)@f$ étant la transformée de Fourier discrète.
 *
 *  Le paramètre @f$N@f$ correspond donc au nombre de points de la transformée de Fourier, et
 *  la résolution fréquentielle sera donc d'autant plus fine que @f$N@f$ est grand.
 *
 *  @param   frequence  Fréquence à  détecter (normalisée à  la fréquence d'échantillonnage).
 *  @param N            Nombre de points de la TFD équivalente (facteur de décimation)
 *  @returns            Filtre générique float -> float
 *
 *  @par Exemple : détection d'une sinusoïde pure
 *  @snippet exemples/src/fourier/ex-fourier.cc exemple_goertzel
 *  @image html goertzel.png
 *
 *  @sa goertzel()
 */
extern sptr<FiltreGen<float>> filtre_goertzel(float frequence, entier N);




/** @brief Configuration pour le calcul de spectre en temps réel */
struct SpectrumConfig
{
  /** @brief Dimension des blocs d'entrée */
  entier BS = 1024;

  /** @brief Nombre de spectres moyennés */
  entier nmeans = 10;

  /** @brief Pour le multi-threading ou le balayage fréquentiel */
  entier nsubs  = 1;

  /** @brief Configuration du balayage fréquentiel */
  struct
  {
    /** @brief Activation */
    bouléen active = non;

    /** @brief Nombre de pas de fréquence */
    entier  step   = 1024;

    /** @brief Nombres de points à masquer en basse fréquence */
    entier masque_bf = 0;

    /** @brief Nombres de points à masquer en haute fréquence */
    entier masque_hf = 0;

  } sweep;

  /** @brief Dim FFT */
  entier Nf() const;

  /** @brief Dim spectre */
  entier Ns() const;

  /** @brief Type de fenêtrage */
  tsd::filtrage::Fenetre fenetre = tsd::filtrage::Fenetre::HANN;

  /** @brief Plan de calcul FFT (optionnel) */
  sptr<FFTPlan> plan;
};


/** @brief Calcul de spectre en temps réel */
sptr<Filtre<cfloat,float,SpectrumConfig>> rt_spectrum(const SpectrumConfig &config);

/** @} */

} // namespace tsd::fourier

namespace tsd::tf {

/** @addtogroup temps-frequence
*  @{ */

/** @brief Calcul d'une matrice temps / fréquence par FFT glissante
 *
 *  @param x Signal à analyser
 *  @param N Dimension des blocs d'analyse */
extern Tabf periodogramme_tfd(const Veccf &x, entier N);

/** @brief Calcul d'un spectrogramme à partir de la CQT (Constant Quality Transform).
 *
 *  @param x Signal à analyser
 *  @param fe Fréquence d'échantillonnage (Hz)
 *  @param fmin Fréquence minimale d'analyse (Hz)
 *  @param fmax Fréquence minimale d'analyse (Hz)
 *  @param γ    Rapport entre deux fréquences d'analyse successives
 *  @param ofs  Fréquence d'échantillonnage de sortie (Hz)
 */
extern Tabf periodogramme_cqt(const Vecf &x, float fe, float fmin, float fmax, float γ, float ofs);

/* @} */

}


ostream_formater(tsd::fourier::Detection)




