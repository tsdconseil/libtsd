#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "dsp/dsp.hpp"
#include "dsp/filter.hpp"
#include "tsd/moniteur-cpu.hpp"
#include "tsd/fourier.hpp"


/** @brief Signal transformations: FFT, CQT, ... */
namespace dsp::fourier
{
  namespace tsdF = tsd::fourier;


/** @addtogroup fourier
 *  @{
 */

  /** @brief Abstract class for a FFT computing algorithm. */
  using tsdF::FFTPlan;

  using tsdF::fftplan_defaut;

  /** @cond undoc */
  //extern std::function<sptr<FFTPlan>()> fftplan_defaut;
  /** @endcond */

/** @brief Creation of a FFT computing plan (to compute efficiently several FFT).
 *
 * <h3>Creation of a FFT computing plan</h3>
 *
 * This function will create a structure which will be efficient to compute sereral FFT
 * with the same data length (the twiddle factors are only computed once).
 *
 * @param n           Vector length (optional parameter)
 * @param forward     If true, direct transform, otherwise inverse transform.
 * @param normalize   If true, then the transform includes a normalization factor of @f$1/\sqrt{N}@f$ so as to have the same total energy for the input and output samples.
 * @return            A poin @ref Filtre cfloat vers cfloat.
 *
 * Note that the parameters n and forward can be left unspecified, or even changed later while calling the FFT plan.
 *
 * @note If the input samples are real,
 * then a more efficient plan is possible: @ref rfftplan_new().
 *
 * @par Example
 * @code
 * int N = 1024; // Bloc length
 * auto plan = fftplan_new();
 * for(i = 0; i < n; i++)
 * {
 *    ArrayXcf x = ArrayXcf::Random(N);
 *    // X is the FFT of x
 *    ArrayXcf X = plan->step(x);
 *    // (equivalent to X = fft(x), but more efficient)
 * }
 *
 * @endcode
 *
 * @sa rfftplan_new(), fft(), ifft()
 **/
inline sptr<FFTPlan> fftplan_new(int n = -1, bool forward = true, bool normalize = true)
{
  return tsdF::fftplan_création(n, forward, normalize);
}

/** @brief Creation of a FFT computing plan for real signals (to compute efficiently several FFT).
 *
 * <h3>Creation of a FFT computing plan (real signals)</h3>
 *
 * This function will create a structure which will be efficient to compute sereral FFT
 * <b>on real signals</b> with the same data length (the twiddle factors are only computed once).
 *
 * @param n       Vector length (optional parameter)
 * @return        A filter block (float to @ref cfloat).
 *
 * Note that the parameters n and forward can be left unspecified, or even changed later while calling the FFT plan.
 *
 *
 * @par Example
 * @code
 * int N = 1024; // Bloc length
 * auto plan = rfftplan_new();
 * for(i = 0; i < n; i++)
 * {
 *    // Input signal must be real
 *    ArrayXf x = ArrayXf::Random(N);
 *    // X is the FFT of x
 *    ArrayXcf X = plan->step(x);
 *    // (equivalent to X = fft(x) or X = rfft(x), but more efficient)
 * }
 *
 * @endcode
 *
 * @sa fftplan_new(), fft(), ifft()
 */
inline sptr<FiltreGen<float, cfloat>> rfftplan_new(int n)
{
  return tsdF::rfftplan_création(n);
}


/** @brief FFT of a real vector.
 *
 *  <h3>FFT of a real vector</h3>
 *
 *  This function computes efficiently the FFT of a real vector <b>with even number of samples</b>
 *  (if the provided vector has an odd number of elements, then the standard fft routine is called).
 *
 *  @note The @ref fft() function will use automatically this function if the input vector is real.
 */
template<typename derived>
ArrayXcf rfft(const Eigen::ArrayBase<derived> &x)
{
  return tsdF::rfft(x);
}

/** @brief DFT based zero-phase resampling.
 *
 * <h3>DFT based zero-phase resampling</h3>
 *
 * This function will resample a signal, without introducing any delay,
 * by working in the frequency domain.
 *
 * @param x     Input vector.
 * @param ratio Desired ratio between output and input sample rates.
 * (a ratio greater than 1 can be used for interpolation,
 * and a ratio less than 1 for a decimation).
 *
 * @warning Note that, because of some assumptions implicit with the DFT,
 * some artefacts may be created at the signal borders.
 * To mitigate this phenomena, one can pre-apply a window, or insert zeroes before and after the actual samples,
 * but these are not perfect solutions.
 *
 * @par Example showing the spurius artefacts
 * @snippet exemples/src/fourier/ex-fourier.cc exemple_resample_freq
 * @image html fourier-resample.png width=1000px
 *
 * @sa resample()
 */
template<typename T>
  Vector<T> resample_freq(const Vector<T> &x, float ratio)
{
  return tsdF::reechan_freq(x, ratio);
}

/** @brief Fast discrete Fourier Transform (FFT).
 *
 *  <h3>Fast discrete Fourier Transform (FFT)</h3>
 *
 *  This function will compute the <b>normalized</b> DFT of a real or complex vector:
 *  @f[
 *    X_n = \frac{1}{\sqrt{N}}\cdot\sum_{k=0}^{N-1} x_k \cdot e^{\frac{-2\pi\mathbf{i}kn}{N}}
 *  @f]
 *
 *  @param x    Input vector (time domain).
 *  @returns    X = DFT(x)
 *
 *  @par Example
 *  @code
 *  ArrayXf  x = randn(100);
 *  ArrayXcf X = fft(x);
 *  @endcode
 *
 *  @sa rfft(), ifft()
 */
template<typename derived>
auto fft(const Eigen::ArrayBase<derived> &x)
{
  return tsdF::fft(x);
}

/** @brief Inverse Fast discrete Fourier Transform
 *
 *  <h3>Inverse Fast discrete Fourier Transform</h3>
 *
 *  This function will compute the <b>normalized</b> DFT inverse of a real or complex vector:
 *  @f[
 *    X_n = \frac{1}{\sqrt{N}}\cdot\sum_{k=0}^{N-1} x_k \cdot e^{\frac{+2\pi\mathbf{i}kn}{N}}
 *  @f]
 *
 *  @par Example
 *  @code
 *  ArrayXcf x1 = ArrayXcf::Random(100);
 *  ArrayXcf X  = fft(x1);
 *  ArrayXf  x2 = ifft(X);
 *
 *  assert((x2-x1).norm2() < 1e-7);
 *  @endcode
 *
 *  @sa fft(), rfft()
 *
 *  */
template<typename derived>
auto ifft(const Eigen::ArrayBase<derived> &X)
{
  return tsdF::ifft(X);
}

/** @brief Spectrum shift so as to be centered on the low frequencies.
 *
 *  <h3>Spectrum shift</h3>
 *
 *  This function shift the input spectrum so as the 0 Hz bin is at the center.
 *  If the number of samples is even:
 *  @f[
 *   X'_k = X_{k+N/2\ [N]}
 *  @f]
 *  That is:
 *  @f[
 *   X' = \left[X_{N/2}, X_{N/2+1}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2-1}\right]
 *  @f]
 *
 *  And if @f$N@f$ is odd:
 *  @f[
 *   X' = \left[X_{N/2+1}, X_{N/2+2}, \dots, X_{N-1}, X_0, X_1, \dots, X_{N/2}\right]
 *  @f]
 *
 *
 *  @param X Input spectrum, uncentered (typicaly ouput of the @ref fft() function)
 *  @returns 0 Hz centered spectrum.
 *
 *  @par Example
 *  @snippet exemples/src/fourier/ex-fourier.cc exemple_fftshift
 *  @image html fftshift.png width=800px
 *
 */
template<typename derived>
Vector<typename derived::Scalar> fftshift(const Eigen::ArrayBase<derived> &X)
{
  return tsdF::fftshift(X);
}

/** @brief Change a vector so as it has the conjugate symetry proporty.
 *
 *  <h3>Forcing conjugate symetry</h3>
 *
 *  @param X Complex vector (input and output parameter)
 *
 *  Over-write the second half of the vector (which contain negative frequency bins),
 *  and also the imaginary part of the first element,
 *  and, if @f$n@f$ is even, the middle element, so that the
 *  resulting vector is conjugate symetrical:
 *  @f[
 *  X_{N-k} = X_k^\star
 *  @f]
 *
 *  @note The inverse DFT of such a signal will necessarily be real.
 */
template<typename derived>
  void force_csym(Eigen::ArrayBase<derived> &X)
{
  return tsdF::force_csym(X);
}




/** @brief Configuration structure for a FFT based filter (see @ref filter_fft()). */
struct FFTFilterConfig: tsdF::FiltreFFTConfig
{
  /** @brief Nombre d'échantillons par bloc temporel (ou zéro si déterminé automatiquement) */
  int &time_blocks_length = dim_blocs_temporel;

  /** @brief Nombre minimum de zéros insérés en début de bloc (<i>"zéro-padding"</i>) */
  int &minimum_zeros_count = nb_zeros_min;

  /** @brief Si vrai, une fenêtre de Hann est utilisée, avec un recouvrement de 1/2
   *  pour un reconstruction parfaite. */
  bool &enable_windowing = avec_fenetrage;

  /** @brief Callback appelée pour le traitement dans le domaine fréquentiel
   *  (à réaliser par l'utilisateur). */
  std::function<void (ArrayXcf &)> &freq_domain_processing = traitement_freq;
};

/** @brief Creation of frequency domain filter (OLA / OverLap-and-Add technique).
 *
 *  <h3>Creation of frequency domain filter</h3>
 *
 *  Pour différentes raisons, il peut être plus intéressant de filtrer dans le domaine fréquentiel
 *  que temporel, par exemple pour alléger la charge de calcul
 *  (les convolutions deviennent de simples produits termes à termes en fréquentiel).
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
 *  @param config La structure de configuration (@ref FFTFilterConfig) doit préciser la dimension @f$N_e@f$
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
 *  @ref filter_fir_fft(), qui est est basée sur cette fonction, et qui plus simple à utiliser.
 *  @endparblock
 *
 *  @par Exemple : Filtrage passe-bas
 *  Dans cet exemple, on filtre passe-bas un signal en mettant tout simplement à zéro la moitié des composantes fréquentielles.
 *  @snippet exemples/src/fourier/ex-ola.cc ex_ola_pb
 *  @image html ex-ola-pb.png width=800px
 *
 *
 *  @sa filtre_rif_fft(), filtre_rfft()
 */
inline std::tuple<sptr<Filtre<cfloat, cfloat, tsdF::FiltreFFTConfig>>, int> filter_fft
  (const FFTFilterConfig &config)
{
  return tsdF::filtre_fft(config);
}



/** @brief Calcul de la complexité d'un filtre OLA, par échantillon d'entrée, en FLOPS.
 *
 * <h3>Calcul de la complexité d'un filtre OLA</h3>
 *
 *  @param M  Taille de filtre (ou motif à détecter)
 *  @param Ne Taille de bloc d'entrée
 *  @param[out] C  Complexité, en FLOPS
 *  @param[out] Nf Dimension de la FFT
 *  @param[out] Nz Nombre de zéros insérés pour chaque bloc de FFT (Nf = Ne + Nz)
 *
 *  @sa ola_complexite_optimise()
 */
inline void ola_complexity(int M, int Ne, float &C, int &Nf, int &Nz)
{
  return tsdF::ola_complexite(M, Ne, C, Nf, Nz);
}

/** @brief Calcul des paramètres optimaux pour un filtre par OLA.
 *
 * <h3>Calcul des paramètres optimaux pour un filtre par OLA</h3>
 *
 * Cette fonction calcul la taille de bloc optimale pour le calcul d'un filtrage par OLA.
 * L'optimum est trouvé pour @f$N_e = 2^k - (M-1)@f$, @f$k@f$ étant déterminé par recherche exhaustive.
 *
 *  @param M  Taille de filtre (ou motif à détecter)
 *  @param[out] C  Complexité, en FLOPS
 *  @param[out] Nf Dimension de la FFT
 *  @param[out] Nz Nombre de zéros insérés pour chaque bloc de FFT (Nf = Ne + Nz)
 *  @param[out] Ne Taille de bloc d'entrée
 *
 *  @sa ola_complexite()
 *
 */
inline void ola_complexity_optimize(int M, float &C, int &Nf, int &Nz, int &Ne)
{
  tsdF::ola_complexite_optimise(M, C, Nf, Nz, Ne);
}


/** @brief Transformée en z-chirp
 *
 *  <h3>Transformée en z-chirp</h3>
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
inline ArrayXcf czt(IArrayXcf x, int m, cfloat W, cfloat z0 = 1.0f)
{
  return tsdF::czt(x, m, W, z0);
}

/** @} */

/** @addtogroup fourier-corr
 *  @{
 */

/** @brief Corrélation circulaire (normalisée) entre deux signaux complexes
 *
 *  <h3>Produit de corrélation circulaire</h3>
 *
 *  Calcul de la corrélation circulaire (normalisée) entre deux signaux complexes (calcul efficace via des FFT) :
 * @f[
 * c_n = \frac{1}{N} \cdot \sum_{k=0}^{N-1} x_k y^{\star}_{k+n[N]},\ \ n=0\dots N-1
 * @f]
 *
 *  @param x Premier vecteur
 *  @param y Deuxième vecteur (éventuellement égal au premier pour avoir l'auto-corrélation)
 *  @returns Premiers vecteur : index des retards @f$n@f$ (soit @f$0, 1, ..., N-1@f$), deuxième vecteur : @f$c_n@f$.
 *
 *
 * @sa xcorr()
 */
inline std::tuple<ArrayXf, ArrayXcf> ccorr(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf())
{
  return tsdF::ccorr(x, y);
}


/** @brief Corrélation (non biaisée) entre deux signaux complexes
 *
 *  <h3>Produit de corrélation (avec correction de biais)</h3>
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
 * @image html xcorr.png width=800px
 *
 * @sa xcorrb(), ccorr(), fft_correlateur()
 */
inline std::tuple<ArrayXf, ArrayXcf> xcorr(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf(), int m = -1)
{
  return tsdF::xcorr(x, y, m);
}


/** @brief Corrélation (biaisée) entre deux signaux complexes.
 *
 *  <h3>Produit de corrélation (sans correction de biais)</h3>
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
 *  @sa xcorr(), ccorr()
 */
inline std::tuple<ArrayXf, ArrayXcf> xcorrb(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf(), int m = -1)
{
  return tsdF::xcorrb(x, y, m);
}

/** @brief Délais entier (décalage du signal) ou fractionnaire (basé sur la FFT).
 *
 *  <h3>Délais</h3>
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
 *  @param x  Signal à retarder ou avancer dans le temps.
 *  @param τ  Délais à appliquer, en nombre d'échantillons (peut être positif ou négatif, et n'est pas forcément un nombre entier).
 *  @returns  Signal décalé dans le temps.
 *
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_delais
 * @image html ex-delais.png width=800px
 *
 **/
template<typename T = float>
  Vector<T> delay(const Vector<T> &x, float τ)
{
  return tsdF::delais(x, τ);
}

/** @brief Estimation du délais (à l'échantillon près) le
 *  plus probable entre deux signaux (via une corrélation) */
inline int delay_estimation_int(IArrayXcf x, IArrayXcf y, float &score)
{
  return tsdF::estimation_delais_entier(x, y, score);
}


/** @brief Alignement de deux signaux */
template<typename T>
  std::tuple<Vector<T>, Vector<T>, int, float> align_int(const Vector<T> &x, const Vector<T> &y)
{
  return tsdF::aligne_entier(x, y);
}


/** @brief Informations computed from the detected pattern. */
struct Detection
{
  /** @brief en nombre d'échantillons (compris entre 0 et Ne-1),
  *    depuis le début du bloc de données en cours. Par exemple :
  *    - @f$0\  \Leftrightarrow@f$    Début du bloc en cours
  *    - @f$1\  \Leftrightarrow@f$    Deuxième échantillon du bloc en cours
  *    - @f$-1\  \Leftrightarrow@f$   Dernier échantillon du bloc précédent
  *    - ... */
  int position;

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

  Detection(const tsdF::Detection &fr)
  {
    memcpy(this, &fr, sizeof(tsdF::Detection));
  }
};


/** @brief Structure de configuration pour un corrélateur par FFT */
struct DetectorConfig: tsdF::DetecteurConfig
{
  /* @brief Nombre minimum d'échantillons / par bloc temporel (si 0 : déterminé automatiquement). */
  uint32_t &Ns        = Ne;

  /* @brief Motif à détecter */
  ArrayXcf &pattern   = motif;

  /* @brief Seuil de détection, entre 0 et 1. */
  float &threshold    = seuil;

  bool &debug_active  = debug_actif;

  /** @brief Mode de calcul du produit de corrélation : FFT avec OLA ou simple filtre RIF */
  tsdF::DetecteurConfig::Mode &mode = mode;

  /** @brief Callback utilisateur appellée à chaque fois que le motif est détecté. */
  std::function<void (const Detection &det)> &on_detection
      = *((std::function<void (const Detection &det)> *) &gere_detection);

  bool &compute_correlation_signal = calculer_signal_correlation;
};



using tsdF::Detecteur;

/** @brief Détecteur par corrélation.
 *
 * <h3>Détecteur par corrélation</h3>
 *
 * Cette fonction calcule, au fil de l'eau, la corrélation normalisée entre un signal @f$(x_k)@f$
 * reçu en continu et
 * un motif fixe @f$(h_k)@f$ de dimension @f$M@f$ :
 * @f[
 * y_n = \frac{\left|\displaystyle{\sum_{k=0}^{M-1} x_{n+k-(M-1)} \cdot h_{k}^\star}\right|}{\displaystyle{ \sqrt{\left(\sum_{k=0}^{M-1}\left|x_{n+k-(M-1)}\right|^2\right) \left(\sum_{k=0}^{M-1}\left|h_k\right|^2\right)}}}
 * @f]
 *
 * Les @f$y_n@f$ sont donc compris entre 0 et 1, et seront égaux à 1 uniquement si le signal est
 * exactement égal au motif (à un facteur d'échelle constant près) sur l'intervalle @f$[n-M+1\dots n]@f$.
 *
 * Le calcul étant fait de manière efficace dans le domaine fréquentiel via la technique Overlap-And-Add
 * (voir @ref filter_fft()), la complexité est de l'ordre de @f$\log_2 M@f$ opérations
 * par échantillon (si @f$M@f$ est une puissance de 2).
 *
 * Ce filtre peut-être utilisé de deux manières :
 *  - Tout d'abord, c'est un filtre qui renvoie, au fil de l'eau, les valeurs de corrélation (signal vert dans l'exemple ci-dessous),
 *  - Pour une utilisation plus simple, une callback utilisateur peut être appelée dès lors que la corrélation dépasse un certain seuil paramétrable. En paramètre de cette callback,
 *    est indiquée la position du motif détecté, ainsi que diverses informations (voir la structure @ref Detection) :
 *      - Valeur de la corrélation normalisée
 *      - Gain réel du signal (rapport d'amplitude entre le signal reçu et le motif théorique attendu)
 *      - Déphasage du signal reçu
 *      - Ecart-type du bruit
 *      - SNR estimé
 *
 * @param config Structure de configuration
 * @return Filtre cfloat @f$\to@f$ float, qui prend en entrée des échantillons, et sort, au fil de l'eau,
 * la corrélation normalisée avec le motif.
 *
 * @par Exemple
 * @snippet exemples/src/fourier/ex-fourier.cc ex_fft_correlateur
 * @image html fft_correlateur.png width=800px
 */
inline sptr<Detecteur>
  detector_new(const DetectorConfig &config = DetectorConfig())
{
  return tsdF::détecteur_création(config);
}



/** @} */



/** @addtogroup fourier-freq
 *  @{
 */

// Calcul des fréquences normalisées associée à une psd
inline ArrayXf psd_freqs(int n, bool complexe = true)
{
  return tsdF::psd_freqs(n, complexe);
}

/** @brief PSD (corrélogramme)
 *
 * <h3>PSD (corrélogramme)</h3>
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
 * @image html ex-psd2.png "PSD signal triangulaire" width=800px
 *
 * @sa psd_welch(), psd_subspace()
 */
template<typename derived>
std::tuple<ArrayXf, ArrayXf> psd(const Eigen::ArrayBase<derived> &x)
{
  return tsdF::psd(x);
}

/** @brief PSD (méthode de Welch - par moyennage).
 *
 * <h3> PSD par la méthode de Welch</h3>
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
 * @image html ex-psd3.png "PSD et moyennage" width=800px
 *
 * @sa psd(), psd_subspace()
 *
 */
inline std::tuple<ArrayXf, ArrayXf> psd_welch(const ArrayXcf &x, int N, const std::string fen = "hn")
{
  return tsdF::psd_welch(x, N, fen);
}


/** @brief Calcul d'un spectre par la méthode des sous-espaces
 *
 *  <h3>Spectre (méthode des sous-espace / MUSIC)</h3>
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
 *  (voir la fonction @ref psd_subspace()).
 *
 *  @warning Cette technique implique la décomposition en valeurs propres d'une matrice @f$m\times m@f$
 *  (la matrice d'auto-corrélation du signal), soit de l'ordre de @f$m^3@f$ opérations.
 *  Par conséquent, si la dimension du signal est importante,
 *  il est conseillé de choisir un @f$m@f$ plus petit.
 *
 *  @par Exemple : détection de 3 exponentielles pures
 *  @snippet exemples/src/fourier/ex-fourier.cc ex_psd_subspace
 *  @image html ex-subspace-freq-spectrum.png "Spectre PSD vs MUSIC" width=800px
 *  Notez comme le spectre MUSIC est plus propre ; il ne faut cependant pas oublier
 *  que cette technique suppose connu le nombre de signaux (ici 3).
 *
 *  @sa subspace_spectrum()
 *
 *  @note Pour la détection de sinusoïdes pures (signaux réels), il faudra choisir @f$N_s@f$ comme étant <b>deux fois</b>
 *  le nombre de signaux à détecter (en effet, un cosinus ou sinus peut se représenter comme la somme de deux exponentielles complexes).
 *
 */
inline std::tuple<ArrayXf, ArrayXf> psd_subspace(const ArrayXcf &x, int Ns, int Nf = 1024, int m = 0)
{
  return tsdF::psd_subspace(x, Ns, Nf, m);
}


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
 *  <h3>Estimation de fréquence</h3>
 *
 *  Cette fonction essaye d'estimer la fréquence d'un signal périodique (exponotielle pure).
 *
 *  @note Si le signal d'entrée est réel (sinusoide pure de fréquence @f$f@f$),
 *  la fonction renvoie une des deux valeurs possibles de fréquence (@f$f@f$ ou @f$-f@f$).
 *
 *  @warning Si le signal d'entrée est composé de plusieurs fréquences distinctes, cette fonction
 *  n'est pas adaptée. Utilisez plutôt un algorithme de type sous-espace (voir @ref psd_subspace()).
 *
 *
 *  @param x Signal dont on veut estimer la fréquence
 *  @param m Méthode d'estimation (voir l'énumération @ref FreqEstimMethode)
 *  @return Une fréquence normalisée (entre -0,5 et 0,5) */
inline float freqestim(IArrayXcf x, FreqEstimMethode m = FreqEstimMethode::CANDAN2)
{
  return tsdF::freqestim(x, (tsdF::FreqEstimMethode) m);
}




/** @brief %Filtre de Goertzel (calcul sur un buffer).
 *
 *  <h3>%Filtre de Goertzel (calcul sur un buffer)</h3>
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
inline float goertzel(const ArrayXf &x, float frequence)
{
  return tsdF::goertzel(x, frequence);
}

/** @brief %Filtre de Goertzel (calcul au fil de l'eau)
 *
 *  <h3>%Filtre de Goertzel (calcul au fil de l'eau)</h3>
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
 *  @image html goertzel.png width=600px
 *
 *  @sa goertzel()
 */
inline sptr<FiltreGen<float>> filter_goertzel(float frequence, int N)
{
  return tsdF::filtre_goertzel(frequence, N);
}




/** @brief Configuration pour le calcul de spectre en temps réel */
using tsdF::SpectrumConfig;



/** @brief Calcul de spectre en temps réel */
inline sptr<Filtre<cfloat,float,SpectrumConfig>> rt_spectrum(const SpectrumConfig &config)
{
  return tsdF::rt_spectrum(config);
}

/** @} */

} // namespace tsd::fourier

namespace tsd::tf {

/** @addtogroup temps-frequence
*  @{ */

/** @brief Calcul d'une matrice temps / fréquence par FFT glissante
 *
 *  <h3>Spectrogramme par FFT glissante</h3>
 *
 *  @param x Signal à analyser
 *  @param N Dimension des blocs d'analyse */
inline ArrayXXf periodogram_dft(const ArrayXcf &x, int N)
{
  return tsd::tf::periodogramme_tfd(x, N);
}

/** @brief Calcul d'un spectrogramme à partir de la CQT
 *
 *  <h3>Spectrogramme CQT (Constant Quality Transform)</h3>
 *
 *  @param x Signal à analyser
 *  @param fe Fréquence d'échantillonnage (Hz)
 *  @param fmin Fréquence minimale d'analyse (Hz)
 *  @param fmax Fréquence minimale d'analyse (Hz)
 *  @param γ    Rapport entre deux fréquences d'analyse successives
 *  @param ofs  Fréquence d'échantillonnage de sortie (Hz)
 */
inline ArrayXXf periodogram_cqt(const ArrayXf &x, float fe, float fmin, float fmax, float γ, float ofs)
{
  return tsd::tf::periodogramme_cqt(x, fe, fmin, fmax, γ, ofs);
}

/* @} */

}



