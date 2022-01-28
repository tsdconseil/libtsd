#pragma once

#include "dsp/dsp.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/filtrage/frat.hpp"
#include "dsp/figure.hpp"

namespace dsp {
/*namespace view {
  struct Figure;
  struct Figures;
}*/

  template<typename T>
   using FRat = tsd::FRat<T>;
}


namespace dsp::filter {

  namespace tsdf = tsd::filtrage;

  /** @addtogroup filtrage-fenetres
    *  @{ */



  /** @cond private  */
  using Window = tsdf::Fenetre;

  inline ArrayXf window(Window type, int n, bool symetrical = true)
  {
    return tsdf::fenetre(type, n, symetrical);
  }
  /** @endcond */

  /** @brief Parameterless window creation (rectangular, Hann, Hamming, triangular or Blackman)
   *
   *  <h3>Parameterless window</h3>
   *
   *  This function will create a simple window, defined without parameter (apart from its dimension and whether it should be symetrical or not).
   *
   *  @param type Window choice: "re" (rectangular), "hn" (Hann), "hm" (Hamming), "tr" (triangular), "bm" (Blackman).
   *  @param n Number of points.
   *  @param symetrical If true, realization of window symetrical around its center point (which should be done for the FIR filter design),
   *  otherwise realization of a periodical window (which is better for spectral analysis applications).
   *  @return Column vector with the window coefficients (dimension of vector is @p n)
   *  @sa window_kaiser(), window_chebychev()
   *
   *  @par Example: creation of a Von Hann window (also called Hanning)
   *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_fenetre
   *  @image html filtrage-fenetre.png width=600px
   */
  inline ArrayXf window(const std::string &type, int n, bool symetrical = true)
  {
    return tsdf::fenetre(type, n, symetrical);
  }

  /** @brief Chebychev window creation.
   *
   *  <h3>Chebychev window creation</h3>
   *
   *  The Chebychev window has the property to give rise to constant ondulations in the frequency domain (si example below).
   *  The design is also very straigthforward for FIR filter design, because one has just to select the window order (number of coefficients),
   *  and the desired attenuation,
   *  and it is the main lobe width that play the role of adjustment variable.
   *
   *  @param n            Number of points.
   *  @param atten_db     Desired attenuation (in dB) on the stop-band.
   *  @param symetrical   If true, realization of window symetrical around its center point (which should be done for the FIR filter design),
   *  otherwise realization of a periodical window (which is better for spectral analysis applications).
   *  @return Column vector with the window coefficients (dimension of vector is @p n)
   *  @par Example: creation of a window with 60 dB of attenuation
   *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_fenetre_cheby
   *  @image html filtrage-fenetre-cheby.png width=600px
   *
   *  @sa design_rif_fen()
   */
  inline ArrayXf window_chebychev(int n, float atten_db, bool symetrical = true)
  {
    return tsdf::fenêtre_chebychev(n, atten_db, symetrical);
  }


  /** @brief Computation of the parameters of a Kaiser filter.
   *
   *  <h3>Computation of the parameters of a Kaiser filter</h3>
   *
   *  This function will compute, for a Kaiser filter, the @f$\beta@f$ parameter
   *  and the required order for a specified attenuation and transition bandwidth.
   *
   *  @param atten_db Attenuation in dB in the stop-band (should be a positive value).
   *  @param df Transition bandwidth (normalized to the sampling frequency)
   *  @return @f$\beta,n@f$
   *
   *  @f[
   *    n = (A - 7.95) / (2.285 \cdot 2 \pi d_f)
   *  @f]
   *
   *  @f[
   *   \beta = \begin{cases}
   *   0.1102 (A - 8.7) & \mbox{si } A > 50\textrm{ dB} \\
   *   0.5842 \cdot (A - 21)^{0.4} + 0.07886 \cdot (A - 21) & \mbox{si } 21 \textrm{ dB} \leq A \leq 50\textrm{ dB}\\
   *   0 & \mbox{sinon}
   *   \end{cases}
   *  @f]
   *  @sa window_kaiser(), window_kaiser1()
   *
   *  @par Example
   *  @code
   *  // Parameters for 60 dB of transition,
   *  // and a transition bandwidth of  one tenth of the sampling frequency.
   *  auto [β, n] = kaiser_param(60, 0.1);
   *  @endcode
   */
  inline std::tuple<float, int> kaiser_param(float atten_db, float df)
  {
    return tsdf::kaiser_param(atten_db, df);
  }

  /** @brief Creation of a Kaiser window.
   *
   *  <h3>Creation of a Kaiser window</h3>
   *
   *  @param atten_db Attenuation in dB in the stop-band (should be a positive value).
   *  @param df       Transition bandwidth (normalized to the sampling frequency)
   *  @param symetrical If true, realization of window symetrical around its center point (which should be done for the FIR filter design),
   *  otherwise realization of a periodical window (which is better for spectral analysis applications).
   *
   *  @return Column vector with window coefficients.
   *  @sa fenetre_kaiser1(), kaiser_param()
   *
   *  @par Example: creation of a window with 60 dB of attenuation
   *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_fenetre_kaiser
   *  @image html filtrage-fenetre-kaiser.png width=600px
   *   */
  inline ArrayXf window_kaiser(float atten_db, float df, bool symetrical = true)
  {
    return tsdf::fenêtre_kaiser(atten_db, df, symetrical);
  }

  /** @brief Creation of a Kaiser window (from the shape parameter @f$\beta@f$).
   *
   * <h3>Creation of a Kaiser window</h3>
   *
   * This fonction will compute a Kaiser window with a specified shape parameter @f$\beta@f$.
   *
   * @param n     Number of points.
   * @param β     Shape parameter.
   * @param symetrical If true, realization of window symetrical around its center point (which should be done for the FIR filter design),
   *  otherwise realization of a periodical window (which is better for spectral analysis applications).
   *
   * @return Column vector with window coefficients.
   *
   *  @sa fenetre_kaiser(), kaiser_param() */
  inline ArrayXf window_kaiser1(int n, float β, bool symetrical = true)
  {
    return tsdf::fenêtre_kaiser1(n, β, symetrical);
  }



  /** @} */



  /** @addtogroup filtrage-analyse
   *  @{ */

/** @brief Magnitude of a FIR or IIR filter.
 *
 *  <h3>Magnitude of a FIR or IIR filter</h3>
 *
 *  This fonction computes the magnitude (e.g., the absolute value) of the frequency response for a linear filter:
 *  @f[
 *  y_k = \left|H(f_k)\right|
 *  @f]
 *
 *  @param h Transfert function to be analyzed.
 *  @param npts Frequency resolution.
 *  @return A tuple of 2 vectors : the frequencies @f$f_k@f$ (normalized, between 0 and 0.5), and the magnitudes @f$y_k@f$.
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage-en.cc ex_frmag_en
 *  @image html frmag.png width=600px
 *
 *  @sa repfreq(), frgroup(), frphase()
 */
template<typename T>
  std::tuple<ArrayXf, ArrayXf> frmag(const FRat<T> &h, unsigned int npts = 1024)
{
  return tsdf::frmag(h, npts);
}

/** @brief Frequency response of a FIR or IIR filter.
 *
 * <h3>Frequency response of a FIR or IIR filter</h3>
 *
 *  This function computes the (complex) frequency response:
 *  @f[
 *  y_k = H(f_k)
 *  @f]
 *
 *  @param h   Transfert function to be analyzed.
 *  @param fr  Frequencies vector (normalized to sampling frequency).
 *  @return    Frequency response @f$H@f$
 *
 *  @sa frmag(), frphase(), frgroup()
 */
template<typename T>
  ArrayXcf repfreq(const FRat<T> &h, const ArrayXf &fr)
{
  return tsdf::repfreq(fr, fr);
}

inline ArrayXcf repfreq(const ArrayXf &h, const ArrayXf &fr)
{
  return tsdf::repfreq(h, fr);
}

/** @brief Phase of FIR and IIR filters.
 *
 *  <h3>Phase of FIR and IIR filters</h3>
 *
 *  This function computes the phase response of a FIR or IIR filter:
 *  @f[
 *  y_k = \arg H(f_k)
 *  @f]
 *
 *  @param h Transfert function to be analyzed.
 *  @param npts Frequency resolution.
 *  @return A tuple of 2 vectors : the frequencies @f$f_k@f$ (normalized, between 0 and 0.5), and the phases @f$y_k@f (in radians)$.
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frphase(const FRat<T> &h, unsigned int npts = 1024)
{
  return tsdf::frphase<T>(h, npts);
}

/** @brief Calcul du temps de groupe (délais en fonction de la fréquence).
 *
 *  <h3>Temps de groupe</h3>
 *
 *  Calcul du délais du filtre en fonction de la fréquence :
 *  @f[
 *  G(\omega) = \frac{d\arg H(\omega)}{d\omega}
 *  @f]
 *
 *  @param h Fonction de transfert à analyser
 *  @param npts Résolution fréquentielle
 *  @return Un tuple de deux vecteurs : le vecteur de fréquences (normalisées, entre 0 et 0,5),
 *  et le temps de groupe (en nombre d'échantillons).
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_frgroup
 *  @image html frgroup.png width=600px
 *
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frgroup(const FRat<T> &h, unsigned int npts = 1024)
{
  return tsdf::frgroup<T>(h, npts);
}




/** @brief Analyse d'un filtre linéaire (tracé des différentes réponses).
 *
 *  <h3>Analyse d'un filtre linéaire</h3>
 *
 * Cette fonction créée une nouvelle figure, et trace les réponses fréquentielles
 * et temporelles du filtre,
 * ainsi que le diagramme des zéros et des pôles.
 *
 *  @param h      Fonction de transfert à analyser (ou vecteur de coefficients).
 *  @param fe     Fréquence d'échantillonnage (optionnel).
 *  @return       La figure créée.
 *
 * @par Exemple
 * @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_analyse
 * @image html filtrage-analyse.png width=1000px
 *
 * */
template<typename T>
  sptr<dsp::view::Figures> filter_analysis(const FRat<T> &h, float fs = 1.0f)
{
  return tsdf::analyse_filtre(h, fs);
}


/** @cond undoc */
template<typename T> std::tuple<ArrayXf, ArrayXf> frmag(const Vector<T> &h, unsigned int npts = 1024)
{
  return tsdf::frmag<T>(h, npts);
}
template<typename T> std::tuple<ArrayXf, ArrayXf> frphase(const Vector<T> &h, unsigned int npts = 1024)
{
  return tsdf::frphase<T>(h, npts);
}
template<typename T> std::tuple<ArrayXf, ArrayXf> frgroup(const Vector<T> &h, unsigned int npts = 1024)
{
  return tsdf::frgroup<T>(h, npts);
}
inline dsp::view::Figures analyse_filtre(const ArrayXf &h, float fe = 1.0f)
{
  return tsdf::analyse_filtre(h, fe);
}
/** @endcond */





/** @brief Tracé des pôles et zéros.
 *
 *  <h3>Tracé des pôles et des zéros</h3>
 *
 *  La fonction de transfert passée en paramètre est factorisée sous la forme :
 *  @f[
 *  H(z) = \frac{\prod z - z_i}{\prod z - p_i}
 *  @f]
 *
 *  où les @f$z_i@f$ et les @f$p_i@f$ sont ce que l'on appelle respectivement
 *  les zéros et les pôles de la fonctions de transfert.
 *
 *  @param h Fonction de transfert
 *  @param fig Figure sur laquelle sera tracé le diagramme des pôles et zéros.
 *
 *
 *  @par Exemple
 *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_plz
 *  @image html filtrage-plz.png width=600px
 *  */
template<typename T>
  void plot_plz(dsp::view::Figure &fig, const FRat<T> &h);


/**  @}
  *  @addtogroup filtrage-design
  *  @{ */



struct SpecFreqIntervalle
{
  float fb, fh;
  float atten = 1.0f;
  float poids = 1.0f;
};

/** @brief Approximation RIF d'un filtre de Hilbert.
 *
 *  <h3>%Filtre de Hilbert (approximation RIF)</h3>
 *
 *  Cette fonction calcule le fenêtrage de la réponse temporelle théorique d'un filtre de Hilbert :
 *  @f[
 *      h_k = \frac{2}{k\pi} \cdot \sin(k \pi / 2)^2 \cdot w_k;
 *  @f]
 *
 *  @param n ordre du filtre
 *  @param fenetre Type de fenetre (par défaut, fenêtre de Hann)
 *  @returns Tableau des coefficients
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_hilbert
 *  @image html design-rif-hilbert.png width=600px
 *
 *  @sa hilbert(), hilbert_transformeur()
 */
inline ArrayXf design_fir_hilbert(int n, const std::string &window = "hn")
{
  return tsdf::design_rif_hilbert(n, window);
}

/** @brief Design RII d'après un prototype analogique classique.
 *
 * <h3>Design RII d'après un prototype analogique classique</h3>
 *
 * Cette fonction renvoie une fonction de transfert <b>discrète</b>, sous forme de pôles et de zéros (idéal pour une implémentation sous la forme de sections du second ordre, voir @ref filtre_sois()).
 * Les prototypes supportés sont les suivants :
 *  - <b>Butterworth</b> (pas d'ondulation, bande de transition large)
 *  - <b>Chebychev type I</b> (ondulations dans la bande passante)
 *  - <b>Chebychev type II</b> (ondulations dans la bande coupée)
 *  - <b>Elliptique</b> (ondulations partout, mais bande de transition la plus étroite)
 *
 *  Le filtre est d'abort conçu dans le domaine analogique, puis converti en filtre digital grâce à la
 *  transformée bilinéaire.
 *
 * @param n           Ordre du filtre
 * @param type        Type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...)
 * @param prototype   "butt", "cheb1", "cheb2" ou "ellip"
 * @param fc          Fréquence de coupure normalisée (entre 0 et 0,5)
 * @param δ_bp        Ondulation maximale en decibels dans la bande passante (utilisé seulement pour un filtre de Chebychev de type I ou un filtre elliptique).
 * @param δ_bc        Atténuation minimale en decibels dans la bande coupée (utilisé seulement pour un filtre de Chebychev type II ou un filtre elliptique).
 * @return h          Fonction de transfert (digitale)
 *
 * @par Exemple
 * @snippet exemples/src/ex-filtrage.cc ex_design_riia
 * @image html design-riia.png width=800px
 *
 * @sa trf_bilineaire(), @ref filtre_sois(), filtre_rii()
 */
inline FRat<cfloat> design_iira(int n, const std::string &type,
    const std::string &prototype, float fc, float δ_bp = 0.1f, float δ_bc = 60)
{
  return tsdf::design_riia(n, type, prototype, fc, δ_bp, δ_bc);
}



/** @brief Type de filtre biquad (voir @ref design_biquad()) */
using BiquadSpec = tsdf::BiquadSpec;

/** @brief Design filtre biquad.
 *
 * D'après https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
 *
 * @param Q Facteur de qualité
 *
 *  */
inline FRat<float> design_biquad(const BiquadSpec &spec)
{
  return tsdf::design_biquad(spec);
}

/** @cond private
 */
enum class PrototypeAnalogique
{
  BUTTERWORTH = 0,
  TCHEBYCHEV_I,
  TCHEBYCHEV_II,
  ELLIPTIQUE
};

enum class TypeFiltre
{
  PASSE_BAS = 0,
  PASSE_HAUT,
  PASSE_BANDE,
  COUPE_BANDE
};

extern FRat<cfloat> design_riia(unsigned int n, TypeFiltre type,
    PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc);

extern FRat<cfloat> design_riia_laplace(int n, TypeFiltre type, PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc);

/** @endcond */


/** @brief Design par échantillonnage fréquentiel.
 *
 *  @param n  Ordre du filtre
 *  @param d  Vector définissant la réponse fréquentielle souhaitée (sur les fréquences positives)
 *  @returns  Vector des coefficients du filtre (dimension = n)
 */
inline ArrayXf design_fir_freq(int n, const ArrayXf &d)
{
  return tsdf::design_rif_freq(n, d);
}

/** @brief Design équiripple / Remez.
 *
 *  <h3>Design équiripple / Remez</h3>
 *
 *  @param n      Ordre du filtre
 *  @param d      Réponse souhaitée
 *  @param w      Poids de pondération (doit être de même longueur que d)
 *  @param debug  Si vrai, affiche des plots des calculs intermédiaires.
 *  @returns      Vector des coefficients du filtre (vecteur de dimension n)
 */
inline ArrayXf design_fir_eq(unsigned int n, IArrayXf d, IArrayXf w, bool debug = false)
{
  return tsdf::design_rif_eq(n, d, w, debug);
}

extern ArrayXf design_rif_eq(unsigned int n, const std::vector<SpecFreqIntervalle> &spec, bool debug = false);

/** @brief Sinus cardinal normalisé avec fréquence de coupure paramétrable.
 *
 *  <h3>Sinus cardinal normalisé avec fréquence de coupure paramétrable</h3>
 *
 *  @f[
 *  y(t) = \frac{\sin(2 \pi t f_c)}{\pi  t}
 *  @f]
 *
 *  La transformée de Fourier de cette fonction est une porte de largeur @f$\pm f_c@f$, ce qui en fait donc le protype idéal pour un filtre passe-bas.
 *
 *  @param t Point d'échantillonnage (en nombre d'échantillons)
 *  @param fc Fréquence de coupure
 *  @returns Valeur de @f$y@f$
 *
 *  @sa sinc(), design_rif_fen() */
inline float sinc2(float t, float fc)
{
  return tsdf::sinc(t);
}

/** @brief Sinus cardinal normalisé avec fréquence de coupure 0,5.
 *
 *  <h3>Sinus cardinal normalisé</h3>
 *  Calcul d'un sinus cardinal (fréquence de coupure 0,5) :
 *  @f[
 *  y(t) = \frac{\sin(\pi t)}{\pi  t}
 *  @f]
 *
 *  @param t Point d'échantillonnage (en nombre d'échantillons)
 *  @returns Valeur de @f$y@f$
 *
 *  @sa sinc2()
 *
 *  */
inline float sinc(float t)
{
  return tsdf::sinc(t);
}

/** \~french  @brief Design RIF par sinus-cardinal fenêtré.
 *
 * @param n       ordre du filtre
 * @param type    type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...)
 * @param fc      fréquence de coupure normalisée
 * @param fen     type de fenêtre ("hn", "hm", "tr" ou "re"), voir fonction @ref fenetre()
 * @param fc2     deuxième fréquence de coupure (uniquement pour les filtres passe-bande ou stoppe-bande)
 * @returns vecteur des coefficients du filtre
 * @sa design_rif_eq(), design_rif_freq()
 *
 *  \~english @brief Windowed sinc filter */
extern ArrayXf design_rif_fen(unsigned int n, const std::string &type, float fc, const std::string &fen = "hn", float fc2 = 0);

/** @brief Design RIF par sinus-cardinal fenêtré (fenêtre de Kaiser).
 *
 *  <h3>Design RIF fenétré (Kaiser)</h3>
 *
 *  L'utilisation d'une fenêtre de Kaiser permet de choisir à la fois l'atténuation du filtre
 *  et la largeur de la bande de transition, la variable d'ajustement étant le nombre de coefficients.
 *
 * @param type    type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...)
 * @param fc      fréquence de coupure normalisée
 * @param atten_db Atténuation souhaitée dans la bande coupée (en décibels)
 * @param df      Largeur de la bande de transition (en fréquence normalisée)
 * @param fc2     deuxième fréquence de coupure (uniquement pour les filtres passe-bande ou stoppe-bande)
 *
 * @sa design_rif_fen(), design_rif_fen_chebychev()
 */
extern ArrayXf design_rif_fen_kaiser(const std::string &type, float fc, float atten_db,
    float df, float fc2 = 0);

/** @brief Design RIF par sinus-cardinal fenêtré (fenêtre de Chebychev).
 *
 *  <h3>Design RIF fenétré (Chebychev)</h3>
 *
 *  L'utilisation d'une fenêtre de Chebychev permet de choisir à la fois l'atténuation du filtre
 *  et le nombre de coefficients, la variable d'ajustement étant la bande de transition.
 *
 * @param n       ordre du filtre
 * @param type    type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...)
 * @param fc      fréquence de coupure normalisée
 * @param atten_db Atténuation souhaitée dans la bande coupée (en décibels)
 * @param fc2     deuxième fréquence de coupure (uniquement pour les filtres passe-bande ou stoppe-bande)
 *
 * @sa design_rif_fen(), design_rif_fen_kaiser()
 */
extern ArrayXf design_rif_fen_chebychev(int n, const std::string &type,
    float fc, float atten_db, float fc2 = 0);

/** @brief Design d'un filtre en cosinus sur-élevé.
 *
 *  <h3>Design d'un filtre en cosinus sur-élevé</h3>
 *
 *  @param n  ordre du filtre
 *  @param α  facteur de dépassement
 *  @param fc fréquence de coupure normalisée
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_cs
 *  @image html design-rif-cs.png width=800px
 *
 *  @sa design_rif_rcs() */
extern ArrayXf design_rif_cs(int n, float α, float fc);

/** @brief Design d'un filtre en racine de cosinus sur-élevé.
 *
 *  <h3>Design d'un filtre en racine de cosinus sur-élevé</h3>
 *
 *  @param n  ordre du filtre
 *  @param α  facteur de dépassement
 *  @param fc fréquence de coupure normalisée
 *  @sa design_rif_cs()
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_rcs
 *  @image html design-rif-rcs.png width=800px
 *
 */
extern ArrayXf design_rif_rcs(int n, float α, float fc);

extern ArrayXf design_rif_rcs1(int n, float α, float osf, char nrm = 's');

/** @brief Coefficients d'une approximation RIF d'un filtre gaussien (non fenêtré), pour une modulation GFSK.
 *
 * <h3>%Filtre Gaussien (approximation RIF)</h3>
 *
 *  Les paramètres de cette fonction sont adaptés à la réalisation
 *  d'un filtre de mise en forme pour une modulation GFSK.
 *
 *  @param n   Nombre de coefficients du filtre
 *  @param BT  Produit Bande-passante - Temps
 *  @param osf Facteur de sur-échantillonnage (rapport entre la fréquence d'échantillonnage et la fréquence symbole)
 *  @return Vector des coefficients du filtre
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_gaussien
 *  @image html design-rif-gaussien.png width=600px
 *
 *
 */
extern ArrayXf design_rif_gaussien(unsigned int n, float BT, unsigned int osf);


/** @~french  @brief Fonction de transfert théorique d'un filtre CIC
    @~english @brief CIC theorical transfert function
    @~

    <h3>Fonction de transfert théorique d'un filtre CIC</h3>

    @param R Decimation ratio
    @param N Number of integrators / differentiators stages
    @param M Optionnal CIC design parameter (default is M=1)
    @returns CIC transfert function (not taking into account the decimation)
    This function computes the theorical transfert function of a CIC filter, when one does not look at the decimation effect. The CIC filter responses is defined as:
           @f[H(z) = \frac{1}{R^N}\left(1+z^{-1}+\cdots+z^{-(R-1)}\right)^N@f]
     (e.g. a cascade of @f$N@f$ moving average filters, each of identical length @f$R@f$.
 **/
extern FRat<float> design_cic(int R, int N, int M = 1);

/** @brief Design of a compensation FIR filter for a CIC filter.
 *
 * <h3>Design of a compensation FIR filter for a CIC filter</h3>
 *
 * Calling Sequence
 * cfir = cic_comp_design(R,N,M,Fin,R2,Fcut,ntaps)
 *
 * @param R       CIC decimation ratio
 * @param N       number of CIC integrators / differentiators
 * @param M       CIC design parameter (1 ou 2)
 * @param Fin     input frequency
 * @param R2      compensation filter decimation ratio
 * @param fc      cutoff frequency for the compensation FIR filter
 * @param ncoefs  number of taps of the FIR filter
 * @returns       compensation FIR filter coefficients
 *
 * This function will try to design a compensation filter for the specified CIC filter. The global decimation chain will be composed of two stages:
 *   - The first stage is the CIC filter, and decimate by the ratio @f$R@f$.
 *   - The second stage is the compensation filter, and decimate by the ratio @f$R_2@f$.
 *
 * So, the global decimation ratio will be @f$R\cdot R_2@f$.
 *
 * The compensation filter is generated through the frequency sampling technique.
 * This function will do the FIR compensation design, plot the frequency responses, and
 * output (output parameter @p cfir) the compensation filter coefficients.
 * @todo déplacer la code dans exemple cic
 * @code
    // Input signal frequency: Fin = 6.4 KHz
    // CIC decimation factor: R = 1/16 (freq at output of CIC: 400 Hz)
    // Compensation FIR decimation factor: R2 = 1/2
    // Output signal frequency: 6400/(R*R2) = 200 Hz
    // Cut-off frequency for FIR filter: 80 Hz
    int R = 16, N = 4, M = 1;
    float fin = 6400, fcut=80;
    int ntaps = 61, R2 = 2;
    ArrayXf h = design_cic_comp(R,N,M,fin,R2,fcut,ntaps);
    // h is the array of coefficents of the FIR compensation filter.
 * @endcode
 * @todo ex_cic_comp_rimp.png
 *    Impulse response of the FIR compensation filter
 * @todo ex_cic_comp_glob.png
 *     Global response of CIC + compensation filters (spectrum) */
extern ArrayXf design_cic_comp(unsigned int R, unsigned int N, unsigned int M, float Fin, unsigned int R2, float fc, unsigned int ncoefs);


/** @brief Fonction de transfert pour un bloqueur DC.
 *
 *  <h3>Fonction de transfert bloqueur DC</h3>
 *
 *  Calcul de la fonction de transfert d'un bloqueur DC :
 *  @f[
 *  H(z) = \frac{1 - z^{-1}}{1 - \alpha \cdot z^{-1}}
 *  @f]
 *  Avec
 *  @f[
 *  \alpha = \frac{\sqrt{3} - 2 \sin(\pi f_c)}{\sin(\pi f_c) + \sqrt{3} \cos(\pi f_c)}
 *  @f]
 *
 *  d'après The DC Blocking Filter, J M de Freitas, 2007.
 *
 *  @param fc Fréquence de coupure normalisée (0 - 0,5)
 *  @sa filtre_dc()
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc exemple_design_bloqueur_dc
 *  @image html bloqueur-dc-resp.png width=1000px
 *
 *  @par Bibliographie
 *  @ref biblio "[RI1]"
 *
 *
 */
extern FRat<float> design_bloqueur_dc(float fc);

/** @brief Fonction de transfert d'un filtre exponentiel
 *
 *  <h3>%Filtre exponentiel</h3>
 *  Cette fonction renvoie la fonction de transfert d'un filtre exponentiel de fréquence de coupure
 *  @f$f_c@f$ :
 *  @f[
 *  y_n = \gamma x_n + (1-\gamma) y_{n-1}
 *  @f]
 *  Soit :
 *  @f[
 *  H(z) = \frac{\gamma z}{z-(1-\gamma)}
 *  @f]
 *  avec
 *  @f$\gamma = 1 - e^{-2\pi f_c}@f$.
 *
 *  @param fc Fréquence de coupure à -3 dB
 *  @return Fonction de transfert digitale
 *
 *  @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rii1
 *  @image html design-rii1.png width=1000px
 */
FRat<float> design_iir1(float fc){return tsdf::design_rii1(fc);}

/** @brief Generate first order IIR filter forget factor (unique coefficient) from cut-off frequency
 *
 *  <h3>Design filtre RII d'ordre 1</h3>
 *
 *  Cette fonction calcule le coefficient (facteur d'oubli) d'un filtre RII du premier ordre
 *  (filtre exponentiel), en fonction de la fréquence de coupure souhaitée.
 *
 *  @f[
 *    y_{n+1} = (1-\gamma) * y_n + \gamma * x_n
 *            = y_n + \gamma * (x_n - y_n)
 *  @f]
 *
 *  @f[
 *    \gamma = 1 - e^{-2\pi f_c}
 *  @f]
 *
 *  @param fc Normalized cut-off frequency, in [0, 0.5]
 *  @returns Coeficient @f$\gamma@f$ for the IIR filter
 *
 */
inline float iir1_coef(float fc){return tsdf::rii1_coef(fc);}

/** @brief Same as previous function, but take as input the time constant (in samples).
 *
 *  <h3>Design filtre RII d'ordre 1</h3>
 *
 *  Cette fonction calcule le coefficient (facteur d'oubli) d'un filtre RII du premier ordre
 *  (filtre exponentiel), en fonction de la constante de temps souhaitée.
 *
 *  @f[
 *    y_{n+1} = (1-\gamma) * y_n + \gamma * x_n
 *            = y_n + \gamma * (x_n - y_n)
 *  @f]
 *
 *  @f[
 *    \gamma = 1 - e^{-1/t_c}
 *  @f]
 *
 *  @param tc Constante de temps, en nombre de symboles.
 *  @returns Coeficient @f$\gamma@f$ for the IIR filter
 *
 */
extern float rii1_tc_vers_coef(float tc);

extern float rii1_coef_vers_tc(float γ);

/** @brief Compute cut-off frequency, from forget factor of first order IIR filter.
 *
 *  @param γ Facteur d'oubli
 *  @returns Fréquence de coupure
 *
 *  @sa rii1_coef */
extern float rii1_fcoupure(float γ);


/** @brief Creation of the polyphase representation of a signal
 *
 *  <h3>Représentation polyphase d'un signal</h3>
 *
 *  Creation of the polyphase matrix X, with zero padding if necessary (so as the length is a multiple of M):
 *
 *  @f[
 *  X = \left(\begin{array}{cccc}
 *  x_0 & x_M & x_{2M} & \dots\\
 *  x_1 & x_{M+1} & x_{2M+1} & \dots\\
 *  \vdots & \vdots & \vdots & \vdots\\
 *  x_{M-1} & x_{M+M-1} & x_{2M+M-1} & \dots\\
 *  \end{array}\right)
 *  @f]
 *
 *  @param x  1d signal (1d vector, n elements)
 *  @param M  Number of polyphase branches, e.g. decimation ratio.
 *  @returns  Polyphase array (M rows, (n+M-1)/M columns)
 *  @sa iforme_polyphase()
 */
template<typename T>
  Tableau<T> forme_polyphase(const Eigen::Ref<const Vector<T>> x, unsigned int M);

/** @brief Compute the standard form from the polyphase representation
 *
 *  <h3>Calcul de la forme standard à partir de la forme polyphase</h3>
 *
 *  @param X Forme polyphase (tableau 2d)
 *  @returns Forme normale (tableau 1d)
 *
 *  @sa forme_polyphase() */
template<typename T>
  Vector<T> iforme_polyphase(const Eigen::Ref<const Tableau<T>> X);

/** @brief Transformée bilinéaire : conversion transformée de Laplace vers transformée en z.
 *
 * <h3>Transformée bilinéaire</h3>
 *
 *  La transformée bilinéaire permet de convertir une fonction de transfert analogique (transformée
 *  de Laplace) en fonction de tranfert digitale (transformée en z).
 *
 *  Le calcul de la transformée en z est fait en remplaçant la variable s par :
 *  @f[
 *  s \to 2 f_e \cdot \frac{1 - z^{-1}}{1 + z^{-1}}
 *  @f]
 *
 *  @param ha Transformée de Laplace du système (fraction rationnelle)
 *  @param fe Fréquence d'échantillonnage
 *  @returns  Transformée en z
 *
 *  @sa fd_vers_fa(), fa_vers_fd()
 */
extern FRat<cfloat> trf_bilineaire(const FRat<cfloat> &ha, float fe);


extern float wa_vers_wd(float wa, float fe);
extern float wd_vers_wa(float wd, float fe);

/** @brief Conversion fréquence digitale vers fréquence analogique
 *
 * <h3>Fréquence digitale vers fréquence analogique</h3>
 *
 *  Cette fonction convertit une fréquence digitale vers fréquence analogique (pré-warping pour la transformée bilinéaire).
 *
 *  @param fd Fréquence digitale (entre 0 et 0.5)
 *  @return   Fréquence analogique (non bornée, positive)
 *
 *  La fréquence analogique est calculée ainsi :
 *  @f[
 *  f_a = \frac{\tan(\pi  f_d)}{\pi}
 *  @f]
 *
 *  Ce mapping est celui associé à la transformée bilinéaire.
 *  */
extern float fd_vers_fa(float fd);

/** @brief Fréquence analogique vers fréquence digitale
 *
 *  <h3>Fréquence analogique vers fréquence digitale</h3>
 *
 *  @param fa Fréquence analogique (non bornée, positive)
 *  @return   Fréquence digitale (entre 0 et 0.5)
 *
 * Fonction inverse de @ref fd_vers_fa() */
extern float fa_vers_fd(float fa);

/** @} */


/**  @addtogroup filtrage-tr
  *  @{ */


/** @brief Retarde le signal d'entrée de @f$n@f$ échantillons.
 *
 *  <h3>Ligne à retard</h3>
 *
 *  La fonction produit autant d'échantillons en sortie qu'il y en a en entrée.
 *  Les échantillons précédant le premier échantillon d'entrée sont supposés nuls.
 *
 *  @param n Délais entier (@f$n \geq 0@f$)
 *
 *
 *  @par Exemple
 *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_ligne_a_retard
 *  @image html filtrage-ligne-a-retard.png width=800px
 *  */
template<typename T>
  sptr<FiltreGen<T>> delay_line(unsigned int n)
  {
    return tsdf::ligne_a_retard<T>(n);
  }

struct HilbertTransformeurConfig
{
  unsigned int ntaps = 63;
  std::string fenetre = "hn";
};

/** @brief Définit un transformeur de Hilbert, qui convertit un signal réel en un signal analytique (complexe).
 *
 * <h3>Transformeur de Hilbert</h3>
 *
 * Ce bloc calcul un signal analytique (complexe) à partir d'un signal réel,
 * par recomposition du signal initial retardé avec le signal filtré avec le filtre de Hilbert :
 *
 * @image  html fig-hilbert-transfo.PNG width=200px
 *
 * @param  n        Ordre du filtre
 * @param  fenetre  Choix de la fenêtre (voir @ref fenetre())
 * @return          Filtre float -> cfloat
 *
 * @sa design_rif_hilbert(), hilbert(), hilbert_tfd()
 */
extern sptr<Filtre<float, cfloat, HilbertTransformeurConfig>>
  hilbert_transformeur(unsigned int n = 31, const std::string &fenetre = "hn");




/** @cond private */
extern void fsfir_plot(const ArrayXcf &f, unsigned int ntaps = 0);
/** @endcond */

/** @brief Implémentation directe d'une filtre RIF (Réponse Impulsionnelle Finie)
 *
 * <h3>Implémentation directe d'une filtre RIF (Réponse Impulsionnelle Finie)</h3>
 *
 *  Implémenté suivant l'équation :
 *  @f[
      y_n = \sum h_k x_{n-k}
    @f]

 *  @param h Vector des coefficients du filtre
 *  @tparam T Type des données à traiter (float, cfloat, etc.)
 *  @tparam Tc Type des coefficients
 *  @return Bloc T -> T
 *
 *  @sa filtre_rif_fft()
 *  @note Si le filtre a un nombre important de coefficients, utilisez plutôt @ref filtre_rif_fft(), qui sera plus efficace
 *  (filtrage dans le domaine fréquentiel).
 */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rif(const Eigen::Ref<const Vector<Tc>> h);


/** @brief Implémentation efficace d'un filtre RIF par FFT
 *
 *  <h3>Implémentation fréquentielle d'un filtre RIF</h3>
 *
 *  Ce bloc réalise un filtrage RIF via la technique OLA (Ovelap-And-Add).
 *  La complexité est donc bien moindre que l'implémentation standard si le nombre de coefficients @f$M@f$
 *  est important (de l'ordre de @f$\log M@f$ opérations par échantillons au lieu de @f$M@f$).
 *
 *  @param h Vector des coefficients du filtre
 *  @tparam T Type de données à traiter
 *  @return Filtre T -> T
 *
 *  Notez que cette technique introduit un délais un peu plus important que l'implémentation temporelle.
 *
 *  @sa filtre_rif()
 */
template<typename T>
  sptr<FiltreGen<T>> filtre_rif_fft(const ArrayXf &h);



/** @brief Filtre RII (forme directe I), non recommandée (utiliser plutôt @ref filtre_sois() à la place).
 *
 *  <h3>Filtrage RII (forme directe I)</h3>
 *
 *  Ce bloc implémente un filtre à Réponse Impulsionnelle Infinie, sous la forme la plus simple (directe I),
 *  c'est-à-dire que le filtre est décomposé ainsi :
 *  @f[
 *  H(z) = \frac{b_0 + b_1 z^{-1} + \dots}{a_0 + a_1 z^{-1} + \dots}
 *       = \left(b_0 + b_1 z^{-1} + \dots\right) \cdot \frac{1}{a_0 + a_1 z^{-1} + \dots}
 *  @f]
 *  (le filtre RIF correspondant au numérateur est calculé en premier, ensuite le filtre purement récursif est calculé).
 *
 *  @param h    Fonction de transfert
 *  @tparam T   Type de données à traiter
 *  @tparam Tc  Type des coefficients
 *  @return     Filtre T -> T
 *
 *  @sa @ref filtre_sois()
 *
 *  @attention Si l'ordre du filtre est important, du fait des erreurs de troncature, cette implémentation a de fortes chances
 *  de diverger. L'implémentation sous forme de cascade de filtres RII du second ordre (@ref filtre_sois()) est alors recommendée.
 *
 * */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rii(const FRat<Tc> &h);


/** @brief Création d'un filtre CIC, opérant sur des vecteurs de type T,
 *  et travaillant en interne avec le type Ti.
 *
 *  <h3>Création d'un filtre CIC</h3>
 *
 *  @param R Decimation ratio
 *  @param N Number of integrators / differentiators stages
 *  @param M Optionnal CIC design parameter (default is M=1)
 *  @param mode 'd' for decimation or 'u' for upsampling.
 *  @tparam T Type d'entrée / sortie
 *  @tparam Ti Type pour les calculs interne (int, int64_t, ...)
 *
 *  @note Pour le type interne, il faut absolument choisir un type entier, car la façon dont le filtre est
 *  implémenté fait que les calculs ne fonctionneront pas avec un type flottant.
 *
 *
 *  @par Exemple pour l'interpolation
 *  @snippet exemples-tsd-pub/src/ex-cic.cc exemple_cic_upsampling
 *  @image html filtrage-cic-interpolation.png width=800px
 *  Notez les repliements du signal utile sur le spectre.
 *  <br>
 *  @par Exemple pour la décimation
 *  @snippet exemples-tsd-pub/src/ex-cic.cc exemple_cic_decimation
 *  @image html filtrage-cic-decimation.png width=800px
 */
template<typename T, typename Ti>
  sptr<FiltreGen<T>> filtre_cic(int R, int N, int M = 1, char mode = 'd');

/** @brief This function computes and shows the frequency response of a CIC filter and
 *  then analyse the aliasing that occurs after decimation.
 *
 * @param R Decimation ratio
 * @param N Number of integrators / differentiators
 * @param M Design parameter (typically M=1)
 * @param Fin Input sample frequency
 * @param Fint Fréquence où mesurer l'atténuation
 * @return A tuple with :
 *  - Theorical transfert function (not taking into account the decimation),
 *  - Number of additionnal bits needed to implement the CIC filter
 *
 * This functions draws two plots:
 *   - The first plot shows the frequency response before decimation (both the global one, and one centered on the passband).
 *   - The second plot shows the effect of the decimation (aliasing in the baseband).
 *
 * Also, this function computes the number of additionnal bits needed to implement the filter in
 * fixed point / integer arithmetic. This is computed as
 * (see http://www.tsdconseil.fr/log/scriptscilab/cic/cic-en.pdf):
 * @f[n_{bits} = N\cdot\log_2(R) - 1@f]
 * @par Example
 * @code
 *   // 10 MHz input sample frequency
 *   Fin = 10e6;
 *   // Decimation ratio = 16, 5 CIC stages
 *   // (sample rate at output of CIC is 10 MHz / 16 = 625 KHz)
 *   R = 16, N = 5, M = 1;
 *   cic_analysis(R, N, M, Fin);
 * @endcode
 * @todo ex_cic_analysis1.png
 * Frequency response of a CIC filter, before decimation
 *
 * @todo ex_cic_analysis0.png
 * Frequency response of a CIC filter, and aliasing, after decimation
 * @sa design_cic()
 **/
extern std::tuple<FRat<float>, int, ArrayXf, ArrayXf> cic_analyse(int R, int N, int M, float Fin, float Fint = 0);



/** @brief Frequency response of a CIC filter
 *
 * @param R decimation / interpolation factor
 * @param N number of integrators / differenciators
 * @param M design parameter
 * @param f normalized frequency, between -0.5 and 0.5 (can be a 1d vector)
 * @returns mag Output magnitude, computed at each frequency point
 *
 * This function computes the CIC frequency response (in magnitude) at
 * the specified normalized frequencies @f$f_k@f$.
 * The magnitude is computed as:
 * @f[
 *   \left|H(f)\right| = \left|\frac{\sin(RM \pi f)}{RM\sin \pi f}\right|^N
 * @f]
 *
 * @par Example
 * @code
 *  R = 4, N = 8, M = 1;
 *  f = linspace(0,0.5,512);
 *  mag = cic_freq(R,N,M,f);
 *  clf(); plot(f,20*log10(mag+1e-10));
 * @endcode
 * @todo ex_cic_freq.png
 * CIC filter frequency response (SINC^N)
 *
 *  @sa cic_transfert
 */
extern ArrayXf cic_freq(int R, int N, int M, const ArrayXf &f);




using RIIStructure = tsd::filtrage::RIIStructure;


/** @brief Création d'un filtre RII sous forme d'une chaine de sections du second ordre.
 *
 *  <h3>Filtrage RII (chaine de sections du second ordre)</h3>
 *
 *  La fonction de transfert passée en entrée est factorisée sous la forme d'une cascade
 *  de sections du second ordre
 *  (plus éventuellement une section du premier ordre si l'ordre est impair), permettant une
 *  implémentation efficace d'une filtre RII :
 *
 *  @f[
 *  H(z) = G\cdot\prod_{i=1}^{M} H_i(z) = G\cdot\prod_{i=1}^{M} \frac{1 + b_1^{(i)} z^{-1} + b_2^{(i)} z^{-2} }{1 - a_1^{(i)} z^{-1} - a_2^{(i)} z^{-2}}
 *  @f]
 *
 *  @param h Fonction de transfert, normalement complexe
 *  @param structure RIIStructure::FormeDirecte1 ou RIIStructure::FormeDirecte2 (voir notes ci-dessous)
 *  @return Filtre réel vers réel
 *
 *  La forme directe 2 est légérement plus efficace que la forme 1. Cependant, si les coefficients du filtre
 *  sont amenés à changer en cours de fonctionnement, il vaut mieux utiliser la première forme, car la forme 2
 *  risque de générer des discontinuités.
 *
 *  @attention Il est recommandé que la fonction de transfert passée en entrée soit sous la forme pôles / zéros,
 *  et non pas sous la forme de la liste des coefficients, car avec cette dernière représentation, la position
 *  des pôles et zéros (et donc la réponse fréquentielle) est très instable (très sensible aux erreurs de troncature
 *  lors de la quantification des coefficients). Notez que la fonction @ref design_riia() renvoie bien une fonction
 *  de transfert sous la forme pôles / zéros.
 *
 *  @sa filtre_rii(), design_riia()
 */
template<typename T>
  sptr<FiltreGen<T>> filter_sois(const FRat<cfloat> &h, RIIStructure structure = RIIStructure::FormeDirecte2)
  {
    return tsdf::filtre_sois<T>(h, structure);
  }


template<typename T>
  sptr<FiltreGen<T>> filter_sois(const FRat<float> &h, RIIStructure structure = RIIStructure::FormeDirecte2)
  {
    return tsdf::filtre_sois<T>(h, structure);
  }

/** @brief %Filtre RII du premier ordre (dit "RC numérique")
 *
 *  <h3>%Filtre RII du premier ordre</h3>
 *
 *  Ce filtre, dit "RC numérique", ou "filtre exponentiel", est un des filtres les plus simples, puisqu'il est
 *  entiérement caractérisé par un seul coefficient.
 *
 *  @param γ Facteur d'oubli
 *
 *  Le filtre RII du premier ordre peut être défini par l'équation :
 *  @f[
 *   y_{n+1} = (1-\gamma) * y_n + \gamma * x_n = y_n + \gamma * (x_n - y_n)
 *  @f]
 *
 *  Le paramètre @f$\gamma@f$ peut être réglé facilement en fonction du temps de réponse ou de la
 *  fréquence de coupure souhaitée (voir @ref rii1_coef()).
 *
 *  @sa rii1_fcoupure(), rii1_coef()
 *   */
template<typename T>
  sptr<FiltreGen<T>> filter_iir1(float γ)
  {
    return tsdf::filtre_rii1<T>(γ);
  }


/** @brief %Filter for DC (low frequencies) suppression.
 *
 * <h3>%Filter for DC (low frequencies) suppression</h3>
 *
 *  Implémente le filtre :
 *  @f[
 *  y_k = x_k - x_{k-1} + \alpha y_{k-1}
 *  @f]
 *  d'après The DC Blocking Filter, J M de Freitas, 2007
 *
 *  @param fc Fréquence de coupure normalisée (entre 0 et 0.5)
 *
 *  @sa design_bloqueur_dc()
 *
 * @par Exemple
 *  @snippet exemples/src/ex-filtrage.cc exemple_filtre_dc
 *  @image html bloqueur-dc-ex.png width=1000px
 *
 */
template<typename T>
  sptr<FiltreGen<T>> filter_dc(float fc)
  {
    return tsdf::filtre_dc<T>(fc);
  }


/** @brief Moving average filter.
 *
 *  <h3>Moving average filter</h3>
 *
 *
 *  Ce filtre est une implémentation optimisée d'une moyenne glissante.
 *
 *  @param K profondeur de la moyenne
 *
 *  @f[
 *  y_n = \frac{1}{N}\cdot\sum_{k=0}^{K-1} x_{n-k}
 *  @f]
 *  @remark Ce filtre est implémenté efficacement grâce à un intégrateur et à un peigne de profondeur @f$K@f$ :
 *  @f[
 *   y_n = y_{n-1} + \frac{x_n - x_{n-K}}{K}
 *  @f] */
template<typename T>
  sptr<FiltreGen<T>> filter_ma(unsigned int K)
  {
    return tsdf::filtre_mg<T>(K);
  }



/** @} */



/** @addtogroup filtrage-fini
  *  @{ */

/** @brief Filtrage d'un signal par un filtre défini par sa fonction de transfert.
 *
 *  Filtrage d'un signal par un filtre défini par sa fonction de transfert.
 *
 *  Cette  fonction ne fonctionne que sur un signal de durée finie. Pour filtrer des données
 *  reçues au fil de l'eau, il faut utiliser une des structure avec contexte (voir
 *  @ref filtre_rif(), @ref filtre_sois(), etc.).
 *
 *  @param h Fonction de transfert ou coefficients du filtre à appliquer (peut être un filtre RIF ou RII)
 *  @param x Signal d'entrée à filtrer
 *  @returns Signal filtré
 *
 *  Pour un filtre RIF, cette fonction calcule le produit de convolution :
 *  @f[
 *    y_n = (h \star x)_n = \sum_{k=0}^{K-1} h_k\cdot x_{n-k}
 *  @f]
 *
 *  Calcul autant d'échantillons de sortie qu'il y a d'échantillons d'entrée
 *  (calcul de @f$y_n@f$ suivant la formule ci-dessus pour @f$n=0\dots N-1@f$),
 *  en introduisant des zéros avant le signal.
 *
 *  */
template<typename T, typename Tc>
Vector<T> filter(const FRat<Tc> &h, const Vector<T> &x)
{
  return tsdf::filtrer(h, x);
}

template<typename T, typename Tc>
Vector<typename T::Scalar> filter(const Vector<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  return tsdf::filtrer(h, x);
}

/** @brief Filtrage zéro-phase (bi-directionnel)
 *
 *  Cette fonction permet de filtrer un signal de durée finie sans décalage temporel, à partir d'un filtre RIF quelquonque,
 *  qui est appliqué deux fois sur le signal : une fois normalement, et une fois dans le sens inverse du temps :
 *  @f[
 *    y = \left((x \star h)_{-.} \star h\right)_{-.}
 *  @f]
 *
 *  @param h Coefficients du filtre
 *  @param x Signal à filtrer
 *  @return Signal filtré, sans décalage temporel.
 *
 *  */
template<typename T, typename Tc>
  Vector<typename T::Scalar> filtfilt(const Vector<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  return tsdf::filtfilt(h, x);
}


template<typename T, typename Tc>
  Vector<T> convol(const Eigen::Ref<const Vector<Tc>> h, const Eigen::Ref<const Vector<T>> x)
{
  return tsdf::convol(h, x);
}

/** @brief Calcul du signal analytique (via un filtrage RIF).
 *
 *  <h3>Calcul du signal analytique (via un filtrage RIF)</h3>
 *
 *  @param x Signal d'entrée
 *  @param ncoefs Nombre de coefficients pour l'approximation RIF
 *  @return Signal analytique (complexe)
 *  @sa design_rif_hilbert(), hilbert_tfd(), hilbert_transformeur() */
inline ArrayXcf hilbert(IArrayXf x, unsigned int ncoefs = 127)
{
  return tsdf::hilbert(x, ncoefs);
}

/** @brief Calcul du signal analytique (via la TFD)
 *
 *  <h3>Calcul du signal analytique (via la TFD)</h3>
 *
 *  Cette fonction met tout simplement à zéro les fréquences négatives du signal,
 *  dans le domaine fréquentiel.
 *  Contrairement à @ref hilbert(), cette technique n'introduit aucun délais dans le domaine temporel.
 *
 *  @attention
 *  Du fait des hypothèses sous-jacentes à la la TFD, des artefacts peuvent apparaitre
 *  aux bords du signal.
 *
 *  @param x Signal d'entrée
 *  @return Signal analytique (complexe)
 *  @sa design_rif_hilbert(), hilbert(), hilbert_transformeur() */
inline ArrayXcf hilbert_dft(IArrayXf x)
{
  return tsdf::hilbert_tfd(x);
}

/** @} */


/** @addtogroup filtrage-rythme-itrp
  *  @{ */

/** @brief Interface générique pour un interpolateur */
template<typename T>
  using Interpolator = tsdf::Interpolateur<T>;

/** @brief %Interpolateur générique à base filtre RIF */
template<typename T>
  using InterpolatorFIR = tsdf::InterpolateurRIF<T>;


/* Interpolateur demi-bande */
template<typename T> sptr<FiltreGen<T>> filtre_allpass_ups();

/* Décimateur demi-bande */
template<typename T> sptr<FiltreGen<T>> filtre_allpass_decim();

/** @brief %Interpolation par splines cardinales */
template<typename T>
  sptr<InterpolatorFIR<T>> itrp_cspline();

/** @brief %Interpolation linéaire (équivalente à Lagrange degré 1) */
template<typename T>
  sptr<InterpolatorFIR<T>> itrp_lineaire();

/** @brief %Interpolateur de Lagrange (la fonction sinus cardinal est interpolée par un polynôme) */
template<typename T>
  sptr<InterpolatorFIR<T>> itrp_lagrange(unsigned int degre);

/** @brief %Interpolateur à sinus cardinal fenêtré */
template<typename T>
  sptr<InterpolatorFIR<T>> itrp_sinc(int ncoefs, float fcut = 0.5, const std::string &type_fenetre = "hn");


/** @} */

/** @addtogroup filtrage-rythme
  *  @{ */


/** @brief Filtrage RIF avec post-décimation (ne calcule qu'une échantillon sur @f$R@f$).
 *
 * <h3>Filtrage RIF et décimation</h3>
 *
 * Cette structure permet l'application successive d'un filtrage RIF et d'une décimation,
 * de manière efficace (sans calcul des échantillons supprimés).
 * Le flux de sortie est décimé d'un facteur @f$R@f$ par rapport au flux d'entrée.
 *
 * @param R Facteur de décimation
 * @param h Coefficients du filtre RIF anti-repliement.
 * @tparam Tc Type de coefficient
 * @tparam T  Type d'entrée / sortie */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filter_fir_decim(const Eigen::Ref<const Vector<Tc>> h, unsigned int R)
  {
    return tsdf::filtre_rif_decim(h, R);
  }

/** @brief Filtrage RIF avec pré-insertion de zéro (implémentation polyphase)
 *
 * <h3>Filtrage RIF et interpolation</h3>
 *
 * Cette structure permet le calcul efficace de la succession d'un sur-échantillonnage (insertion de R-1 zéros entre chaque échantillon d'entrée),
 * et d'un filtrage (typiquement anti-repliement) RIF, grâce à une décomposition polyphase
 * du filtre RIF.
 *
 *
 *  @param h  Coefficients du filtre RIF anti-repliement.
 *  @param R  Facteur d'upsampling.
 *  @tparam Tc Type des coefficients (float, double, etc.)
 *  @tparam T  Type des échantillons d'entrée / sortie (float, double, cfloat, etc.)
 *  @return un filtre (type T vers T), le flux de sortie étant upsamplé d'un facteur @f$R@f$ par rapport au flux d'entrée. */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filter_fir_ups(const Eigen::Ref<const Vector<Tc>> h, unsigned int R)
  {
    return tsdf::filtre_rif_ups(h, R);
  }

/** @brief Rythm adapter, with arbitrary ratio.
 *
 *  <h3>Rythm adapter, with arbitrary ratio</h3>
 *
 *  Cette fonction renvoi un bloc de ré-échantillonnage, permettant de changer le rythme d'un signal reçu au fil de l'eau.
 *  L'implémentation est basée sur une cascade de décimateurs (si le ratio est inférieur à 1) ou d'interpolateurs (si le ratio est supérieur à 1) demi-bandes,
 *  suivis d'un interpolateur de ratio arbitraire.
 *
 *  @param ratio Ratio de décimation / interpolation (rapport entre les fréquences d'échantillonnage de sortie et d'entrée).
 *  @sa resample(), filtre_rif_ups(), filtre_rif_decim()
 **/
template<typename T>
  sptr<FiltreGen<T>> filter_resample(float ratio)
  {
    return tsdf::filtre_reechan<T>(ratio);
  }

/** @brief Interpolation d'un ratio arbitraire (calcul au fil de l'eau)
 *
 * <h3>Interpolation (ratio arbitraire)</h3>
 *
 *  @param ratio Ratio de décimation / interpolation (rapport entre les fréquences d'échantillonnage de sortie et d'entrée).
 *  @param itrp Interpolateur générique
 *  @sa itrp_cpline(), itrp_lineaire(), itrp_lagrange(), itrp_sinc() */
template<typename T> sptr<FiltreGen<T>> filter_itrp(float ratio, sptr<Interpolator<T>> itrp = itrp_cspline<T>())
    {
    return tsdf::filtre_itrp<T>(ratio, itrp);
    }


/** @} */




}






