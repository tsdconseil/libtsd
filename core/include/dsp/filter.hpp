#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

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


  /**  @addtogroup filtrage
    *  @{ */

  /** @brief Vérification de la validité d'une fréquence normalisée.
   *
   * <h3>Vérification de la validité d'une fréquence normalisée</h3>
   *
   * Vérifie que la fréquence est comprise entre 0 et @f$0{,}5@f$.
   *
   * En cas d'échec, léve une exception.
   *
   * @param f Fréquence normalisée.
   * @param msg Descriptif optionnel (affiché en cas d'échec).
   *
   * @sa echec(), msg_erreur()
   */
  inline void verifie_frequence_normalisee(float f, const std::string &msg = "")
  {
    return tsdf::verifie_frequence_normalisee(f, msg);
  }


  /** @} */


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



  /** @brief Slepian window creation.
   *
   *
   */
  static inline ArrayXf window_slepian(int N, float B)
  {
    return tsdf::fenêtre_slepian(N, B);
  }


  using FenInfos = tsdf::FenInfos;

  static inline FenInfos window_analysis(const std::string &nom, const ArrayXf &x)
  {
    return tsdf::fenetre_analyse(nom, x);
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
 *  @return A tuple of 2 vectors: the frequencies @f$f_k@f$ (normalized, between 0 and 0.5), and the magnitudes @f$y_k@f$.
 *
 *  @par Example
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


/** @brief Impulse response. */
template<typename T>
  ArrayXf repimp(const FRat<T> &h, int npts = -1)
{
  return tsdf::repimp(h, npts);
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

/** @brief Computes the group delay.
 *
 *  <h3>Group delay</h3>
 *
 *  @f[
 *  G(\omega) = \frac{d\arg H(\omega)}{d\omega}
 *  @f]
 *
 *  @param h Transfert function to be analyzed.
 *  @param npts Frequency resolution.
 *  @return A tuple of 2 vectors : the frequencies @f$f_k@f$ (normalized, between 0 and 0.5),
 *  et and the group delay (expressed in number of samples).
 *
 *  @par Example
 *  @snippet exemples/src/ex-filtrage.cc ex_frgroup
 *  @image html frgroup.png width=600px
 *
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frgroup(const FRat<T> &h, unsigned int npts = 1024)
{
  return tsdf::frgroup<T>(h, npts);
}




/** @brief Analysis of a linear filter (plot the different responses).
 *
 *  <h3>Analysis of a linear filter</h3>
 *
 * This function build a new figure, and plots the frequency
 * and time response of the filter,
 * alongside with the zeros / poles diagram.
 *
 *  @param h      Transfert function to be analyzed.
 *  @param fs     Sampling frequency.
 *  @return       The new figure.
 *
 * @par Example
 * @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_analyse
 * @image html filtrage-analyse.png width=1000px
 *
 * @sa filter_display()
 */
template<typename T>
  sptr<dsp::view::Figures> filter_analysis(const FRat<T> &h, float fs = 1.0f)
{
  return tsdf::analyse_filtre(h, fs);
}

/** @brief Display impulse and frequency response of a filter.
 *
 *  <h3>Display impulse and frequency response of a filter.</h3>
 *
 * This function build a new figure, and plots the frequency
 * and time response of the filter.
 *
 *  @param h      Transfert function to be analyzed.
 *  @param fs     Sampling frequency.
 *  @return       The new figure.
 *
 * @par Exeaple
 * @snippet exemples/src/filtrage/ex-filtrage.cc exemple_affiche
 * @image html filtrage-affiche.png width=1000px
 *
 * @sa filter_analysis() (for more plots)
 */
template<typename T>
  tsd::vue::Figures filter_display(const FRat<T> &h, float fs = 1.0f)
{
  return tsdf::affiche_filtre(h, fe);
}


/** @cond undoc */


template<typename T> ArrayXf repimp(const Vector<T> &h, int npts = -1)
{
  return tsdf::repimp(h, npts);
}
inline ArrayXcf repfreq(const ArrayXf &h, const ArrayXf &fr)
{
  return tsdf::repfreq(h, fr);
}
template<typename T> ArrayXf repfreq(const Vector<T> &h, int npts = 1024)
{
  return tsdf::repfreq<T>(h, npts);
}

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

inline dsp::view::Figures filter_analysis(const ArrayXf &h, float fe = 1.0f)
{
  return tsdf::analyse_filtre(h, fe);
}

inline dsp::view::Figures filter_display(const ArrayXf &h, float fe = 1.0f)
{
  return tsdf::affiche_filtre(h, fe);
}

/** @endcond */





/** @brief Poles / zeros diagram.
 *
 *  <h3>Poles / zeros diagram</h3>
 *
 *  The transfert function is factored as:
 *  @f[
 *  H(z) = \frac{\prod z - z_i}{\prod z - p_i}
 *  @f]
 *
 *  where the @f$z_i@f$ and @f$p_i@f$ are called respectively
 *  the zeros et poles.
 *
 *  @param h transfert function
 *  @param fig Figure where to plot the zeros and poles.
 *
 *
 *  @par Example
 *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_plz
 *  @image html filtrage-plz.png width=600px
 *  */
template<typename T>
  void plot_plz(dsp::view::Figure &fig, const FRat<T> &h)
{
  tsdf::plot_plz(fig, h);
}



/** @brief Amplitude response of a symetrical or antisymetrical FIR filter (linear phase).
 *
 *  <h3>Amplitude response of a symetrical or antisymetrical FIR filter</h3>
 *
 *  This function computes the amplitude response @f$A(\omega)@f$, for a real FIR filter
 *  with linear phase (the coefficients must be symetrical or antisymetrical around the central point).
 *
 *  For a symetrical filter, the response is:
 *  @f[
 *  A(\omega) = h_M + 2 \sum_{n=0}^{M-1} h_n \cos((M-n)\omega)
 *  @f]
 *
 *  @param L Frequency resolution.
 *  @param h Filter coefficients.
 *  @param symetrical True if the coefficients are symetrical around @f$(N-1)/2@f$ (FIR filter of type I or II).  Otherwise, a filter of type III or IV is assumed.
 *  @returns A tuple of two vectors: frequency vector (normalized, between 0 and 0.5),
 *  and the amplitude response vector.
 *
 *  @warning The real symetry (or antisymetry) of the coefficients is not checked checked!
 *
 *  @sa frmag(), frphase()
 */
inline std::tuple<ArrayXf, ArrayXf> firamp(const Eigen::ArrayXf &h, int L = 1024, bool symetrique = true)
{
  return tsdf::rifamp(h, L, symetrique);
}



/** @brief Computes the delay of a linear phase FIR filter.
 *
 *  <h3>Computes the delay of a linear phase FIR filter</h3>
 *
 *  This function returns the delay, expressed in number of samples,
 *  due to a linear phase FIR filter
 *  (e.g. with symetrical or anti-symetrical coefficients):
 *  @f[
 *  \tau = \frac{N-1}{2}
 *  @f]
 *
 *  @param N  Number of coefficients
 *  @returns Delay @f$\tau@f$ of the filter.
 *
 */
inline float fir_delay(int N)
{
  return tsdf::rif_delais(N);
}


/**  @}
  *  @addtogroup filtrage-design
  *  @{ */



struct SpecFreqIntervalle
{
  float fb, fh;
  float atten = 1.0f;
  float poids = 1.0f;
};

/** @brief Hilbert filter FIR approximation.
 *
 *  <h3>Hilbert filter FIR approximation</h3>
 *
 *  This function computes the windowed theoretical temporal response of a Hilbert filter:
 *  @f[
 *      h_k = \frac{2}{k\pi} \cdot \sin(k \pi / 2)^2 \cdot w_k;
 *  @f]
 *
 *  @param n Filter order
 *  @param fenetre Window type (by default, Hann window)
 *  @returns FIR filter coefficients
 *
 *  @par Example
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_hilbert
 *  @image html design-rif-hilbert.png width=600px
 *
 *  @sa hilbert(), hilbert_transformeur()
 */
inline ArrayXf design_fir_hilbert(int n, const std::string &window = "hn")
{
  return tsdf::design_rif_hilbert(n, window);
}

using BiquadSpec = tsdf::BiquadSpec;



/** @brief Biquad filter design.
 *
 * <h3>Biquad filter design</h3>
 *
 * For a more complete description, see @ref design_biquad().
 *
 * @param spec Specification (type, cutoff frequency, quality factor, etc.)
 *
 * @sa design_biquad()
 */
inline FRat<float> design_biquad(const BiquadSpec &spec)
{
  return tsdf::design_biquad(spec);
}

/** @brief Biquad filter design.
 *
 * <h3>Biquad filter design</h3>
 *
 * These second order RII filter are adapted from analog prototype through the bilinear transform.
 *
 * The analog  prototypes are the following ones (defined for a cutoff pulsation of 1 radian/s, and @f$Q@f$ being the quality factor):
 *
 * - Low-pass filter:
 * @f[
 * H(s) = \frac{1}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - High-pass filter:
 * @f[
 * H(s) = \frac{s^2}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - Band-pass filter:
 * @f[
 * H(s) = \frac{s/Q}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - Band-stop filter:
 * @f[
 * H(s) = \frac{s^2+1}{s^2+\frac{1}{Q}s+1}
 * @f]
 *
 *
 * @param type    Filter type ("lp", "hp", "bp", "sb", ...).
 * @param f       Normalized cut-off frequncy (or central frequency for band-pass and band-stop), between 0 and 0.5.
 * @param Q       Quality factor (note: for @f$Q>1/\sqrt(2)\sim 0{,71}@f$, there will be some resonnance).
 * @param gain_dB Gain, in dB, used only for the filter of type resonnance or plateau.
 *
 * @par Example: low-pass filters, with different values for the quality factor
 * @snippet exemples/src/filtrage/ex-filtrage.cc ex_biquad_lp
 * @image html ex-biquad-pb.png width=800px
 *
 * @par Bibliography
 * - <i>Cookbook formulae for audio equalizer biquad filter coefficients,</i> Robert Bristow-Johnson,
 *      https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html,
 * - <i>F0 and Q in filters, Mini tutorial,</i> Analog Devices,
 */
inline FRat<float> design_biquad(const std::string type, float f, float Q, float gain_dB = 0)
{
  return tsdf::design_biquad(type, f, Q, gain_dB);
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

extern FRat<cfloat> design_riia(int n, TypeFiltre type,
    PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc);

extern FRat<cfloat> design_riia_laplace(int n, TypeFiltre type, PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc);

/** @endcond */


/** @brief IIR design from a classical analog prototype.
 *
 * <h3>IIR design from a classical analog prototype</h3>
 *
 * This function returns a <b>discrete</b> transfert function,
 * in the form of poles and zeros
 * (ideal for an implementation with second order sections, see @ref filter_sois()).
 * The following prototypes are possible:
 *  - <b>Butterworth</b> (no ondulation, wide transition band)
 *  - <b>Chebychev type I</b> (ondulations in the pass-band)
 *  - <b>Chebychev type II</b> (ondulations in the stop-band)
 *  - <b>Elliptique</b> (ondulations everywhere, but shortest transition band)
 *
 *  The filter is first designed in the analog domain (Laplace transform),
 *  then converted in digital form (z transform) with the
 *  bilineair transform.
 *
 * @param n           Fitler order.
 * @param type        Filter type ("lp" for low-pass, "hp" for high-pass, ...)
 * @param prototype   "butt", "cheb1", "cheb2" or "ellip"
 * @param fc          Normalized cut-off frequncy (or central frequency for band-pass and band-stop), between 0 and 0.5.
 * @param δ_bp        Maximal ondulation (in decibels) in the pass-band (used only for Chebychev type I or elliptic filters).
 * @param δ_bc        Minimal attenuation (in decibels) in the sto-band (used only for Chebychev type II or elliptic filters).
 * @return h          Transfert function (digital)
 *
 * @par Example
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




/** @brief Frequency sampling design.
 *
 *  <h3>Frequency sampling design</h3>
 *
 * This method enables to approximate with a FIR filter with @f$n@f$ coefficients (@f$n@f$ being odd)
 * an arbitrary frequency response, given as a specification.
 * The frequency specification must be given as a vector @f$d@f$ of @f$m=\frac{n+1}{2}@f$ real elements, of type:
 * @f[
 * d_k = H(f_k),\ f_k = k \cdot \frac{1}{2m-1},\ k = 0,\dots, m-1
 * @f]
 *
 *
 * The frequencies @f$f_k@f$ can be computed with the function @ref design_fir_freq_freqs().
 *
 * @note
 * Note that the elements @f$m@f$ of the frequency response enable to specify uniquely  a real filter
 * with @f$n=2m-1@f$ coefficients (indeed, each element of the response, except for the first,
 * is used twice, for positive and negative frequencies).
 * Thus, if @f$n\neq 2m-1@f$, the desired frequency response is, before computing the filter coefficients, resampled (by linear interpolation)
 * with @f$m'@f$ values such as
 * @f$n=2m'-1@f$ (in every case, the number of coefficients must be odd).
 *
 *
 *  @param n     Filter order (must be odd).
 *  @param d     Vector with desired frequency response (positives frequencies).
 *  @returns     Vector with FIR filter coefficients (dimension = @f$n@f$).
 *
 * @par Example
 * @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_freq
 * @image html design-rif-freq.png width=800px
 *
 * @sa design_fir_freq_freqs()
 */
inline ArrayXf design_fir_freq(int n, const ArrayXf &d)
{
  return tsdf::design_rif_freq(n, d);
}


/** @brief Computes the @f$m@f$ frequencie values used for the sampling frequency technique design.
 *
 * <h3>Computes the @f$m@f$ frequencie values used for the sampling frequency technique design</h3>
 *
 * This function returns the following frequency vector:
 * @f[
 *  f_k = k \cdot \frac{1}{2m-1},\ k = 0,\dots, m-1
 * @f]
 * where @f$m=(n+1)/2@f$ is the number of frequency points to specify.
 *
 * @param n Filter order (must be odd).
 *
 * @sa design_fir_freq()
 */
inline ArrayXf design_fir_freq_freqs(int n)
{
  return tsdf::design_rif_freq_freqs(n);
}

/** @brief Equiripple / Remez FIR design.
 *
 *  <h3>Equiripple / Remez FIR design</h3>
 *
 *  @param n      Filter order
 *  @param d      Desired frequency response
 *  @param w      Weighting coefficients vector (must have the same length as d)
 *  @param debug  If true, display some plots with intermediate computings.
 *  @returns      Coefficients vector of the FIR filter (vector of length n)
 */
inline ArrayXf design_fir_eq(unsigned int n, IArrayXf d, IArrayXf w, bool debug = false)
{
  return tsdf::design_rif_eq(n, d, w, debug);
}

/** @brief Normalized cardinal sine with configurable cut-off frequency
 *
 *  <h3>Normalized cardinal sine with configurable cut-off frequency</h3>
 *
 *  @f[
 *  y(t) = \frac{\sin(2 \pi t f_c)}{\pi  t}
 *  @f]
 *
 *  The Fourier transform of this fonction is a door of width @f$\pm f_c@f$,
 *  which is an ideal (yet not realizable) prototype for a low-pass filter.
 *
 *  @param t   Temporal sampling point (as fractionnal number of samples).
 *  @param fc  Normalized cut-off frequency, between 0 and 0.5.
 *  @returns   @f$y@f$ value
 *
 *  @sa sinc(), design_fir_fen() */
inline float sinc2(float t, float fc)
{
  return tsdf::sinc2(t, fc);
}

/** @brief Normalized cardinal sine with configurable cut-off frequency at 0.5.
 *
 *  <h3>Normalized cardinal sine</h3>
 *  Compute a cardinal sine with configurable cut-off frequency at 0.5:
 *  @f[
 *  y(t) = \frac{\sin(\pi t)}{\pi  t}
 *  @f]
 *
 *  @param t Temporal sampling point (as fractionnal number of samples).
 *  @returns Normalized cut-off frequency, between 0 and 0.5.
 *
 *  @sa sinc2()
 *
 *  */
inline float sinc(float t)
{
  return tsdf::sinc(t);
}

/** @brief Windowed cardinal sine FIR design.
 *
 * <h3>Windowed cardinal sine FIR design</h3>
 *
 * @param n       Filter order.
 * @param type    Filter type ("lp" for low-pass, "hp" for high-pass, ...)
 * @param fc      Normalized cut-off frequency.
 * @param fen     Window type ("hn", "hm", "tr" ou "re"), see function @ref window()
 * @param fc2     Second cut-off frequency (only band-pass or band-stop type filters)
 * @returns       Vector with filter coefficients.
 * @sa design_fir_eq(), design_fir_freq()
 */
inline ArrayXf design_fir_wnd(unsigned int n, const std::string &type, float fc, const std::string &fen = "hn", float fc2 = 0)
{
  return tsdf::design_rif_fen(n, type, fc, fen, fc2);
}

/** @brief Windowed cardinal sine FIR design (Kaiser window).
 *
 *  <h3>Windowed cardinal sine FIR design (Kaiser window)</h3>
 *
 *  The use of a Kaiser window enable to choose both the filter  attenuation
 *  and the transition band width, the adjusting variable being the number of coefficients.
 *
 * @param type    Filter type ("lp" for low-pass, "hp" for high-pass, ...).
 * @param fc      Normalized cut-off frequency.
 * @param atten_db Desired attenuation in the stop-band (decibels).
 * @param df      Transition band width (normalized frequency units).
 * @param fc2     Second cut-off frequency (only band-pass or band-stop type filters).
 *
 * @sa design_fir_fen(), design_fir_wnd_chebychev()
 */
inline ArrayXf design_rif_wnd_kaiser(const std::string &type, float fc, float atten_db,
    float df, float fc2 = 0)
{
  return tsdf::design_rif_fen_kaiser(type, fc, atten_db, df, fc2);
}

/** @brief Windowed cardinal sine FIR design (Chebychev window).
 *
 *  <h3>Windowed cardinal sine FIR design (Chebychev window)</h3>
 *
 *  The use of a Chebychev window enable to choose both the filter attenuation
 *  and the number of coefficients, the adjusting variable being the transition band width.
 *
 * @param n       Filter order.
 * @param type    Filter type ("lp" for low-pass, "hp" for high-pass, ...).
 * @param fc      Normalized cut-off frequency.
 * @param atten_db Desired attenuation in the stop-band (decibels).
 * @param fc2     Second cut-off frequency (only band-pass or band-stop type filters).
 *
 * @sa design_fir_wnd(), design_fir_wnd_kaiser()
 */
inline ArrayXf design_fir_wnd_chebychev(int n, const std::string &type,
    float fc, float atten_db, float fc2 = 0)
{
  return tsdf::design_rif_fen_chebychev(n, type, fc, atten_db, fc2);
}

/** @brief Raised cosine filter design.
 *
 *  <h3>Raised cosine filter design</h3>
 *
 *  @param n  Filter order.
 *  @param β  Roll-off factor.
 *  @param fc Normalized cut-off frequency.
 *
 *  @par Example
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_cs
 *  @image html design-rif-cs.png width=800px
 *
 *  @sa design_fir_srrc() */
inline ArrayXf design_fir_rc(int n, float β, float fc)
{
  return tsdf::design_rif_cs(n, β, fc);
}

/** @brief Square-root raised cosine filter design.
 *
 *  <h3>Square-root raised cosine filter design</h3>
 *
 *  @param n  Filter order.
 *  @param β  Roll-off factor.
 *  @param fc Normalized cut-off frequency.
 *
 *  @par Example
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rif_rcs
 *  @image html design-rif-rcs.png width=800px
 *
 *  @sa design_fir_rc(), design_fir_srrc1()
 *
 */
inline ArrayXf design_fir_srrc(int n, float β, float fc)
{
  return tsdf::design_rif_rcs(n, β, fc);
}

/** @brief Square-root raised cosine filter design (1).
 *
 *  <h3>Square-root raised cosine filter design (1)</h3>
 *
 *  This function is equivalent to @ref design_fir_srrc(), only instead
 *  of the cut-off frequency @f$f_c@f$, it is the over-sampling factor which is specified
 *  (@f$\textrm{OSF} = \frac{1}{2 f_c}@f$).
 *
 *  @param n  Filter order.
 *  @param β  Roll-off factor.
 *  @param osf Oversampling factor (@f$\frac{f_e}{f_{symb}}@f$).
 *  @param nrm Type of normalization for the coefficients ('s' for sum = 1, 'e' for sum of squares = 1).
 *
 *  @sa design_fir_srrc()
 */
inline ArrayXf design_fir_srrc1(int n, float β, float osf, char nrm = 's')
{
  return tsdf::design_rif_rcs1(n, β, osf, nrm);
}

/** @brief Coefficients for FIR approximation of a Gaussian filter (non windowed).
 *
 * <h3>Coefficients for FIR approximation of a Gaussian filter</h3>
 *
 *  @param n   Number of taps.
 *  @param σ   Standard deviation (in number of samples).
 *
 *  @par Example
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_gaussien
 *  @image html design-rif-gaussien.png width=800px
 *
 *  @sa design_fir_gaussian_telecom()
 */
inline ArrayXf design_fir_gaussian(int n, float σ)
{
  return tsdf::design_rif_gaussien(n, σ);
}


/** @brief Coefficients for FIR approximation of a Gaussian filter, for a GFSK modulation.
 *
 *  <h3>Gaussian filter convolved with NRZ filter</h3>
 *
 *  Computes the convolution of Gaussian filter and a moving average filter, with depth equal to
 *  the over-sampling factor.
 *  The standard deviation of the Gaussian filter is computed according to the Bandwidth-Time product (BT parameter).
 *
 *  This function is targeted for the design of a shaping filter for GFSK modulation.
 *
 *  @param n   Number of taps.
 *  @param BT  Bandwidth - time product.
 *  @param osf Oversampling factor (ratio of sampling frequency versus symbol frequency).
 *  @return    Vector with filter coefficients.
 *
 *  @par Example
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_gaussien_telecom
 *  @image html design-rif-gaussien-telecom.png width=800px
 *
 *  @sa design_fir_gaussian()
 */
inline ArrayXf design_fir_gaussian_telecom(int n, float BT, int osf)
{
  return tsdf::design_rif_gaussien_telecom(n, BT, osf);
}

/** @brief Computes the standard deviation (relative to the sample period) as a function of the Bandwidt - Time product (BT). */
inline float design_fir_gaussian_telecom_BT2sigma(float BT)
{
  return tsdf::design_rif_gaussien_telecom_BT_vers_sigma(BT);
}


/** @brief Computes the FIR filter equivalent to the serie concatenation of two other FIR filters.
 *
 * <h3>Serie concatenation of two FIR filters</h3>
 *
 * This function computes the coefficients of a FIR filter that has the same response
 * as the serie cascading of two FIR filters given as parameters:
 * @f[
 * h = h_1 \star h_2
 * @f]
 *
 * The number of output coefficients is computed as: @f$n = n_1 + n_2 - 1@f$.
 *
 * @param h1 Coefficients of the first filter (@f$n_1@f$ taps)
 * @param h2 Coefficients of the second filter (@f$n_2@f$ taps)
 * @return   Coefficients of the serie concatenation (@f$n_1+n_2-1@f$ taps)
 *
 */
inline ArrayXf design_fir_prod(const ArrayXf &h1, const ArrayXf &h2)
{
  return tsdf::design_rif_prod(h1, h2);
}


using tsdf::CICConfig;





/** @brief CIC filter theorical transfert function.

    <h3>CIC filter theorical transfert function</h3>

    @param config CIC filter specifications (see @ref CICConfig)
    @returns CIC transfert function (not taking into account the decimation)
    This function computes the theorical transfert function of a CIC filter, when one does not look at the decimation effect. The CIC filter responses is defined as:
           @f[H(z) = \frac{1}{R^N}\left(1+z^{-1}+\cdots+z^{-(R-1)}\right)^N@f]
     (e.g. a cascade of @f$N@f$ moving average filters, each of identical length @f$R@f$.
 */
inline FRat<float> design_cic(const CICConfig &config)
{
  return tsdf::design_cic(config);
}

using tsdf::CICComp;

/** @brief Design of a compensation FIR filter for a CIC filter.
 *
 * <h3>Design of a compensation FIR filter for a CIC filter</h3>
 *
 *
 * @param config Paramètres principaux (voir @ref CICConfig)
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
inline CICComp design_cic_comp(const CICConfig &config, float Fin, int R2, float fc, int ncoefs)
{
  return design_cic_comp(config, Fin, R2, fc, ncoefs);
}



/** @brief Transfert function for a DC blocker.
 *
 *  <h3>Transfert function for a DC blocker</h3>
 *
 *  Computes the following transfert function:
 *  @f[
 *  H(z) = \frac{1 - z^{-1}}{1 - \alpha \cdot z^{-1}}
 *  @f]
 *  With
 *  @f[
 *  \alpha = \frac{\sqrt{3} - 2 \sin(\pi f_c)}{\sin(\pi f_c) + \sqrt{3} \cos(\pi f_c)}
 *  @f]
 *
 *  @param fc Normalized cut-off frequency (0 - 0.5)
 *
 *  @par Example
 *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_design_bloqueur_dc
 *  @image html bloqueur-dc-resp.png width=1000px
 *
 *  @par Bibliography
 *  <i>The DC Blocking Filter,</i> J.M. de Freitas, 2007
 *
 *  @sa filter_dc()
 */
inline FRat<float> design_dc_blocker(float fc)
{
  return tsdf::design_bloqueur_dc(fc);
}

/** @brief Transfert function of an exponential filter.
 *
 *  <h3>Transfert function of an exponential filter</h3>
 *  This function returns the transfert function for the following system:
 *  @f[
 *  y_n = \gamma x_n + (1-\gamma) y_{n-1}
 *  @f]
 *  That is:
 *  @f[
 *  H(z) = \frac{\gamma z}{z-(1-\gamma)}
 *  @f]
 *  where @f$\gamma@f$ is defined as a function of the desired cut-off frequency @f$f_c@f$:
 *  @f[
 *    \gamma = 1 - e^{-2\pi f_c}
 *  @f]
 *
 *  @param fc -3 dB normalized cut-off frequency.
 *  @return Digital transfert function.
 *
 *  @par Example
 *  @snippet exemples/src/ex-filtrage.cc ex_design_rii1
 *  @image html design-rii1.png width=1000px
 *
 *  @sa iir1_coef()
 */
inline FRat<float> design_iir1(float fc){return tsdf::design_rii1(fc);}

/** @brief Generate first order IIR filter forget factor (unique coefficient) from cut-off frequency
 *
 *  <h3>First order RII design</h3>
 *
 *  This function computes the forget factor for a first order IIR filter
 *  (that is, an exponential filter), as a function of the desired cut-off frequency.
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
 *  @sa design_iir1(), iir1_tc2coef()
 */
inline float iir1_coef(float fc){return tsdf::rii1_coef(fc);}

/** @brief Same as previous function, but take as input the time constant (in samples).
 *
 *  <h3>First order RII design</h3>
 *
 *  This function computes the forget factor for a first order IIR filter
 *  (that is, an exponential filter), as a function of the desired time constant.
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
 *  @param tc Time constant, in number of samples.
 *  @returns Coeficient @f$\gamma@f$ for the IIR filter
 *
 */
inline float iir1_tc2coef(float tc)
{
  return tsdf::rii1_tc_vers_coef(tc);
}

inline float iir1_coef2tc(float γ)
{
  return tsdf::rii1_coef_vers_tc(γ);
}

/** @brief Compute cut-off frequency, from forget factor of first order IIR filter.
 *
 *  @param γ Forget factor.
 *  @returns Normalized cut-off frequency.
 *
 *  @sa iir1_coef() */
inline float iir1_fcut(float γ)
{
  return tsdf::rii1_fcoupure(γ);
}


/** @brief Creation of the polyphase representation of a signal.
 *
 *  <h3>Polyphase representation of a signal</h3>
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
 *  @sa polyphase_iforme()
 */
template<typename T>
  Tableau<T> polyphase_form(const Eigen::Ref<const Vector<T>> x, unsigned int M)
{
  return tsdf::forme_polyphase(x, M);
}

/** @brief Compute the standard form from the polyphase representation.
 *
 *  <h3>Standard form from the polyphase representation</h3>
 *
 *  @param X Polyphase form (2d array)
 *  @returns Standard dorm (1d array)
 *
 *  @sa polyphase_form()
 */
template<typename T>
  Vector<T> polyphase_iforme(const Eigen::Ref<const Tableau<T>> X)
{
  return tsdf::iforme_polyphase(X);
}

/** @brief Bilinear transform: conversion from Laplace transform to z transform.
 *
 * <h3>Bilinear transform</h3>
 *
 *  The bilinear transform enables to approximate an analog transfert function (Laplace transform)
 *  with a digital tranfert function (z transform).
 *
 *  To do this, the @f$s@f$ variable in Laplace transform is replaced by the following approximation:
 *  @f[
 *  s \mapsto 2 f_e \cdot \frac{1 - z^{-1}}{1 + z^{-1}}
 *  @f]
 *
 *  @param ha Laplace transform of the system (rational fraction)
 *  @param fe Sampling frequency
 *  @returns  Z transform (rational fraction)
 *
 *  @sa fd_vers_fa(), fa_vers_fd()
 */
inline FRat<cfloat> bilinear_transform(const FRat<cfloat> &ha, float fe)
{
  return tsdf::trf_bilineaire(ha, fe);
}


inline float ωa2ωd(float ωa, float fe)
{
  return tsdf::ωa_vers_ωd(ωa, fe);
}
inline float ωd2ωa(float wd, float fe)
{
  return tsdf::ωd_vers_ωa(wd, fe);
}

/** @brief Bilinear transform frequency warping (digital to analog).
 *
 * <h3>Bilinear transform frequency warping (digital to analog)</h3>
 *
 *
 *  @param fd Digital frequency (between 0 et 0.5)
 *  @return   Analog frequency (unbounded, positive)
 *
 *  The analog frequency is computed as:
 *  @f[
 *  f_a = \frac{\tan(\pi  f_d)}{\pi}
 *  @f]
 *
 *  This mapping is the one caused by the bilinear transform.
 *  */
inline float fd2fa(float fd)
{
  return tsdf::fd_vers_fa(fd);
}

/** @brief Bilinear transform frequency warping (analog to digital).
 *
 *  <h3>Bilinear transform frequency warping (analog to digital)</h3>
 *
 *  @param fa Analog frequency (unbounded, positive)
 *  @return   Digital frequency (between 0 et 0.5)
 *
 * Inverse function of de @ref fd2fa() */
inline float fa2fd(float fa)
{
  return tsdf::fa_vers_fd(fa);
}

/** @} */


/**  @addtogroup filtrage-tr
  *  @{ */


/** @brief Delay the input signal by an integer number of samples.
 *
 *  <h3>Delay line</h3>
 *
 *  This filter produce as many output samples as input samples.
 *  The samples preceding the first one are assumed to be null.
 *
 *  @param n Integer delay (@f$n \geq 0@f$)
 *
 *  @par Example
 *  @snippet exemples-tsd-pub/src/ex-filtrage.cc exemple_ligne_a_retard
 *  @image html filtrage-ligne-a-retard.png width=800px
 *  */
template<typename T>
  sptr<FiltreGen<T>> delay_line(unsigned int n)
{
  return tsdf::ligne_a_retard<T>(n);
}

using tsdf::HilbertTransformeurConfig;

/** @brief Defines a Hilbert transformer (conversion from a real signal to a complex, analytic one.
 *
 * <h3>Hilbert transformer</h3>
 *
 * This filter computes an analytic (complex) signal from a real one,
 * by recomposition of the input signal (delayed) with the signal filtered with a Hilbert filter:
 *
 * @image  html fig-hilbert-transfo.PNG width=200px
 *
 * @param  n        Filter order
 * @param  fenetre  Window choice (see @ref window())
 * @return          Filter float -> cfloat
 *
 * @sa design_fir_hilbert(), hilbert(), hilbert_dft()
 */
inline sptr<Filtre<float, cfloat, HilbertTransformeurConfig>>
  hilbert_transformer(unsigned int n = 31, const std::string &fenetre = "hn")
{
  return tsdf::hilbert_transformeur(n, fenetre);
}

/** @brief Direct implementation of a FIR filter.
 *
 * <h3>Direct implementation of a FIR filter</h3>
 *
 *  Implemented with the convolution equation:
 *  @f[
      y_n = \sum h_k x_{n-k}
    @f]

 *  @param h      Vector with filter coefficients.
 *  @tparam T     Scalar data type (float, cfloat, etc.).
 *  @tparam Tc    Coefficients type (float, cfloat, etc.).
 *  @return       Filter T -> T
 *
 *  @note  If the number of coefficients is important,
 *  better usr the @ref filter_fir_fft(), which may be more efficient
 *  (filtering in the frequency domain).
 *
 *  @sa filter_fir_fft()
 */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filter_fir(const Eigen::Ref<const Vector<Tc>> h)
{
  return tsdf::filtre_rif(h);
}


/** @brief Identity filter.
 *
 *  <h3>Identity filter</h3>
 *
 *  This filter let the signal unchanged.
 */
template<typename T>
  sptr<FiltreGen<T>> filter_id()
{
  return tsdf::filtre_id();
}

/** @brief Decimation by a factor 1:R
 *
 *  <h3>Decimation by a factor 1:R</h3>
 *
 *  This "filter" delete @f$R-1@f$ samples every @f$R@f$ samples.
 */
template<typename T>
  sptr<FiltreGen<T>> decimator(int R)
{
  return tsdf::decimateur(R);
}




/** @brief Efficient implementation of a FIR filter in the frequency domain.
 *
 *  <h3>Efficient implementation of a FIR filter in the frequency domain</h3>
 *
 *  This filter implements a FIR filter through the (Ovelap-And-Add) technique.
 *  The complexity is thus much lower than the direct implementation if the number of coefficients @f$M@f$
 *  is important (in the order of @f$\log M@f$ operations by sample instead of @f$M@f$).
 *
 *  @param h Vector with filter coefficients.
 *  @tparam T Scalar type of input / output samples.
 *  @return Filter T -> T
 *
 *  @note Note that this technique introduce a little more delay than the direct, time domain, implementation.
 *
 *  @sa filter_fir()
 */
template<typename T>
  sptr<FiltreGen<T>> filter_fir_fft(const ArrayXf &h)
{
  return tsdf::filtre_rif_fft(h);
}



/** @brief IIR filter, direct form I implementation, not  recommanded (use rather @ref filter_sois() instead).
 *
 *  <h3>IIR filter, direct form I implementation</h3>
 *
 *  This filmter implements an IIR filter, using the most direct form (direct I),
 *  that is, the filter is decomposed as:
 *  @f[
 *  H(z) = \frac{b_0 + b_1 z^{-1} + \dots}{a_0 + a_1 z^{-1} + \dots}
 *       = \left(b_0 + b_1 z^{-1} + \dots\right) \cdot \frac{1}{a_0 + a_1 z^{-1} + \dots}
 *  @f]
 *  and the FIR filter corresponding to the numerator is computed first, then the all-poles RII filter corresponding to the denominator.
 *
 *  @param h    Transfert function (rational function)
 *  @tparam T   Scalar type of input / output samples.
 *  @tparam Tc  Coefficients type.
 *  @return     Filter T -> T
 *
 *  @sa filter_sois()
 *
 *  @warning If the filter order is large, because of truncation errors,
 *  this implementation a a lot of chance to diverge.
 *  The alternative implementation through a cascade of second order RII filters (@ref filter_sois()) is then recommanded.
 *
 * */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filter_iir(const FRat<Tc> &h)
{
  return tsdf::filtre_rii(h);
}


/** @brief Creation of a CIC filter, processing on vectors of type T,
 *  and computing internally with samples of type Ti.
 *
 *  <h3>Creation of a CIC filter</h3>
 *
    @param config Main parameters of the CIC filter (see @ref CICConfig)
 *  @param mode 'd' for decimation or 'u' for upsampling.
 *  @tparam T   Input / output type.
 *  @tparam Ti  Internal sample type (int, int64_t, ...)
 *
 *  @warning For the internal sample type,
 *  it is highly recommanded to choose an integer type,
 *  because the way the filter is implemented can make the filter
 *  diverge with flotting point arithmetic.
 *
 *
 *  @par Example for interpolation
 *  @snippet exemples/src/filtrage/ex-cic.cc exemple_cic_upsampling
 *  @image html filtrage-cic-interpolation.png width=800px
 *  Note the frequency aliasings.
 *  <br>
 *  @par Example for decimation
 *  @snippet exemples/src/filtrage/ex-cic.cc exemple_cic_decimation
 *  @image html filtrage-cic-decimation.png width=800px
 */
template<typename T, typename Ti>
  sptr<FiltreGen<T>> filter_cic(const CICConfig &config, char mode = 'd')
{
  return tsdf::filtre_cic(config, mode);
}

using tsdf::CICAnalyse;

/** @brief This function computes and shows the frequency response of a CIC filter and
 *  then analyse the aliasing that occurs after decimation.
 *
 *  @param config Main parameters of the CIC filter (see @ref CICConfig)
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
inline CICAnalyse cic_analysis(const CICConfig &config, float Fin, float Fint = 0)
{
  return tsdf::cic_analyse(config, Fin, Fint);
}



/** @brief Frequency response of a CIC filter
 *
 * @param config Main parameters of the CIC filter (see @ref CICConfig)
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
inline ArrayXf cic_freq(const CICConfig &config, const ArrayXf &f)
{
  return tsdf::cic_freq(config, f);
}




using RIIStructure = tsd::filtrage::RIIStructure;


/** @brief IIR filter implementation through a cascad of second order sections.
 *
 *  <h3>IIR filter implementation through a cascad of second order sections</h3>
 *
 *  The given transfert function is factored into a cascad
 *  of second order sections, and eventually a first order section if the whole filter order is odd,
 *  enabling an efficient implementation of a IIR filter:
 *
 *  @f[
 *  H(z) = G\cdot\prod_{i=1}^{M} H_i(z) = G\cdot\prod_{i=1}^{M} \frac{1 + b_1^{(i)} z^{-1} + b_2^{(i)} z^{-2} }{1 - a_1^{(i)} z^{-1} - a_2^{(i)} z^{-2}}
 *  @f]
 *
 *  @param h Transfert function.
 *  @param structure RIIStructure::FormeDirecte1 or RIIStructure::FormeDirecte2 (see notes below)
 *  @return real to real filter.
 *
 *  The direct form II is a little more efficient than the form I. However, if the coefficients are to be reloaded while the filter is in use,
 *  the form II can generate discontinuities, which is not the case with the first form.
 *
 *  @warning It is highly recommanded that the given transfert function is specified by its roots and poles rather by its coefficients,
 *  because in the latter case, the
 *  poles and zeros locations (and hence the frequency response) may be very unstable (high sensitivity to troncating errors).
 *  Note that the function @ref design_iira() returns a correct tranfert function (zeros & poles).
 *
 *  @sa filter_iir(), design_iira()
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

/** @brief First order IIR filter (exponential filter).
 *
 *  <h3>First order IIR filter (exponential filter)</h3>
 *
 *  This filter, also called the "numerical RC filter", or "exponential filter",
 *  is one of the simplest filter, as it is completely speficied by only one coefficient.
 *
 *  @param γ Forget factor.
 *
 *  It is defined by the following equation:
 *  @f[
 *   y_{n+1} = (1-\gamma) * y_n + \gamma * x_n = y_n + \gamma * (x_n - y_n)
 *  @f]
 *
 *  The @f$\gamma@f$ parameter (forget factor) can be tuned easily as a function of time constant or of
 *  the desired cut-off frequency (see @ref iir1_coef()).
 *
 *  @sa iir1_fcut(), iir1_coef()
 */
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






