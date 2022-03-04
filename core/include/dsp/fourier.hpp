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
 * @return            A poin @ref Filter cfloat vers cfloat.
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
inline sptr<FilterGen<float, cfloat>> rfftplan_new(int n = -1)
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
 *  For different purposes, it can be convenient to apply filtering in the frequency domain rather
 *  than in the time domain, for instante to reduce the computing load.
 *  (time domain convolutions become simple term by term product in the frequency domain).
 *  By frequency domain filtering, we mean here to be able to adjust term by term the different bins of the DFT of a signal,
 *  and then to build the result in the time domain.
 *
 *  So as to work in the frequency domain, samples are grouped into blocks of
 *  @f$N_e@f$ samples, and for each of these blocks the following operations are performed:
 *   - <b>Windowing and overlap</b> (optional): so as to smooth processings from one block to the next one,
 *      an apodisation window is applied (by default,
 *      a Hann window with an overlap of 50 %, so as to enable perfect reconstruction of the time domain signal).
 *   - <b>Zero-padding</b> : So as to compute efficiently the FFT,
 *   zeros are inserted, making the length of each block equal to some power of 2.
 *   So, after zero-padding, the blocks length is @f$N=2^k\geq N_e@f$.
 *   - <b>Fourier transform:</b> through a FFT.
 *   - <b>Frequency domain processing</b> : through a callback provided by the caller (see below).
 *   - <b>Inverse Fourier transform</b>.
 *   - <b>Recomposition of time domain blocks</b>: through the OLA technique (<i>OverLap-and-Add</i>).
 *
 *  @param config Configuration structure (@ref FFTFilterConfig) which specify:
 *   - Dimension @f$N_e@f$ of input time-domain blocks,
 *   - Minimal number of zeros to be inserted at the beginnning of each block,
 *   - User callback for frequency domain processing.
 *
 *  @return A tuple with two elements:
 *    - the OLA filter,
 *    - an integer which is the dimension @f$N@f$ of the FFT (that is,
 *      the dimension of the vectors that will be given to the callback).
 *
 *
 *  @parblock
 *  @note For the more specific problem of efficient FIR filter computing (using the FFT), use rather
 *  @ref filter_fir_fft(), which is based on this function, and is simpler to use.
 *  @endparblock
 *
 *  @par Example: Low-pass filter
 *  In this example, we low-pass filter a signal by simply zeroing half the spectrum.
 *  @snippet exemples/src/fourier/ex-ola.cc ex_ola_pb
 *  @image html ex-ola-pb.png width=800px
 *
 *
 *  @sa filter_fir_fft(), filter_rfft(), ola_complexity(), ola_complexity_optimize()
 */
inline std::tuple<sptr<Filter<cfloat, cfloat, tsdF::FiltreFFTConfig>>, int> filter_fft
  (const FFTFilterConfig &config)
{
  return tsdF::filtre_fft(config);
}



/** @brief Compute the complexity of an OLA filter, in FLOPS by input sample.
 *
 * <h3>Complexity of an OLA filter</h3>
 *
 *  @param M    Time-domain filter length (in samples),
 *  @param Ne   Input block length,
 *  @param[out] C  Complexity, in FLOPS by input sample,
 *  @param[out] Nf FFT dimension,
 *  @param[out] Nz Number of zeros inserted for each FFT block (@f$N_f = N_e + N_z@f$)
 *
 *  @sa ola_complexity_optimize()
 */
inline void ola_complexity(int M, int Ne, float &C, int &Nf, int &Nz)
{
  return tsdF::ola_complexité(M, Ne, C, Nf, Nz);
}

/** @brief Compute optimal parameters for an OLA filter.
 *
 * <h3>Optimal parameters for an OLA filter</h3>
 *
 * This function computes the optimal input block length for an OLA process (see @ref filter_fft()).
 * The optimum is found for @f$N_e = 2^k - (M-1)@f$, @f$k@f$
 * being found by exhaustive search.
 *
 *  @param M  Time-domain filter length (in samples),
 *  @param[out] C  Complexity, in FLOPS by input sample,
 *  @param[out] Nf FFT dimension,
 *  @param[out] Nz Number of zeros inserted for each FFT block (@f$N_f = N_e + N_z@f$)
 *  @param[out] Ne Optimal input block length.
 *
 *  @sa ola_complexity()
 *
 */
inline void ola_complexity_optimize(int M, float &C, int &Nf, int &Nz, int &Ne)
{
  tsdF::ola_complexité_optimise(M, C, Nf, Nz, Ne);
}


/** @brief Z-chirp transform.
 *
 *  <h3>Z-chirp transform</h3>
 *
 *  Evaluate the z transform on the following points:
 *    @f$z_0 \cdot W^n@f$, pour @f$n=0,1,...,m-1@f$
 *
 *  @param x  Input vector
 *  @param m  Number of evaluation points
 *  @param W  Progression of the geometric serie
 *  @param z0 Initial term
 *  @returns The values of @f$H(z)@f$:
 *  @f[
 *    X\left(z_0 * W^n\right) = \sum_{k=0}^{N-1} x_k \cdot z_0^k \cdot W^{nk},\quad \textrm{pour }n=0,1,...,m-1
 *  @f]
 */
inline ArrayXcf czt(IArrayXcf x, int m, cfloat W, cfloat z0 = 1.0f)
{
  return tsdF::czt(x, m, W, z0);
}

/** @} */

/** @addtogroup fourier-corr
 *  @{
 */

/** @brief Circular correlation (normalized) between two complex vectors.
 *
 *  <h3>Circular correlation product</h3>
 *
 *  FFT-based computing of:
 * @f[
 * c_n = \frac{1}{N} \cdot \sum_{k=0}^{N-1} x_k y^{\star}_{k+n[N]},\ \ n=0\dots N-1
 * @f]
 *
 *  @param x First vector
 *  @param y Second vector (if not specified, the auto-correlation is computed)
 *  @returns First vector: index of time lags @f$n@f$ (that is @f$0, 1, ..., N-1@f$), second vector: @f$c_n@f$.
 *
 *
 * @sa xcorr(), xcoorb()
 */
inline auto ccorr(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf())
{
  return tsdF::ccorr(x, y);
}


/** @brief Unbiased correlation between two complex vectors
 *
 *  <h3>Correlation product (with bias correction)</h3>
 *
 *  Compute the correlation between 2 complex vectors (through FFT),
 *  that is, for negative lags (@f$\tau=-(m-1)\dots -1@f$) :
 *  @f[
 *  c(\tau) = \frac{1}{N-\tau} \cdot \sum_{k=0}^{N-1-\tau} x_{k-\tau} y^{\star}_k
 *  @f]
 *  And, for positive lages (@f$\tau=0\dots m-1@f$):
 *  @f[
 *  c(\tau) = \frac{1}{N-\tau} \cdot \sum_{k=0}^{N-1-\tau} x_k y^{\star}_{k+\tau}
 *  @f]
 *
 *  Lags are sampled uniformly between @f$-(m-1)@f$ et @f$m-1@f$ :
 *  @f[
 *   \tau = \textrm{linspace}(-(m-1), m-1, 2m-1) = -(m-1), -(m-2), \dots, -1, 0, 1, \dots, (m-2), (m-1)
 *  @f]
 *
 *  @param x First vector
 *  @param y Second vector (if not specified, the auto-correlation is computed)
 *  @param m Number of delays to be computed (if negative, then take @f$m=N@f$).
 *  @returns First vector: index of time lags @f$n@f$ (that is @f$0, 1, ..., N-1@f$), second vector: @f$c_n@f$.
 *
 *
 * @note The vectors are implicitly zero-padded so as to be able to compute efficiently the correlation in the frequency domain.
 *
 * @par Example
 * @snippet exemples/src/fourier/ex-fourier.cc ex_xcorr
 * @image html xcorr.png width=800px
 *
 * @sa xcorrb(), ccorr(), detector_new()
 */
inline auto xcorr(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf(), int m = -1)
{
  return tsdF::xcorr(x, y, m);
}


/** @brief Correlation (biased) between two complex vectors.
 *
 *  <h3>Correlation product (without bias correction)</h3>
 *
 *  Computes (through FFT):
 *  @f[
 *  c_n = \frac{1}{N} \cdot \sum_{k=0}^{N-1-k} x_k y^{\star}_{k+n},\ \ n=-m\dots m-1
 *  @f]
 *
 *  @param x First vector
 *  @param y Second vector (if not specified, the auto-correlation is computed)
 *  @param m Number of delays to be computed (if negative, then take @f$m=N@f$).
 *  @returns First vector: index of time lags @f$n@f$ (that is @f$0, 1, ..., N-1@f$), second vector: @f$c_n@f$.
 *
 *  @sa xcorr(), ccorr(), detector_new()
 */
inline auto xcorrb(const ArrayXcf &x, const ArrayXcf &y = ArrayXcf(), int m = -1)
{
  return tsdF::xcorrb(x, y, m);
}

/** @brief FFT-based delaying of a vector.
 *
 *  <h3>FFT-based delaying of a vector</h3>
 *
 *  Apply a fractionnal delay to the input signal:
 *  @f[
 *  y_k = x_{k-\tau}
 *  @f]
 *
 *
 *  The delay can be positive (to effectively delay the signal), or negative
 *  (to have the signal earlier).
 *
 *  @note If the request delay is an integer, a simple shift of the vector is done.
 *  Otherwise, the delay is applied in the frequency domain
 *  (after the input signal has been zero-padded so as to mitigate border effects):
 *  @f[
 *  Y_n = X_n \cdot e^{\frac{-2\pi\mathbf{i}\tau \cdot n}{N}}
 *  @f]
 *
 *  @param x  Vector of samples.
 *  @param τ  Delay to be applied, in number of samples (can be positive or negative, and not necessarily an integer).
 *  @returns  Time shifted vector.
 *
 *
 * @par Example
 * @snippet exemples/src/fourier/ex-fourier.cc ex_delais
 * @image html ex-delais.png width=800px
 *
 **/
template<typename T = float>
  Vector<T> delay(const Vector<T> &x, float τ)
{
  return tsdF::délais(x, τ);
}


/** @brief Delay estimation between two vectors. */
inline std::tuple<float, float> delay_estimation(IArrayXcf x, IArrayXcf y)
{
  return tsdF::estimation_délais(x, y);
}

/** @brief Alignement of two vectors */
template<typename T>
  std::tuple<Vector<T>, Vector<T>, int, float> align_int(const Vector<T> &x, const Vector<T> &y)
{
  return tsdF::aligne_entier(x, y);
}


/** @brief Informations computed from the detected pattern. */
struct Detection
{
  /** @brief Position in number of samples (between @f$-N_e@f$ and @f$N_e-1@f$),
  *    since the beginning of current block. For instance:
  *    - @f$0\  \Leftrightarrow@f$    First sample of current block
  *    - @f$1\  \Leftrightarrow@f$    Second sample of current block
  *    - @f$-1\  \Leftrightarrow@f$   Last sample of previous block
  *    - ... */
  int position;

  /** @brief Same a position field, but using quadratic interpolation for sub-sample accuracy. */
  float position_frac;

  /** @brief Absolute value of the normalized correlation (positive value, between 0 et 1). */
  float score;

  /** @brief Signal gain (comparing to the expected pattern). */
  float gain;

  /** @brief Signal phase (between @f$-\pi/2@f$ and @f$\pi/2@f$) */
  float θ;

  /** @brief Estimated SNR. */
  float SNR_dB;

  /** @brief Estimated standard deviation of the noise. */
  float σ_noise;

  Detection(const tsdF::Detection &fr)
  {
    memcpy(this, &fr, sizeof(tsdF::Detection));
  }
};


/** @brief Configuration structure for a FFT based detector. */
struct DetectorConfig: tsdF::DetecteurConfig
{
  /* @brief Minimum number of samples for each time-domain block (if 0 : computed automatically). */
  uint32_t &Ns        = Ne;

  /* @brief Pattern to be detected. */
  ArrayXcf &pattern   = motif;

  /* @brief Detection threshold, between 0 et 1. */
  float &threshold    = seuil;

  bool &debug_active  = debug_actif;

  /** @brief Computing mode of the correlation product: FFT with OLA or simple FIR filter. */
  tsdF::DetecteurConfig::Mode &mode = mode;

  /** @brief User callback called each time the pattern is detected. */
  std::function<void (const Detection &det)> &on_detection
      = *((std::function<void (const Detection &det)> *) &gere_detection);

  bool &compute_correlation_signal = calculer_signal_correlation;
};



using tsdF::Detecteur;

/** @brief Correlation-based pattern detector.
 *
 * <h3>Correlation-based pattern detector</h3>
 *
 * This block will compute the normalized correlation between a streaming signal @f$(x_k)@f$
 * and a fixed pattern @f$(h_k)@f$ of dimension @f$M@f$ :
 * @f[
 * y_n = \frac{\left|\displaystyle{\sum_{k=0}^{M-1} x_{n+k-(M-1)} \cdot h_{k}^\star}\right|}{\displaystyle{ \sqrt{\left(\sum_{k=0}^{M-1}\left|x_{n+k-(M-1)}\right|^2\right) \left(\sum_{k=0}^{M-1}\left|h_k\right|^2\right)}}}
 * @f]
 *
 * The @f$y_n@f$ are thus between 0 and 1, and will be equal to 1 only if the signal is exactly equal to the pattern (up to some scale factor)
 * on the interval @f$[n-M+1\dots n]@f$.
 *
 * The computing being done in the frequency domain through the Overlap-And-Add technique
 * (see @ref filter_fft()), the complexity is of order @f$\log_2 M@f$ operations
 * per sample (if @f$M@f$ is a power of 2).
 *
 * This block will:
 *  - Return, as output, the correlation signal (green signal in the example below),
 *  - Call a user callback whenever the  correlation is above a configurable threshold.
 *  As parameters of this callback,
 *    the following informations are given (see the structure @ref Detection):
 *      - Peak value of the normalized correlation,
 *      - Real gain of the sinal (amplitude ratio between the received signal and the expected pattern),
 *      - Relative phase of the received signal (compared to the expected pattern),
 *      - Noise standard deviation,
 *      - Estimated SNR.
 *
 * @param config Configuration structure
 * @return %Filter cfloat @f$\to@f$ float, that take as input samples, and produces as output,
 * the normalized correlation with the pattern.
 *
 * @par Example
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

/** @brief PSD (correlogram).
 *
 * <h3>PSD (correlogram)</h3>
 *
 * Computes the correlogram using a Hann window:
 *
 * @f[
 * S(k) = \left|TFD(x \cdot w)(k)\right|
 * @f]
 *
 * The resulting spectrum is 0 Hz centered (frequency vector: linspace(-0.5,0.5-1.0/N,N)).
 *
 * @param x   Signal to be analyzed.
 * @return    A tuple with two elements: frequences vectors (normalized, between -0.5 and 0.5), and the spectrum.
 *
 * @par Example
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

/** @brief PSD (Welch method - spectrum averaging).
 *
 * <h3> PSD using Welch method</h3>
 *
 * This function estimates the PSD using the spectrum averaging method:
 * several PSD are computed using windowed parts of the signal (optionnaly with some overlap), then averaged.
 *
 * @param x   Signal to be analyzed.
 * @param N   Dimension of each window (that is, the resulting PSD resolution), in number of samples.
 * @param fen Type of windowing ("hn", "hm", "re", ...).
 * @return A tuple with two elements: frequences vectors (normalized, between -0.5 and 0.5), and the spectrum.
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
 *  @param x  Signal to be analyzed.
 *  @param Ns Nombre d'exponentielles pures attendues
 *  @param Nf Résolution fréquentielle (en nombre de points)
 *  @param m  Dimension de la matrice d'auto-corrélation (si non précisé = dimension du signal)
 *  @returns A tuple with two elements: frequences vectors (normalized, between -0.5 and 0.5), and the spectrum.
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
  return tsdF::psd_sousesp(x, Ns, Nf, m);
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
inline sptr<FilterGen<float>> filter_goertzel(float frequence, int N)
{
  return tsdF::filtre_goertzel(frequence, N);
}




/** @brief Configuration pour le calcul de spectre en temps réel */
using tsdF::SpectrumConfig;



/** @brief Calcul de spectre en temps réel */
inline sptr<Filter<cfloat,float,SpectrumConfig>> rt_spectrum(const SpectrumConfig &config)
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



