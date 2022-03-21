#pragma once

/** (C) 2022 J. Arzi / GPL V3 - see LICENSE file. */

#include <dsp/view.hpp>
#include "dsp/dsp.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "dsp/filter.hpp"
#include "dsp/fourier.hpp"
#include "tsd/telecom.hpp"
#include <vector>
#include <random>


namespace dsp::telecom {


  using tsd::MoniteursStats;


namespace nfr = tsd::telecom;

using nfr::BitStream; // TODO

/** @addtogroup telecom-ps
 *  @{
 */


/** @brief Specification of a shaping filter */
struct ShapingFilterSpec
{
  ShapingFilterSpec(const nfr::SpecFiltreMiseEnForme &fr)
  {
    memcpy(this, &fr, sizeof(*this));
  }

  auto fr() const
  {
    nfr::SpecFiltreMiseEnForme res;
    memcpy((void *) &res, this, sizeof(*this));
    return res;
  }


  /** @brief No filtering, raw impulses are transmitted. */
  static ShapingFilterSpec none()
  {
    return nfr::SpecFiltreMiseEnForme::aucun();
  }

  /** @brief Gaussian + moving average filter. */
  static ShapingFilterSpec gaussian(float BT)
  {
    return nfr::SpecFiltreMiseEnForme::gaussien(BT);
  }

  /** @brief "NRZ" filtering (moving average). */
  static ShapingFilterSpec nrz()
  {
    return nfr::SpecFiltreMiseEnForme::nrz();
  }

  /** @brief SRRC filtering (Square-Root Raised Cosine). */
  static ShapingFilterSpec srrc(float β)
  {
    return nfr::SpecFiltreMiseEnForme::rcs(β);
  }

  /** @brief Filter type. */
  enum Type
  {
    /** @brief "NRZ" filtering (moving average). */
    NRZ,
    /** @brief No filtering, raw impulses are transmitted. */
    NONE,
    /** @brief Gaussian + moving average filter. */
    GAUSSIAN,
    /** @brief Square root raised cosine filter. */
    SRRC
  };

  /** @brief Filter type. */
  Type type = Type::SRRC;

  /** @brief BT product (Gaussian filter only). */
  float BT            = 0.8;

  /** @brief Roll-off factor for SRRC filter (see @ref design_fir_rc()) */
  float β = 0.2;

  /** @brief Computes the coefficient of FIR shaping filter.
   *
   *  <h3>Coefficient of FIR shaping filter</h3>
   *
   *  @param ncoefs Desired number of taps.
   *  @param osf    Oversampling factor (e.g. numbers of samples per symbol).
   *  @returns      Coefficient vector.
   */
  ArrayXf get_coefs(int ncoefs, int osf) const
  {
    return fr().get_coefs(ncoefs, osf);
  }


  /** @brief Creation of a shaping filter with integrated over-sampling.
   *
   *  <h3>Shaping filter</h3>
   *
   *  Creation of a shaping filter with integrated over-sampling.
   *  This filter block accepts as input symbols,
   *  and generate an oversampled (by a factor of @f$R@f$, as given in parameter)
   *  and pulse shaped signal.
   *
   *  The filter coefficients are energy-normalized:
   *  @f[
   *  \sum h_k^2 = R
   *  @f]
   *  so as the mean power at the input and output of the filter are identical.
   *
   *  @param ncoefs Number of taps for RIF filter.
   *  @param R      Oversampling ratio.
   *
   *  @sa matched_filter()
   */
  sptr<FilterGen<cfloat>> shaping_filter(int ncoefs, int R) const
  {
    return fr().filtre_mise_en_forme(ncoefs, R);
  }


  /** @brief Same as @ref shaping_filter() method,  but without upsampling.
   *
   *  <h3>Matched filter</h3>
   *
   *  @param ncoefs Number of taps for RIF filter.
   *  @param osf    Oversampling ratio.
   *
   *  @sa shaping_filter()
   */
  sptr<FilterGen<cfloat>> matched_filter(int ncoefs, int osf) const
  {
    return fr().filtre_adapté(ncoefs, osf);
  }

  /** @brief Matched filter and downsampling at symbol rate.
   *
   *  <h3>Matched filter and downsampling at symbol rate</h3>
   *
   *  @param ncoefs Number of taps for RIF filter.
   *  @param osf    Input oversampling ratio.
   *
   *  @sa shaping_filter(), matched_filter()
   */
  sptr<FilterGen<cfloat>> matched_filter_with_decimation(int ncoefs, int osf) const
  {
    return fr().filtre_adapté_décimation(ncoefs, osf);
  }

  struct Analyse
  {
    tsd::vue::Figures fig;
  };

  Analyse analyse(int ncoefs, int osf) const;

};

inline std::ostream& operator<<(std::ostream &ss, const ShapingFilterSpec &t)
{
  return ss << t.fr();
}



/** @}
 *  @addtogroup telecom-mods-wf
 *  @{
 */





/** @brief WaveForm specification. */
struct WaveForm
{
  sptr<nfr::FormeOnde> fr;

  WaveForm(sptr<nfr::FormeOnde> fo)
  {
    this->fr = fo;
  }

  /** @brief Generation of I/Q symbols from a binary stream. */
  ArrayXcf make_symbols(const BitStream &bs)
  {
    return fr->génère_symboles(bs);
  }

  /** @brief Generation of I/Q samples from a binary stream (including upsampling and shaping filter).
   *
   * <h3>Generation of I/Q samples from a binary stream</h3>
   *
   * This function will create, from a binary sequence, I/Q samples.
   * The following steps are done:
   *   - Constellation mapping,
   *   - Upsampling,
   *   - Shaping filter.
   *
   * @param bs            Binary sequence.
   * @param ncoefs        Number of taps for shaping filter implementation.
   * @param osf           Oversampling factor.
   * @param[out] delay    Delay, in number of output samples, between the beginning of the output I/Q stream, and the middle of the first symbol (the delay is due to the shaping filter).
   *
   *
   */
  ArrayXcf make_samples(const BitStream &bs, int ntaps, int osf, float &delay)
  {
    return fr->génère_échantillons(bs, ntaps, osf, delay);
  }

# if 0
  // Contexte de démodulation,
  // pour les modulations à mémoire (par exemple, FSK, π/4-QPSK).
  // c'est-à-dire où la constellation n'est pas constante.
  struct Ctx
  {
    virtual void reset() = 0;

    /** Index = -1 si pas d'échantillon à sortir */
    virtual std::tuple<int, cfloat> step(cfloat x) = 0;
  };

  virtual sptr<Ctx> get_ctx(int OSF = 1) const; // Par défaut, contexte sans mémoire

  // Contexte pour la génération de symboles
  struct CtxGen
  {
    virtual void reset() = 0;
    virtual ArrayXcf step(const BitStream &bs) = 0;
  };
# endif

  sptr<nfr::FormeOnde::CtxGen> get_contexte_tx(int ncoefs, int osf)
  {
    return fr->get_contexte_tx(ncoefs, osf);
  }

  /** @brief Decoding of I/Q symbols to binary sequence. */
  void decode_symbols(BitStream &bs, const ArrayXcf &x)
  {
    fr->decode_symboles(bs, x);
  }

  /** @brief Returns the symbol number i from the constellation, as a complex value. */
  cfloat get_symbol(unsigned int i) const
  {
    return fr->lis_symbole(i);
  }

  /** @brief Returns the index to the closest symbol in the constellation. */
  int closest_symbol(const cfloat &point) const
  {
    return fr->symbole_plus_proche(point);
  }

  /** @brief Theoretical binary error rate (ber), for this waveform, as a function of the normalized SNR. */
  float ber(float EbN0_dB)
  {
    return fr->ber(EbN0_dB);
  }

  /** @brief Theoretical binary error rate (ber), for this waveform, as a function of the normalized SNR. */
  ArrayXf ber(const ArrayXf &EbN0_dB)
  {
    return fr->ber(EbN0_dB);
  }

  /** @brief Returns all the points of the constellation. */
  ArrayXcf constellation() const
  {
    return fr->constellation();
  }

  /** @brief Get frequency excursion fréquentielle, in multiple of symbol frequency. */
  float excursion() const
  {
    return fr->excursion();
  }

  /** @brief A short description of this waveform. */
  std::string desc() const
  {
    return fr->desc();
  }

  /** @brief Miscelaneous informations about this waveform. */
  struct Infos
  {
    /** @brief True for modulation such as PSK, ASK, QAM, false for FSK. */
    bool is_linear = true;

    /** @brief True for linear and constant amplitude modulation. */
    bool is_psk = false;

    /** @brief Amplitude only modulation. */
    bool is_ask = false;

    /** @brief Frequency modulation. */
    bool is_fsk = false;

    /** @brief Quadrature modulation (phase & amplitude). */
    bool is_qam = false;

    /** @brief Modulation index (valid only for FSK modulations) */
    float index = 1.0f;

    /** @brief Total number of possible symbols. */
    int M;

    /** @brief Number of bits conveyed by each symbole (@f$k=\log_2(M)@f$) */
    int k;
  };

  /** @brief Get waveform informations. */
  Infos get_infos() const
  {
    Infos res;
    memcpy((void *) &res, &(fr->infos), sizeof(Infos));
    return res;
  }

  /** @brief Change the shaping filter. */
  void set_shaping_filter(const ShapingFilterSpec &spec)
  {
    fr->filtre = spec.fr();
  }

  /** @brief Read shaping filter specifications. */
  ShapingFilterSpec get_shaping_filter() const
  {
    return fr->filtre;
  }
};

inline std::ostream& operator<<(std::ostream &ss, const WaveForm &t)
{
  return ss << t.fr;
}

/** @brief Phase modulation.
 *
 * <h3>Phase modulation</h3>
 *
 * Creation of M-PSK waveform.
 *
 * @param   M         Number of bits / symbol.
 * @param   filtre    Shaping filter (default: NRZ)
 * @return            Abstract waveform pointer (@ref WaveForm).
 *
 * @warning @f$M@f$ must be a power of 2.
 *
 * @par Example 1: Some PSK modulations
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_psk
 * @image html waveform-psk.png "Waveforms examples: BPSK, QPSK, 8PSK, 16PSK" width=800px
 *
 * @par Example 2: Theoretical bit error rates
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_psk2
 * @image html waveform-psk2.png "BER" width=800px
 *
 * @sa waveform_qam(), waveform_qpsk()
 */
inline sptr<WaveForm> waveform_psk(unsigned int M, const ShapingFilterSpec &filtre = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_psk(M, filtre.fr()));
}

/** @brief Creation of a BPSK waveform.
 *
 * <h3>Creation of a BPSK waveform</h3>
 *
 * Same as @ref waveform_psk() for @f$M = 2@f$.
 *
 * @sa waveform_psk(), waveform_qam(), waveform_fsk()
 */
inline sptr<WaveForm> waveform_bpsk(const ShapingFilterSpec &filtre = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_bpsk(filtre.fr()));
}

/** @brief Creation of a M-ASK waveform.
 *
 * <h3>Creation of a M-ASK waveform</h3>
 *
 * @f[
 * x_n = K_1 + \frac{s_n}{M-1}\cdot K_2
 * @f]
 *
 *
 * @sa waveform_psk(), waveform_bpsk(), waveform_qam(), waveform_fsk()
 */
inline sptr<WaveForm> waveform_ask(int M = 2, float K1 = -1, float K2 = 2, const ShapingFilterSpec &filter = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_ask(M, M, M, filter.fr()));
}

/** @brief Creation of a QPSK waveform.
 *
 * <h3>Creation of a QPSK waveform</h3>
 *
 * Same as @ref waveform_psk() for @f$M = 4@f$.
 *
 * @sa waveform_psk(), waveform_qam(), waveform_fsk()
 */
inline sptr<WaveForm> waveform_qpsk(const ShapingFilterSpec &filter = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_qpsk(filter.fr()));
}

/** @brief Creation of a π/4 - QPSK waveform.
 *
 * <h3>Creation of a π/4 - QPSK waveform</h3>
 *
 * Alternate between two QPSK constellations rotated by @f$\pi/4@f$.
 *
 * @sa waveform_qpsk()
 */
inline sptr<WaveForm> waveform_π4_qpsk(const ShapingFilterSpec &filter = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_π4_qpsk(filter.fr()));
}

/** @brief Creation of a QAM waveform.
 *
 * <h3>Creation of a QAM waveform</h3>
 *
 * Constellation points are put on a regular grid.
 *
 * @par Example: QAM constellations plotting
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_qam
 * @image html waveform-qam.png "QAM16, QAM64, QAM256 constellations" width=800px
 *
 * @sa waveform_psk(), waveform_fsk()
 */
inline sptr<WaveForm> waveform_qam(unsigned int M, const ShapingFilterSpec &filter = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_qam(M, filter.fr()));
}

/** @brief Creation of a FSK waveform.
 *
 * <h3>Creation of a FSK waveform</h3>
 *
 * FSK is for Frequency Shift Keying.
 * The main parameter of this modulation is the <b>modulation index</b>,
 * which is the ratio between the maximal excursion (that is, twice the deviation) and the symbol rate:
 * @f[
 * h = \frac{2 \Delta f}{f_{symb}}
 * @f]
 *
 * the instantaneous frequency varying between @f$f_c - \Delta f@f$ and @f$f_c + \Delta f@f$.
 *
 * @param M       Number of possible values for each symbol.
 * @param index   Modulation index.
 * @param filtre  Shaping filter.
 *
 * @par Example
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_fsk
 * @image html waveform-fsk.png width=800px
 *
 * @sa waveform_psk(), waveform_qam()
 */
inline sptr<WaveForm> waveform_fsk(unsigned int M = 2, float index = 0.4, const ShapingFilterSpec &filter = ShapingFilterSpec::nrz())
{
  return std::make_shared<WaveForm>(nfr::forme_onde_fsk(M, index, filter.fr()));
}



/** @} */






/** @addtogroup telecom-mods
 *  @{
 */




// Attention, démodulateur de protocole ! (à renommer)
template<typename TC, typename TR>
  using ProtocoleDemodulateur = nfr::ProtocoleDemodulateur<TC, TR>;

/** @} */


/** @addtogroup telecom-ps
 *  @{
 */


/** @brief Sample and hold oversampling
 *
 *  <h3>Sample and hold oversampling</h3>
 *
 *  Each sample of the input signal is repeated @f$R@f$ times.
 *  For instance, if @f$R = 2@f$, and @f$x = x_0, x_1, \dots@f$,
 *  then @f$y = x_0, x_0, x_1, x_1, \dots@f$
 *
 *  This function can be used to generate NRZ sequence.
 *
 *  @param x Input vector.
 *  @param R Number of samples to hold (e.g. oversampling factor).
 *  @returns Oversampled signal.
 *
 *
 *  @par Example 1: duplicating values
 *  @code
 *    ArrayXf x(5);
 *    x << 0, 1, 2, 3, 4;
 *    ArrayXf y = sah(x, 2);
 *    // y = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]
 *  @endcode
 *
 *  @par Example 2: generating a pseudo-random NRZ sequence
 *  @code
 *    int nsymbs = 5;  // Number of symbol to generate
 *    int osf    = 10; // Over-sampling ratio
 *    ArrayXf y = sah(randb(nsymbs), osf);
 *    // Will generate nsymbs * osf samples (each symbol is repeated osf times)
 *  @endcode
 */
template<typename T>
Vector<T> sah(const Vector<T> &x, int R)
{
  return nfr::sah(x, R);
}

/** @brief Conversion from binary sequence @f$\to@f$ symbol index sequence.
 *
 * <h3>Conversion from binary sequence @f$\to@f$ symbol index sequence</h3>
 *
 * This function convert a binary vector (0 and 1),
 * to symbols indexes, using @f$k@f$ bits / symbol,
 * according to the following binary encoding:
 * @f[
 * y_i = \sum_{j=0}^{k-1} x_{ki+j} 2^j,\quad i = 0,\dots,(n+k-1)/k
 * @f]
 *
 * @param x Input bit stream (zeros and ones).
 * @param k Number of bits / symbol.
 * @return Vector with @f$(n+k-1)/k@f$ samples.
 *
 * @note If the number of samples @f$n@f$ of the input vector is not a multiple of @f$k@f$, zeros are inserted at the end to complete the last symbol.
 *
 * @sa symunmap_binary()
 *
 */
inline ArrayXi symmap_binary(const BitStream &x, int k)
{
  return nfr::symmap_binaire(x, k);
}

/** @brief Conversion index @f$\to@f$ train binaire
 *
 * <h3>Conversion index @f$\to@f$ train binaire</h3>
 *
 * Cette fonction réalise l'inverse de @ref symmap_binaire(), c'est-à-dire qu'à partir d'une suite de symboles
 * entiers compris entre 0 et @f$2^{k}-1@f$, elle produit une chaine de symboles binaires (0 ou 1).
 *
 * @param x     Input vector.
 * @param k     Number of bits / symbol.
 * @param[out]  bs Binary sequence (zeros and ones).
 *
 *
 * @sa summap_binary()
 */
inline void symunmap_binary(BitStream &bs, const ArrayXi &x, int k)
{
  return nfr::symdemap_binaire(bs, x, k);
}

/** @brief Differential encoder (polynomial = @f$1/(1+X)@f$), MSB first.
 *
 * <h3>Differential encoder</h3>
 *
 * This function build a differentially encoded binary sequence:
 * @f[
 * y_n = x_n \oplus y_{n-1}
 * @f]
 *
 * That is, the transfert function is:
 * @f[
 * P = \frac{1}{1+X}
 * @f]
 *
 * This can be used for instance in DBPSK modulation:
 * in this case, the phase is unchanged for @f$x_n=0@f$,
 * and shifted of 180° for @f$x_n=1@f$.
 *
 * @param x Input binary sequence (@f$n@f$ bits)
 * @param y Output binary sequence (@f$n@f$ bits)
 *
 * @sa diff_decode()
 *
 */
inline void diff_encode(BitStream &y, const BitStream &x)
{
  return nfr::diff_encode(y, x);
}

/** @brief Differential decoder (polynomial = 1+X), MSB first.
 *
 * <h3>Differential decoder</h3>
 *
 * Restauration du signal original à partir d'un signal encodé en différentiel :
 * @f[
 * y_n = x_n \oplus x_{n-1}
 * @f]
 *
 * @note Le bit précédent le premier (@f$x_{-1}@f$) n'est pas connu, aussi il n'est pas possible de calculer
 * @f$y_0@f$. Par conséquent, le train binaire de sortie contiendra 1 bit de moins que le train d'entrée.
 *
 * @param x Input binary sequence (@f$n@f$ bits)
 * @param y Output binary sequence (@f$n-1@f$ bits)
 *
 * @sa diff_encode()
 *
 */
inline void diff_decode(BitStream &y, const BitStream &x)
{
  return nfr::diff_decode(y, x);
}

/** @brief Hard decoding of LLR data.
 *
 *  <h3>Hard decoding of LLR data</h3>
 *
 *  @f[
 *  y_k = 1\textrm{ si }L_k \leq 0, 0\textrm{ sinon.}
 *  @f]
 *
 *  @param llr Vecteur de log-vraisemblances
 *  @param[out] y Train binaire de sortie
 *
 */
inline void decode_hard(BitStream &y, const ArrayXf &llr)
{
  return nfr::decode_hard(y, llr);
}

/** @} */


/** @addtogroup telecom-simu
 *  @{
 */

//
// Parameters
// x: input signal
// σ: square root of noise power
//
// Description
// Compute <latex>$y = x + n$</latex>, with <latex>$n: N(0,\sigma)$</latex>.
// If x is a complex signal, or if complex noise is specifically specified,
// then noise (with same energy) is also added on the imaginary axis.
// So be carefull, with complex noise, the noise power is two times more than for real noise.
//
// Examples
// x = nrz(ts01(10),4);
// y = awgn(x,0.1);
// plot(x,'b'); plot(y,'g');
//
// See also
//  thnoise_power
//  chn_simu
//  fading_chn_init

/** @brief Add (complex) Gaussian white noise.
 *
 * <h3>Simulation of a complex AWGN channel (Additive White Gaussian Noise)</h3>
 *
 * @f[
 * y_k = x_k + b_k^{(r)} + \mathbf{i}\cdot b_k^{(i)}, \quad b^{(r)}, b^{(i)} : \mathcal{N}\left(0,\sigma^2\right)
 * @f]
 *
 * @param x Input vector (complex)
 * @param σ Noise standard deviation (for each I/Q component)
 *
 * @warning The signal being complex, the power of added noise is @f$2\sigma^2@f$.
 *
 * @sa thermal_noise()
 */
inline ArrayXcf awgn_noise(IArrayXcf &x, float σ)
{
  return nfr::bruit_awgn(x, σ);
}

/** @brief Add Gaussian white noise.
 *
 * <h3>Simulation of a real AWGN channel (Additive White Gaussian Noise)</h3>
 *
 * @f[
 * y_k = x_k + b_k , \quad b_k : \mathcal{N}\left(0,\sigma^2\right)
 * @f]
 *
 * @param x Input vector (real)
 * @param σ Noise standard deviation
 * @sa thermal_noise()
 *
 */
inline ArrayXf bruit_awgn(IArrayXf &x, float σ)
{
  return nfr::bruit_awgn(x, σ);
}


/** @brief Dispersive channel type (with or without dominant path) */
enum ChannelType
{
  /** @brief Without dominant path. */
  RAYLEIGH = 0,
  /** @brief With dominant path. */
  RICE
};

/** @brief Configuration for a dispersive channel. */
struct DispersiveChannelConfig : nfr::CanalDispersifConfig
{
  /** @brief Channel type (with or without dominant path) */
  ChannelType &channel_type = *((ChannelType *) &type);

  /** @brief Maximum Doppler frequency. */
  float &max_Doppler = fd;

  /** @brief Sampling frequency. */
  float &fs = fe;

  /** @brief Rician factor. */
  float &Rician_factor = K;
};


/** @brief Creation of a (baseband) dispersive channel simulator.
 *
 * <h3>Dispersive channel simulator (baseband)</h3>
 *
 * This function enable to simulate a Rayleigh (without dominant path)
 * or Rice (with dominant path) channel, in baseband.
 *
 * @param config Configuration (channel type, maximum Doppler and Rician factor).
 * @return %Filtre Baseband signal (cfloat) @f$\to@f$ baseband signal, after fading.
 *
 *
 * @par Example
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_canal_dispersif
 *  @image html canal-dispersif.png width=800px
 */
inline sptr<Filter<cfloat, cfloat, nfr::CanalDispersifConfig>> dispersive_channel(const DispersiveChannelConfig &config)
{
  return nfr::canal_dispersif(config);
}

/** @} */



/** @addtogroup telecom-canalisation
 *  @{
 */

/** @brief Frequency Hopping Spread Sequence configuration */
struct FHSSConfig: nfr::FHSSConfig
{
  /** @brief Séquence de sauts de fréquence */
  Eigen::ArrayXi &hop_sequence = seq;

  /** Facteur de sur-échantillonnage en entrée */
  int &osf_in2 = osf_in;

  /** Facteur de sur-échantillonnage en sortie */
  int &osf_out2 = osf_out;

  /** Fréquence individuelle de chaque saut, normalisée par rapport à la fréquence d'éch d'entrée */
  float &df2 = df;

  /** Durée de chaque slot, en nombre de symboles de sortie */
  int &slot_duration = duree_slot;
};

/** @brief Instanciation of a FHSS (Frequency Hopping Spread Sequence) spreader
 *
 *  <h3>FHSS (Frequency Hopping Spread Sequence)</h3>
 *
 *  @param config  Configuration structure
 *  @return Filtre cfloat @f$\to@f$ cfloat
 */
inline sptr<Filter<cfloat,cfloat,nfr::FHSSConfig>> fhss_modulation(const FHSSConfig &config)
{
  return nfr::fhss_modulation(config);
}


/** @brief DSSS configuration */
struct DSSSConfig: nfr::DSSSConfig
{
  ArrayXf &chips_en = chips;
  // Facteur de sur-échantillonnage en entrée
  int &osf_in_en = osf_in;
};


/** @brief Instanciation of a DSSS (Direct Sequence Spread Sequence) spreader
 *
 *  <h3>DSSS (Direct Sequence Spread Sequence)</h3>
 *
 *  @param config  Configuration structure
 *  @return Filtre cfloat @f$\to@f$ cfloat
 */
inline sptr<Filter<cfloat,cfloat,nfr::DSSSConfig>> dsss_modulation(const DSSSConfig &config)
{
  return nfr::dsss_modulation(config);
}


/** @brief Baseband downconversion configuration. */
struct TranspoBBConfig: nfr::TranspoBBConfig
{
  /** @brief Intermediate frequency (normalized, between 0 and 0.5). */
  float &intermediate_frequency = fi;

  /** @brief Requested downsampling (1 = nonde, 0.5 = 1/2, etc). */
  float &ra_ratio = ratio_ar;
};


/** @brief Baseband downconversion.
 *
 * <h3>Baseband downconversion</h3>
 *
 * Ce bloc permet de convertir un signal radio reçu avec une certaine fréquence intermédiaire, vers un signal
 * bande de base (centré à 0 Hz).
 * Pour cela, les étapes suivantes sont effectuées :
 *  -# Transposition en fréquence du signal entrant via un oscillateur harmonique.
 *  -# Si le signal d'entrée est réel, filtrage du signal image grâce à un filtre RIF.
 *  -# Eventuellement (en fonction de la configuration), décimation du signal bande de base
 *  (en général, une réduction de la fréquence d'échantillonnage est en effet possible, car le signal utile est alors centré à 0 Hz).
 *
 *  La fréquence de coupure du filtre image est réglée ainsi :
 *  @f[
 *    \begin{cases} f & \mbox{si } f < \frac{1}{4}\\
 *    \frac{1}{2} - f & \mbox{sinon} \end{cases}
 *  @f]
 * (dans tous les cas, à mi-chemin entre le signal bande de base et le signal image).
 *
 *  @param config Structure de configuration (fréquence intermédiaire, etc.).
 *
 *  @sa TranspoBBConfig
 */
template<typename T>
  sptr<Filter<T,cfloat,nfr::TranspoBBConfig>> transpo_bb(const TranspoBBConfig &config)
{
  return nfr::transpo_bb<T>(config);
}



/** @} */


/** @addtogroup telecom-crec
 *  @{
 */

using nfr::Ted;

/*struct Ted
{
  // npts : nombre de points nécessaires
  // osf  : Index de sur-échantillonage
  unsigned int npts, osf;
  virtual float calcule(cfloat x0, cfloat x1, cfloat x2) = 0;
};*/


/** @brief Clock recovery configuration structure */
struct ClockRecConfig: nfr::ClockRecConfig
{
  /** @brief Input signal oversampling factor (e.g. ratio of input signal frequency vs symbol frequency) */
  int &osf_en = osf;

  /** @brief Timing error detector object (default is Gardner detector). A ted can be created with the @ref ted_init() function. */
  sptr<Ted> &ted_en = ted;

  /** @brief Interpolator object (default is cardinal cubic spline interpolator).
   *  A interpolator can be created for instance with the @ref itrp_sinc() function.*/
  sptr<tsd::filtrage::Interpolateur<cfloat>> &itrp_en = itrp;

  /** @brief Time constant of the loop, in symbols (default is 5 symbols) */
  float &tc_en = tc;

  /** @brief Enable debug mode (generation of plots) */
  bool &debug_active = debug_actif;

  // Coefficients du filtre adapté
  ArrayXf &matched_filter_coefs = h_fa;
};





/** @brief Creation of a clock recovery object
 *
 * This function will create a clock recovery object, that can be used
 * afterwards to resynchronize and resample an incoming data stream with implicit clock.
 *
 * Optionnaly, the user can choose a specific Timing Error Detector (default is Gardner),
 * and a specific interpolator (default is cardinal cubic spline).
 *
 * Example
 * @code
 * // Create a clock recovery object for an oversampling ratio of 8
 * cr = clock_rec_init(8);
 * @endcode
 *
 *
 * @sa ted_init, itrp_sinc, itrp_cspline, itrp_lineaire **/
inline sptr<FilterGen<cfloat>> clock_rec_new(const ClockRecConfig &config)
{
  return nfr::clock_rec_init(config);
}

inline sptr<FilterGen<cfloat>> clock_rec2_new(const ClockRecConfig &config)
{
  return nfr::clock_rec2_init(config);
}

/** @brief Les différents types de détecteur d'erreur d'horloge */
enum class TedType
{
  GARDNER = 0,
  MM,
  EARLY_LATE
};

/** @brief Création d'un détecteur d'erreur d'horloge (ted / timing error detector) */
inline sptr<Ted> ted_new(TedType type)
{
  return nfr::ted_init((nfr::TedType) type);
}


/** @brief Interface pour un détecteur de phase */
using Ped = nfr::Ped;


/** @brief Les différents types de détecteurs d'erreur de phase */
enum class PedType
{
  AUTO = 0,
  COSTA,
  POWER_LOOP,
  TAN_LOOP,
  DEC_LOOP
};

// M : 2^nb bits / symboles

/** @brief Création d'un détecteur d'erreur de phase (ped / phase error detector) */
inline Ped ped_new(PedType type, sptr<WaveForm> wf)
{
  return nfr::ped_init((nfr::PedType) type, wf->fr);
}


inline Ped ped_costa(int M){return nfr::ped_costa(M);}
inline Ped ped_ploop(int M){return nfr::ped_ploop(M);}
inline Ped ped_tloop(int M){return nfr::ped_tloop(M);}
inline Ped ped_decision(sptr<WaveForm> wf){return nfr::ped_decision(wf->fr);}



/** @brief Interface abstraite pour un filtre de boucle. */
struct LoopFilter
{
  sptr<nfr::FiltreBoucle> fr;

  LoopFilter(sptr<nfr::FiltreBoucle> fr)
  {
    this->fr = fr;
  }

  /** @brief Calcul du prochain déphasage à appliquer à partir de l'erreur de phase courante */
  float step(float err_phase)
  {
    return fr->step(err_phase);
  }

  /** @brief Redémarrage de la boucle */
  void reset()
  {
    fr->reset();
  }
};

/** @brief %Filtre de boucle du premier ordre
 *
 *  <h3>%Filtre de boucle du premier ordre</h3>
 *
 *   Ce filtre consiste tout simplement en un filtre RII du premier ordre :
 *   @f[
 *     \theta_k = \theta_{k-1} + \alpha \cdot e_k
 *   @f]
 *
 *   Le facteur @f$\alpha@f$ étant calculé d'après la constante de temps spécifiée (voir @ref rii1_tc_vers_coef()).
 *
 *  @param τ Constante de temps du filtre (en nombre d'échantillons).
 *
 *  @sa filter_loop_order_2()
 */
inline sptr<LoopFilter> filter_loop_first_order(float τ)
{
  return std::make_shared<LoopFilter>(nfr::filtre_boucle_ordre_1(τ));
}

/** @brief %Filtre de boucle du seconde ordre
 *
 *  <h3>%Filtre de boucle du second ordre</h3>
 *
 *   @f[
 *     \theta_k = \theta_{k-1} + \mu_{k-1}\\
 *     \mu_k    = \mu_{k-1} + \gamma\cdot\left((1+\rho) e_k - e_{k-1}\right)
 *   @f]
 *
 *  Avec :
 *  @f[
 *  \gamma = \frac{16 \eta^2 \cdot B}{A \cdot (1 + 4 \eta^2)}\\
    \rho  = \frac{4 B}{1 + 4 \eta^2}
 *  @f]
 *
 *  @param BL  Bande passante (normalisée à la fréquence d'échantillonnage) de la boucle
 *  @param η   Facteur d'amortissement
 *
 *  @par Bibliography
 *  <i>DVBS2 : Carrier phase synchronization techniques for broadband satellite transmissions, ESA, 2003</i>
 *
 *  @sa filter_loop_order_1()
 */
inline sptr<LoopFilter> filter_loop_second_order(float BL, float η)
{
  return std::make_shared<LoopFilter>(nfr::filtre_boucle_ordre_2(BL, η));
}


/** @} */


/** @addtogroup telecom-mods
 *  @{
 */

/** @brief Paramétrage d'un modulateur numérique */
struct ModConfig
{
  /** @brief Spécifications de la forme d'onde */
  sptr<WaveForm> wf;

  /** @brief Fréquence d'échantillonnage (Hz) */
  float fe = 1;

  /** @brief Fréquence intermédiaire (Hz) */
  float fi = 0;

  /** @brief Fréquence symbole (Hz) */
  float fsymb = 1;

  bool sortie_reelle = true;
  bool debug_actif   = false;

  int ncoefs_filtre_mise_en_forme = 0;

  auto fr() const
  {
    nfr::ModConfig res;
    res.forme_onde = wf->fr;
    res.fe = fe;
    res.fi = fi;
    res.fsymb = fsymb;
    res.sortie_reelle = sortie_reelle;
    res.debug_actif = debug_actif;
    res.ncoefs_filtre_mise_en_forme = ncoefs_filtre_mise_en_forme;
    return res;
  }
};


/** @brief Interface abstraite vers un modulateur */
struct Modulator
{
  sptr<nfr::Modulateur> mod;
  Modulator(sptr<nfr::Modulateur> mod)
  {
    this->mod = mod;
  }


  /** @brief Modulation.
   *
   *  <h3>Modulation</h3>
   *
   *  @param       bs  Train binaire
   *  @return      x   Flot d'échantillons I/Q
   */
  virtual ArrayXcf step(const BitStream &bs)
  {
    return mod->step(bs);
  }


  /** Compléte l'émission avec des échantillons à zéros, filtrés proprement */
  virtual ArrayXcf flush(int nech)
  {
    return mod->flush(nech);
  }

  /** @brief Délais, en nombre d'échantillons.
   *
   * <h3>Délais, en nombre d'échantillons</h3>
   *
   * Nombre d'échantillons entre le premier sorti et le début du premier symbole transmis.
   */
  virtual float delay() const
  {
    return mod->delais();
  }


  // Modifie la forme d'onde
  // Intérêt : si filtre partagé par 2 modulateurs
  // (exemple : modulation différente en-tête et données)
  virtual void set_waveform(sptr<WaveForm> wf)
  {
    mod->def_forme_onde(wf->fr);
  }
};

/** @brief Interface abstraite vers un démodulateur */
struct Demodulator
{
  sptr<nfr::Démodulateur> fr;
  Demodulator(sptr<nfr::Démodulateur> fr)
  {
    this->fr = fr;
  }

  /** @brief Démodulation.
   *
   *  <h3>Démodulation</h3>
   *
   *  @param      x   Flot I/Q à démoduler
   *  @param[out] bs  Train binaire (hard decision)
   */
  virtual void step(const ArrayXcf &x, BitStream &bs)
  {
    fr->step(x, bs);
  }

  /** @brief Démodulation, avec calcul des LLR.
   *
   *  <h3>Démodulation (avec LLR)</h3>
   *
   *  @param      x   Flot I/Q à démoduler
   *  @param[out] bs  Train binaire (hard decision)
   *  @param[out] llr Log-vraisemblances de chaque symbole (une ligne par symbole possible, une colonne par échantillon)
   */
  virtual void step(const ArrayXcf &x, BitStream &bs, ArrayXXf &llr)
  {
    fr->step(x, bs, llr);
  }


  /** @brief Délais, en nombre d'échantillons.
   *
   * <h3>Délais, en nombre d'échantillons</h3>
   *
   * Nombre d'échantillons entre le premier sorti et le début du premier symbole transmis.
   */
  virtual float delay()
  {
    return fr->delais();
  }

  /** Régle le décalage d'horloge, avec un délais compris entre -1 et 1
   */
  virtual void tune_clock(float delais)
  {
    fr->regle_horloge(delais);
  }

  virtual void reset(int cnt = 0)
  {
    fr->reset(cnt);
  }
};


/** @brief Création d'un modulateur numérique.
 *  <h3>Création d'un modulateur numérique</h3>
 *
 * Un modulateur consiste ici à convertir un train binaire en un signal bande de
 * base (ou déjà transposé à une fréquence intermédiaire), mis en forme et sur-échantilloné
 * (de manière à être prêt à être transmis à un ADC).
 *
 * @par Exemple 1 : modulation BPSK (avec filtre NRZ)
 * @snippet exemples/src/sdr/ex-sdr.cc ex_modulateur
 * @image html ex-modulateur.png "Exemples de modulation BPSK" width=800px
 *
 * @par Exemple 2 : modulation QPSK (avec filtre SRRC)
 * @snippet exemples/src/sdr/ex-sdr.cc ex_modulateur2
 * @image html ex-modulateur2.png "Exemples de modulation QPSK" width=800px
 *
 * @sa démodulateur_création()
 */
inline sptr<Modulator> modulator_new(const ModConfig &config)
{
  return std::make_shared<Modulator>(modulateur_création(config.fr()));
}




enum class ItrpType
{
  CSPLINE = 0,
  LINEAIRE,
  LAGRANGE
};




/** @brief Paramétrage d'un démodulateur numérique */
struct DemodConfig
{
  /** @brief Choix de l'architecture du démodulateur (basé ou non sur la décision symbole). */
  enum
  {
    /** @brief Architecture avec détecteurs d'erreurs (horloge et phase) basés sur la décision symbole. */
    ARCHI_AVEC_DECISION = 0,
    /** @brief Architecture avec boucles de correction d'horloge et de phase indépendantes. */
    ARCHI_SANS_DECISION
  } architecture = ARCHI_AVEC_DECISION;

  /** @brief Paramètres utilisés uniquement pour un démodulateur avec architecture <b>avec décision</b>. */
  struct DemodDecConfig
  {
    /** @brief Paramètrage du recouvrement d'horloge. */
    struct
    {
      /** @brief Activation ou non du recouvrement d'horloge */
      bool actif = true;

      /** @brief Constante de temps du filtre de boucle, en nombre de symboles */
      float tc = 100;
    } clock_rec;

    /** @brief Paramètrage du recouvrement de porteuse. */
    struct
    {
      /** @brief Activation ou non du recouvrement de porteuse */
      bool actif = true;

      /** @brief Bande-passante de la boucle (normalisée à la fréquence symbole). */
      float BL = 0.01;

      /** @brief Facteur d'ammortissement */
      float η = 1;
    } carrier_rec;
  } dec;

  /** @brief Paramètres utilisés uniquement pour un démodulateur avec architecture <b>sans décision</b>. */
  struct DemodNDecConfig
  {
    /** @brief Paramètrage calcul RSSI */
    float tc_rssi_coarse = 10;
    float tc_rssi_fine   = 3;

    /** @brief Paramètrage du recouvrement d'horloge */
    struct
    {
      bool actif = true;
      bool mode_ml = false;
      TedType ted = TedType::GARDNER;
      float tc = 5.0f;
      ItrpType itrp = ItrpType::CSPLINE;
      unsigned int itrp_lagrange_degre = 3;
    } clock_rec;

    /** @brief Paramètrage du recouvrement de porteuse */
    struct
    {
      bool actif = true;
      PedType ped = PedType::AUTO;
      float BL = 0.01, η = 1;
    } carrier_rec;
  } ndec;

  /** @brief Affichage des signaux intermédiaires */
  bool debug_actif = false;

  auto fr() const
  {
    nfr::DemodConfig res;
    memcpy((void *) &res, this, sizeof(*this));
    return res;
  }

};


/** @brief Création d'un démodulateur numérique.
 *  <h3>Création d'un démodulateur numérique</h3>
 *
 * Un démodulateur consiste ici à convertir un signal en bande de
 * base (ou transposé à une fréquence intermédiaire) vers un train binaire (ou un train de LLR symboles,
 * pour le cas où un code correcteur est utilisé en aval).
 *
 * La première étape (optionnelle) de la démodulation consiste à transposer un signal centré
 * autour d'une fréquence intermédiaire donnée vers un signal bande de base.
 * Cette étape comprends éventuellement une décimation, afin d'alléger les traitemens ultérieurs
 * (une fois en bande de base, la fréquence d'échanitillonnage peut être réduite) :
 *
 * @image html demod-glob.png "Démodulation - transposition en bande de base" width=400px
 *
 * La deuxième étape est la démodulation des signaux, qui est possible suivant deux architectures décrites ci-après :
 *  - Architecture dite <b>basée sur la décision</b> : c'est l'architecture permettant la meilleure sensibilité,
 *    et c'est donc celle-ci qui est recommandée.
 *
 *  @warning Un inconvénient de cette architecture est qu'il est nécessaire d'avoir un en-tête en début de trame
 *    afin de pré-acrocher les boucles de correction.
 *
 *  - Architecture <b>non basée sur la décision</b>. Elle présente l'avantage de pouvoir s'accrocher sur un signal en cours de route, même si les erreurs initiales de
 *    fréquence ou d'horloge sont importantes, alors que la première nécessite une "pré-accroche".
 *
 *  @warning Contrairement à la première architecture, celle-ci fonctionnera très mal pour des modulations d'ordre élevé (8PSK, QAM, etc.),
 *  du fait des détecteurs d'erreur d'horloge et de phase utilisés, très sensibles à la modulation.
 *
 * <h4>1. Architecture basée sur la décision</h4>
 *
 *  Cette architecture est dite basée sur la décision, car les détections d'erreurs d'horloge et de phase
 *  sont déduite après la "décision" sur chaque symbole, c'est-à-dire le démapping.
 *  Elle est constituée des blocs suivants (dans l'ordre) :
 *    -# Filtrage adapté,
 *    -# Correction de phase / fréquence,
 *    -# Correction d'horloge (interpolation),
 *    -# Démapping (décision symbole le plus proche),
 *    -# Calcul des erreurs d'horloge et de phase / mise à jour des corrections.
 *  @image html demod-archi1.png "Démodulation bande de base - architecture basée sur la décision" width=400px
 *
 *
 * <h5>Correction d'horloge</h5>
 * La correction d'horloge est basée sur le détecteur d'erreur suivant :
 * @f[
 *    \epsilon_k = \frac{\left<d_k - d_{k-1},\ y_{k-1/2} - \frac{d_k + d_{k-1}}{2}\right>}{\left|d_{k-1} - d_k\right|}
 * @f]
 *
 * où @f$d_k,d_{k-1}@f$ sont les derniers symboles décodés (après décision),
 * et @f$y_{1/2}@f$ est l'avant dernier symbole interpolé (après correction d'horloge, mais avant décision),
 * sachant que l'interpolation se fait à deux fois le rythme symbole (autrement dit, @f$y_{k-1/2}@f$ est la valeur reçue à mi-chemin entre les symboles @f$d_k@f$ et @f$d_{k-1}@f$).
 *
 * Ce détecteur, adapté du détecteur de Gardner, présente d'avantage de pouvoir fonctionner pour la plupart des types de modulation.
 *
 * Le filtre de boucle est un simple filtre du premier ordre, dont la constante de temps est réglable (voir @ref DemodConfig).
 *
 * <h5>Correction de phase</h5>
 *  Le détecteur d'erreur de phase est ici très simple, puisque l'on dispose des symboles après décodage :
 *  @f[
 *  \epsilon_k = \widehat{d_k - y_k}
 *  @f]
 *
 *  Le filtre de boucle est un filtre du second ordre, dont les paramètres sont réglables (voir @ref DemodConfig et @ref filtre_boucle_ordre_2()).
 *
 *
 *  <h4>2. Architecture non basée sur la décision</h4>
 *
 *  Cette architecture est constituée des blocs suivants :
 *    -# Filtrage adapté,
 *    -# Boucle de recouvrement d'horloge,
 *    -# Boucle de recouvrement de porteuse,
 *    -# Correction automatique de gain,
 *    -# Démapping des symboles.
 *
 * Dans cette architecture, contrairement à la première, chaque bloc fonctionne indépendemment des autres (il sont juste concaténés en série).
 *
 *
 *
 * @param modconfig Structure de configuration, permettant de choisir la forme d'onde, les fréquences, ... (voir @ref ModConfig)
 * @param demodconfig Structure de configuration, permettant de choisir l'architecture et les paramètres du démodulateur  (voir @ref DemodConfig)
 *
 * @par Exemple : démodulation QPSK
 * @snippet exemples/src/sdr/ex-sdr.cc ex_demodulateur
 * @image html ex-demodulateur.png "Exemples de démodulation QPSK" width=800px
 * Notez que le train binaire démodulé est décalé dans le temps, ceci est du aux filtres utilisés en réception (ainsi qu'à l'interpolation utilisée pour le recouvrement d'horloge).
 *
 *
 *
 * @sa modulateur_création()
 */
inline sptr<Demodulator> demodulator_new(const ModConfig &modconfig,
    const DemodConfig &demodconfig = DemodConfig())
{
  return std::make_shared<Demodulator>(nfr::démodulateur_création(modconfig.fr(), demodconfig.fr()));
}


/** @brief Définition du format d'une trame */
struct TrameFormat
{
  /** Paramètres de modulation */
  ModConfig modulation;

  /** @brief Définition du motif de synchronisation. */
  BitStream entete;

  /** @brief Optionnel : forme d'onde spécifique pour l'en-tête */
  sptr<WaveForm> fo_entete;

  /** @brief Dimension des trames (nombre de bits utiles, après le motif de synchronisation). */
  int nbits = 0;

  auto fr() const
  {
    nfr::TrameFormat res;
    res.modulation = modulation.fr();
    res.entete     = entete;
    res.fo_entete  = fo_entete->fr;
    res.nbits      = nbits;
    return res;
  }
};



/** @brief Structure de configuration d'un récepteur générique */
struct PacketReceiverConfig
{
  // TODO : gestion de plusieurs en-têtes...

  /** @brief Format des trames */
  TrameFormat format;

  /** @brief Configuration du démodulateur. */
  DemodConfig config_demod;

  /** @brief Dimension des blocs d'entrée (?). */
  int BS = 0;

  /** @brief Seuil de détection pour l'en-tête (entre 0 et 1), voir @ref détecteur_création(). */
  float seuil = 0.7;

  /** @brief Minimal SNR (in dB) detected on the header to decode the frame. */
  float SNR_mini = 0;

  /** @brief Si vrai, calcul de la corrélation avec les motifs via des FFT / OLA (sinon filtres RIF). */
  bool correl_fft = false;

  /** @brief Callback optionnelle appelée avec le signal de corrélation normalisée (peut servir pour de la mise au point). */
  std::function<void (const ArrayXf &c)> callback_corr;

  /** @brief Activation ou non des plots de mise au point */
  bool debug_actif = false;

  /** @brief Nombre de coefficient du filtre d'interpolation RIF utilisé avant le démodulateur pour corriger l'horloge. */
  int ncoefs_interpolateur = 15;

  auto fr() const
  {
    nfr::RécepteurConfig res;

    res.format = format.fr();
    res.config_demod = config_demod.fr();
    res.BS  = BS;
    res.seuil = seuil;
    res.SNR_mini = SNR_mini;
    res.correl_fft = correl_fft;
    res.callback_corr = callback_corr;
    res.debug_actif = debug_actif;
    res.ncoefs_interpolateur = ncoefs_interpolateur;

    return res;
  }
};

/** @brief Trame décodée par un récepteur */
struct RécepteurTrame
{
  RécepteurTrame(const nfr::RécepteurTrame &fr)
  {
    det = fr.det;
    bs  = fr.bs;
    EbN0 = fr.EbN0;
    x = fr.x;
    x1 = fr.x1;
  }

  /** @brief Paramètres RF calculés à partir du motif de synchronisation */
  tsd::fourier::Detection det;

  /** @brief Données démodulées */
  BitStream bs;

  /** @brief Rapport signal à bruit normalisé (SNR / bit) */
  float EbN0;

  /** @brief Raw data, before demodulation */
  ArrayXcf x;

  /** @brief Raw data, with clock & phase corrected, before demodulation */
  ArrayXcf x1;
};

struct RécepteurEtat
{
  RécepteurEtat(){}
  RécepteurEtat(const nfr::RécepteurEtat &fr)
  {
    memcpy((void *) this, &fr, sizeof(*this));
  }

  // Dim des blocs d'entrée
  int Ne;
};

/** @brief Interface abstraite vers un récepteur de trames.
 *
 *  Documentation détaillée : @ref récepteur_création() */
struct PacketReceiver
{
  sptr<nfr::Récepteur> fr;

  PacketReceiver(sptr<nfr::Récepteur> fr)
  {
    this->fr = fr;
  }

  /** @brief Fonction de configuration */
  virtual int configure(const PacketReceiverConfig &config)
  {
    return fr->configure(config.fr());
  }

  /** @brief Traitement d'un buffer de données. */
  virtual std::vector<RécepteurTrame> step(const ArrayXcf &x)
  {
    std::vector<RécepteurTrame> res;
    auto r = fr->step(x);
    for(auto &t: r)
      res.push_back(t);
    return res;
  }

  /** @brief Lecture des moniteurs CPU. */
  virtual MoniteursStats moniteurs()
  {
    return fr->moniteurs();
  }

  virtual RécepteurEtat get_etat()
  {
    return fr->get_etat();
  }
};

/** @brief Création d'un récepteur de trame.
 *
 * <h3>Création d'un récepteur de trames</h3>
 *
 * Un récepteur est constitué de trois sous-blocs suivant :
 *  -# Une détecteur d'en-tête de synchronisation.
 *  -# Un interpolateur, permettant de passer au rythme symbole et dont le retard est réglé d'après celui mesuré par le détecteur.
 *  -# Un démodulateur fonctionnant à la fréquence symbole (pas de de correction d'horloge),
 *     dont la boucle de correction de phase (ainsi que la CAG) est
 *     initialisée en fonction des paramètres RF détectés sur le motif de synchronisation.
 *
 * <h4>Description détaillée du fonctionnement</h4>
 *  L'implémentation (src/telecom/recepteur.cc) réalise les opérations suivantes :
 *  1. Les données reçues sont d'abord découpées en bloc de @f$N_e@f$ échantillons
 *     (méthode "step()" appelée depuis l'extérieur), traités ensuite
 *     par la méthode "step_Ne()".
 *  2. Le détecteur de motif (voir @ref détecteur_création()) est appelé sur chaque bloc
 *  3. Pour chaque motif de synchronisation trouvé par le détecteur, les échantillons I/Q correspondant au début des données utiles sont extraits, puis dans la méthode step_demod() :
 *    1. Appel de l'interpolateur (correction d'horloge)
 *    2. Appel du filtre adapté
 *    3. Décimation au rythme symbole
 *    4. Calage du premier échantillon sur le milieu du premier symbole
 *    5. Appel du démodulateur (voir @ref démodulateur_création())
 *
 *  @note Dans une future version, les étapes 3a, 3b et 3c seront fusionnées en un seul filtre d'interpolation polyphase, basé sur le filtre de mise en forme.
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_récepteur
 * @image html ex-recepteur.png width=1000px
 *
 *
 *  @sa émetteur_création(), détecteur_création()
 */
inline sptr<PacketReceiver> packet_receiver_new(const PacketReceiverConfig &rc)
{
  return std::make_shared<PacketReceiver>(nfr::récepteur_création(rc.fr()));
}



/** @brief Structure de configuration d'un récepteur_création générique */
struct PacketEmitterConfig
{
  /** @brief Format des trames */
  TrameFormat format;

  /** @brief Activation ou non des plots de mise au point */
  bool debug_actif = false;

  auto fr() const
  {
    nfr::ÉmetteurConfig res;
    res.format      = format.fr();
    res.debug_actif = debug_actif;
    return res;
  }
};


/** @brief Interface abstraite vers un générateur de trames. */
struct PacketEmitter
{
  sptr<nfr::Émetteur> fr;

  PacketEmitter(sptr<nfr::Émetteur> fr)
  {
    this->fr = fr;
  }

  /** @brief Fonction de configuration */
  virtual int configure(const PacketEmitterConfig &config)
  {
    return fr->configure(config.fr());
  }

  /** @brief Traitement d'un buffer de données. */
  virtual ArrayXcf step(const BitStream &x)
  {
    return fr->step(x);
  }

  /** @brief Lecture des moniteurs CPU. */
  virtual MoniteursStats moniteurs()
  {
    return fr->moniteurs();
  }

  /** @brief Retard, en nombre d'échantillons */
  virtual float delay() const
  {
    return fr->retard();
  }
};

/** @brief Création d'un générateur de trames.
 *
 * <h3>Création d'un générateur de trames.</h3>
 *
 * La structure de configuration (@ref ÉmetteurConfig) indique le format de la trame, c'est-à-dire :
 *   - L'en-tête de synchronisation,
 *   - Le nombre de bits utiles,
 *   - Les paramètres de modulation
 *
 * Ce bloc va concaténer l'en-tête avec les bits utiles de manière à générer des échantillons I/Q à partir des bits utiles,
 * en s'occupant des éventuels problème de padding (si plusieurs bits / symbole).
 * Par exemple, avec une modulation QPSK (2 bits / symboles), si l'en-tête fait 127 bits,
 * alors un zéro est inséré de manière à former un en-tête de 128 bits.
 *
 * En fin de trame, le filtre de mise en forme est appliqué un peu plus loin que nécessaire, afin que le signal I/Q revienne proprement à zéro sans discontinuité (voir exemple ci-dessous, dernière courbe).
 *
 * Notez que la forme d'onde n'est pas forcément identique pour l'en-tête et pour les données utiles,
 * si le champs ÉmetteurConfig::format.fo_entete est renseigné.
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_émetteur
 * @image html ex-emetteur.png width=1000px
 *
 *  @sa récepteur_création(), modulateur_création(), démodulateur_création()
 */
inline sptr<PacketEmitter> packet_emitter_new(const PacketEmitterConfig &ec)
{
  return std::make_shared<PacketEmitter>(nfr::émetteur_création(ec.fr()));
}




/** @} */


/** @addtogroup telecom-simu
 *  @{
 */

//  Modèle statistique Gans, d'après le modèle de Clarke. ?
//  A partir de <i>"A MATLAB-based Object-Oriented Approach to Multipath Fading Channel Simulation"</i>, équation 10

/** @brief Densité spectrale due au Doppler (modèle statistique).
 *
 * <h3>Densité spectrale due au Doppler</h3>
 *
 * Pour un canal avec multiples trajets et un objet mobile, les fréquences reçues pour chaque trajet possible sont décalées
 * suivant le Doppler @f$f_d^{max} \cos \theta@f$, @f$\theta@f$ étant l'angle d'incidence du trajet,
 * et @f$f_d^{max}@f$ le Doppler maximal correspondant à un trajet aligné avec la direction de déplacement
 * du mobile.
 *
 * Si on suppose toutes les directions incidentes équiprobables, sans qu'il n'existe de trajet dominant
 * (modèle de Rayleigh), alors on peut montrer que le décalage Doppler est dristribué suivant la loi (modèle de JAKES) :
 * @f[
 *    P(f) =  \frac{1}{\pi f_d \sqrt{1-\left(\frac{f-fc}{fd}\right)^2}}
 * @f]
 *
 * @param f  Tableau de fréquence où calculer le spectre
 * @param fd Doppler maximum (d'après la vitesse relative max. du mobile)
 * @param fc Fréquence porteuse
 * @return @f$P(f)@f$ (en densité de probalité / Hz)
 *
 * @warning En passant en paramètre le tableau des fréquences (paramètre f), soyez
 * vigilant qu'en général le doppler est une petite valeur par rapport à la fréquence porteuse.
 * De ce fait, il faut absolument représenter ces fréquences
 * en double précision (comme dans l'exemple ci-dessous).
 *
 *  @par Exemple
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_doppler_psd
 *  @image html doppler_psd.png width=600px
 *
 */
inline ArrayXf doppler_distri(ArrayXd f, float fd, double fc)
{
  return nfr::doppler_distri(f, fd, fc);
}


// @brief Compute thermal noise power
// bw: Noise bandwidth, in Hz
/** @brief Calcul de la puissance du bruit thermique.
 *
 * <h3>Puissance du bruit thermique</h3>
 *
 * D'après https://fr.wikipedia.org/wiki/Bruit_thermique :
 * @f[
 * P = k_B \cdot T_k \cdot \Delta_f
 * @f]
 * @f$T_k@f$ étant la température absolue et @f$\Delta_f@f$ la bande passante.
 *
 * @param bp Bande passante, en Hz.
 * @param T  Température ambiante, en dégrés Celsius.
 * @return Puissance du bruit, en Watt
 *
 * @note Pour avoir la puissance du bruit en dBm (1 dBm = 1mW) : @f$P_{dB} = 10 \log_{10}(P \cdot 1000)@f$.
 *
 * @sa bruit_awgn()
 */
inline float thermal_noise(float bp, float T = 25)
{
  return nfr::bruit_thermique(bp, T);
}

/** @brief Paramétrage d'un émulateur de canal de propagation */
struct CPEConfig
{
  auto fr() const
  {
    nfr::ECPConfig res;
    memcpy((void *)&res, this, sizeof(*this));
    return res;
  }


  /** @brief Normalized signal to noise ratio (in dB) */
  float Eb_N0 = 0;

  /** @brief Constant phase offset of the carrier (Hz) */
  float décalage_phase = 0;

  /** @brief Constant frequency offset of the carrier (Hz)  */
  float décalage_fréquence = 0;

  /** @brief Level of phase noise on the carrier (dB/Hz) */
  float phase_noise = 0;

  /** @brief Clock delay, in samples */
  float délais_horloge = 0;

  /** @brief Sampling frequency, Hz */
  float fe = 1;

  /** @brief Symbol frequency, Hz (used to compute SNR) */
  float fsymb = 1;

  /** @brief Data rate, bit/s (used to compute SNR from Eb/N0) */
  float fbit  = 1;

  /** @brief Affichage de courbes intermédiaires */
  bool debug_actif = false;
};

/** @brief Création d'un émulateur de canal de propagation
 *
 *  <h3>Emulateur de Canal de Propagation (ECP)</h3>
 *
 *  @param config Paramétrage
 *  @return Un Filtre cfloat @f$\to@f$ cfloat
 */
inline sptr<Filter<cfloat, cfloat, nfr::ECPConfig>> ecp_new(const CPEConfig &config)
{
  return nfr::ecp_création(config.fr());
}




/** @} */


/** @addtogroup telecom-eq
 *  @{
 */

/** @brief Création d'un égaliseur
 *
 * <h3>Création d'un égaliseur</h3>
 *
 *
 *
 *
 * @sa égaliseur_zfe()
 */
inline sptr<FilterGen<cfloat>> equalizer_fir_new(sptr<WaveForm> wf, const std::string &structure, const std::string &errf,
    float osf, float gain, int N1, int N2)
{
  return nfr::égaliseur_rif_création(wf->fr, structure, errf, osf, gain, N1, N2);
}

/** @brief Compute an inverse filter by zero-forcing.
 *
 * <h3>Compute an inverse filter by zero-forcing</h3>
 *
 * Given the channel impulse response @f$h@f$, this function computes the coefficients
 * of a FIR filter @f$g@f$, trying to approximate:
 * @f[
 * g\star h \sim \delta_d
 * @f],
 *
 * where @f$d@f$ is the global delay of the channel plus equalization filter pair.
 * That is, @f$g@f$ is the inverse filter (up to some delay @f$d@f$) of @f$h@f$.
 *
 * @note Cette fonction requiert de pouvoir mesurer la réponse du canal (par exemple en envoyant un signal de type impulsionnel côté émetteur).
 *
 * @warning
 *  - The inversion is only approximative, as the exact filter inverse would have an infinite number taps.
 *  - Si la réponse du canal présente des zéros (ou des magnitudes faibles) dans le domaine fréquentielle, ce type d'égalisation n'est pas recommandée (amplification du bruit).
 *
 * @param h Réponse impulsionnelle du canal,
 * @param n Nombre de coefficients souhaités pour le filtre inverse.
 * @returns %Filtre RIF inverse (coefficients).
 *
 *  @par Example
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_eg_zfe
 *  @image html zfe-0.png "Impulse responses (channel and equalization filter)" width=800px
 *  <br/>
 *  @image html zfe-1.png "Frequency responses (channel and equalization filter)" width=800px
 *  <br/>
 *  @image html zfe-2.png "Equalization example on a NRZ signal" width=800px
 *
 *
 * @sa égaliseur_création()
 */
inline ArrayXf equalizer_zfe(IArrayXf h, int n)
{
  return nfr::égaliseur_zfe(h, n);
}


/** @} */


/** @addtogroup telecom-simu
 *  @{
 */

// ~english @brief Channel capacity.

/**
 *  @brief Capacité d'un canal AWGN.
 *
 * <h3>Capacité d'un canal</h3>
 *
 * Computes the ideal AWGN channel capacity, in bits/s:
 * @f[
 *   c = B\cdot \log_2(1+\textrm{SNR})
 * @f]
 *
 * @param snr     Signal to noise power ratio (linear scale)
 * @param B       Channel bandwidth, in Hz (default is 1 Hz)
 * @return        Channel capacity, in bits/s
 *
 *  @par Exemple
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_capa
 *  @image html capa.png width=800px
 *
 */
inline float awgn_channel_capacity(float snr, float B = 1)
{
  return nfr::capacite_canal_awgn(snr, B);
}

/** @} */





/** @addtogroup telecom-crec
 *  @{
 */

/** @brief Structure de configuration pour une PLL */
struct PLLConfig
{
  auto fr() const
  {
    nfr::PLLConfig res;

    res.sortie_porteuse = sortie_porteuse;
    res.freq = freq;
    res.loop_filter_order = 2;
    res.bp = bp;
    res.tc = tc;
    res.ped = ped;
    res.debug = debug;

    return res;
  }

  /** Sortie de la PLL :
   *   - si sortie_porteuse = false : le signal résiduel, sinon la porteuse reconstruite
   */
  bool sortie_porteuse = false;

  /** @~french  @brief Fréquence attendue (normalisée)
   *  @~english @brief Expected frequency (normalized) */
  float freq = 0.0;

  int loop_filter_order = 2;

  /** @~french  @brief Bande passante normalisée
   *  @~english @brief Loop bandwidth (normalized) */
  float bp = 0.01;

  /** @brief Time constant, in number of samples
   *  Only used if loop_filter_order = 1
   */
  float tc = 100;

  /** @brief Détecteur d'erreur de phase (optionnel).
   *
   *  Par défaut, le détecteur d'erreur de phase est directement l'argument (y = std::arg(x)).
   *  Cette fonction doit retourner un angle en radians.
   *
   *  @par Exemple 1 (détecteur par défaut) :
   *  @code
   *  config.detecteur_erreur_phase = [](cfloat x)
   *    {
   *      return std::arg(x);
   *    };
   *  @endcode
   *
   *  @par Exemple 2 (boucle quadratique) :
   *  Dans cet exemple, le signal est mis au carré avec de la calculer l'argument,
   *  ce qui permet de détecter la porteuse résiduelle, même en présence d'une modulation BPSK (ou AM).
   *  @code
   *  config.detecteur_erreur_phase = [](cfloat x)
   *    {
   *      return std::arg(x*x) / 2;
   *    };
   *  @endcode
   *
   *   */
  Ped ped;

  /** @~french  @brief Activation du mode de mise au point (tracé des figures)
   *  @~english @brief Activation of the debug mode (plot figures) */
  bool debug = false;
};


/** @brief Structure de configuration pour une PLL à sortie réelle */
struct RPLLConfig
{
  auto fr() const
  {
    nfr::RPLLConfig res;

    res.pll_interne = pll_interne.fr();
    res.freq = freq;
    res.ncoefs_bb = ncoefs_bb;
    res.bp = bp;
    res.debug = debug;

    return res;
  }


  /** @brief Paramètrage commun avec une PLL à sortie complexe */
  PLLConfig pll_interne;

  // TODO : duplicata ?
  float freq;

  /** @brief Activation d'un filtrage autour de la porteuse */
  bool filtre_bb_actif = false;

  /** @brief Nombre de coefficients du filtre */
  int ncoefs_bb = 127;

  /** @~french  @brief Bande passante normalisée
   *  @~english @brief Loop bandwidth (normalized) */
  float bp = 0.01;

  bool debug = false;
};

/** @brief Création d'une PLL (boucle à vérouillage de phase) à sortie réelle
 *
 *  <h3>PLL (réelle)</h3>
 *
 *  @param config Structure de configuration
 *  @returns Filtre générique à entrées / sorties réelles
 *
 *  Cette PLL est capable de se vérouiller sur une sinusoide réelle (ou un signal modulé si
 *  un détecteur d'erreur de phase adéquat est fournit).
 *  Pour un signal complexe (exponentielle complexe), voir la fonction @ref cpll_création().
 *
 *  @~english
 *  @brief Creation of a real-output PLL (phase-locked loop)
 *  @param config Configuration structure
 *  @returns Generic filter with real input / output
 *
 *  This PLL is able to lock on a sinusoidal signal
 *  (or a modulated signal provided an adequat phase error detector is provided).
 *  For a complex carrier (complex exponential), see the function @ref cpll_création().
 */
inline sptr<Filter<float, float, nfr::RPLLConfig>> rpll_new(const RPLLConfig &config)
{
  return nfr::rpll_création(config.fr());
}

/** @brief Création d'une PLL (boucle à vérouillage de phase) à sortie complexe
 *
 *  <h3>PLL (complexe)</h3>
 *
 *  @param config Structure de configuration
 *  @returns Filtre générique à entrées / sorties complexes
 *
 *  Cette PLL (boucle à vérouillage de phase) est capable de se vérouiller
 *  sur une exponentielle complexe
 *  (ou un signal modulé si un détecteur d'erreur de phase adéquat est fournit).
 *  Pour une porteuse réelle (sinusoide), voir la fonction @ref rpll_création().
 *
 *  @~english
 *  @brief Creation of a complex-output PLL (phase-locked loop)
 *  @param config Configuration structure
 *  @returns Generic filter with complex input / output
 *
 *  This PLL is able to lock on a complex carrier
 *  (or a modulated signal provided an adequat phase error detector is provided).
 *  For a real carrier (sinusoidal signal), see the function @ref creation_pll().  */
inline sptr<Filter<cfloat, cfloat, nfr::PLLConfig>> cpll_new(const PLLConfig &config)
{
  return nfr::cpll_création(config.fr());
}

/** @} */


/** @addtogroup telecom-ber
 *  @{
 */


#if 0
/** @brief Résultat de la comparaison de deux chaines binaires */
struct CmpBitsRes
{
  /** @brief Premier vecteur ré-aligné */
  ArrayXf b0;

  /** @brief Deuxième vecteur ré-aligné */
  ArrayXf b1;

  /** @brief Nombre total d'erreurs détectées */
  unsigned int nerr = 0;

  /** @brief Taux d'erreur binaire */
  float ber = 0;

  /** @brief Décalage temporel détecté */
  int decalage = 0;

  /** @brief Déphasage détecté (pour les modulations de type M-PSK) */
  int dec_phase = 0;

  /** @brief Score de corrélation (entre -1 et 1) */
  float score = 0;
};

/** @brief Comparaison de chaines binaires et calcul de taux d'erreur binaire
 *
 * <h3>Comparaison de chaines binaires</h3>
 *
 *  Try to find the best correlation between the 2 bit vectors and
 *  count the number of errors (ignoring the 2 first bits and 2 last bits).
 *
 *  @param b0 Premiere chaine (vecteur de 0 et 1)
 *  @param b1 Deuxième chaine (idem)
 *
 *  @par Exemple
 *  @code
 *  b1 = [0 1 0 0 0 1];
 *  b2 =   [1 0 0 0 1];
 *  auto res = cmp_bits(b1,b2);
 *  // res.ber is the bit error rate
 *  @endcode
 *
 */
inline CmpBitsRes cmp_bits(const BitStream &b0, const BitStream &b1);

/** @brief Idem @ref cmp_bits(), avec gestion des ambiguité de phase M-PSK
 *
 *  <h3>Comparaison de chaines binaires (PSK)</h3>
 *
 */
extern CmpBitsRes cmp_bits_psk(const BitStream &b0, const BitStream &b1, int k);
#endif

/** @} */


/** @addtogroup telecom-plots
 *  @{
 */

// Plot the eye diagram
//
// Parameters
// x: input sequence
// T: symbol period (in samples)
//
// Description
// Plot the eye diagram of a synchronous data signal, which is a scatter plot of the signal
// where the time domain is considered modulo the symbol period (actually using a trigger
// on the signal, to account for symbol period variations).
//
// This diagram is useful to view the impact of ISI (Inter-Symbols Interferences).
// <refsection><title>Example</title></refsection>
// <programlisting>
//T = 128; // Symbol period
//x = nrz(prbs(500),T); // 500 symbols, NRZ shape
//x = ma(x, osf); // moving average
//x = awgn(x, 0.1); // AWGN noise
//clf();
//plot_eye(x, T);
// </programlisting>
// <imageobject><imagedata fileref="ex_eyediagram.png" format="PNG"/></imageobject>
//

/** @brief Diagramme de l'oeil.
 *
 * <h3>Diagramme de l'oeil</h3>
 *
 * Plot the eye diagram of a synchronous data signal, which is a scatter plot of the signal
 * where the time domain is considered modulo the symbol period (actually using a trigger
 * on the signal, to account for symbol period variations).
 *
 * This diagram is useful to view the impact of ISI (Inter-Symbols Interferences).
 *
 * @param f Figure sur laquelle tracer le diagramme.
 * @param x Signal à analyser.
 * @param T Période symbole (en nombre d'échantillons).
 *
 * @par Exemple
 *
 */
inline void plot_eye(tsd::vue::Figure &f, const ArrayXf &x, float T)
{
  return nfr::plot_eye(f, x, T);
}

/** @} */

/** @addtogroup telecom-snr
 *  @{
 */

/** @brief Interface abstraite pour un estimateur de SNR. */
struct SNREstimator
{
  sptr<nfr::EstimateurSNR> fr;

  SNREstimator(sptr<nfr::EstimateurSNR> fr){this->fr = fr;}


  /** @brief Calcul deux vecteurs (S et N) correspondant resp. aux énergies du signal et du bruit, à partir d'un signal bruité x. */
  virtual void step(const ArrayXcf &x, ArrayXf &S, ArrayXf &N)
  {
    return fr->step(x, S, N);
  }
};


//*  Hard-coded coefficients for M-PSK or FSK (constant amplitude) modulation,
//*  and gaussian white noise.
// Matzner algorithm for S and @f$N_0@f$ estimation.

/** @brief Algorithm de Matzner pour l'estimation du niveau de signal et de bruit.
 *
 *  <h3>Algorithm de Matzner pour l'estimation du niveau de signal et de bruit.</h3>
 *
 *
 *  Cette estimateur est basé sur le calcul des moments d'ordres 2 et 4 du signal :
 *  @f[
 *  M_2 = \mathbb{E}[\Vert x\Vert^2], \quad M_4 = \mathbb{E}[\Vert x\Vert^4]
 *  @f]
 *
 *  Notez qu'en pratique ces espérances sont estimées au fil de l'eau via un filtrage exponentiel
 *  de coefficient d'oubli paramétrable.
 *
 *  Alors :
 *
 *  @f[
 *      \hat{S} = \sqrt{2 M_2^2 - M_4},\quad \hat{N} = M_2 - S
 *  @f]
 *
 *  @param γ Coefficient d'oubli du filtre de lissage (pour l'estimation des espérances)
 *
 *  @par Référence :
 *  <i>An SNR estimation algorithm for complex baseband signals using
 *  higher order statistics. R. Matzner, 1993.</i>
 *
 */
inline sptr<SNREstimator> snr_Matzner(float γ = 0.1)
{
  return std::make_shared<SNREstimator>(nfr::snr_Matzner(γ));
}

/** @} */

/** @addtogroup telecom-mods-analog
 *  @{
 */

/** @brief Congiguration modulateur / démodulateur AM */
struct AMConfig
{
  /** @brief Type de modulation AM */
  enum Mode
  {
    /** @brief Double-side band, with carrier */
    DSB = 0,
    /** @brief Double-side band, no carrier */
    DSB_SUPPRESSED_CARRIER,
    /** @brief Single-side band (bande latérale unique), lower side band */
    LSB,
    /** @brief Single-side band (bande latérale unique), upper side band */
    USB
  } mode;

  /** @brief Est-ce une modulation à bande latérale unique ? */
  bool est_BLU() const {return (mode == Mode::LSB) || (mode == Mode::USB);}


  /** @brief Indice de modulation (utilisé seulement en mode DSB) */
  float indice = 1.0;

  /** @brief ? */
  float fe_sortie = 1;

  /** @brief ? */
  float fe_rf = 1;

  /** @brief Fréquence IF ou RF, Hz */
  float f_rf = 1;

  /** @brief Fréquence de coupure du filtre audio passe-bas */
  float fcut_audio_low  = 200;

  /** @brief Fréquence de coupure du filtre audio passe-haut */
  float fcut_audio_high = 8000;

  /** @brief Tracé des signaux intermédiaires */
  bool debug_actif = false;
};


#if 0
/** @brief Congiguration discriminateur FM */
struct FMDiscriConfig
{
  /** @brief Fréquence d'échantillonnage */
  float freq_ech  = 0;
};
#endif

/** @brief Congiguration démodulateur FM */
struct FMDemodConfig
{
  float fe;
  bool genere_img_debug = false;
};

/** @brief Modulation d'amplitude (analogique).
 *
 * <h3>Modulation d'amplitude (analogique)</h3>
 *
 *
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_modulateur_AM
 * @image html ex-modulateur-am.png width=800px
 *
 * @sa demodulateurAM(), modulateurFM()
 *
 */
extern sptr<Filter<float, float, AMConfig>> modulateurAM();

/** @brief TODO */
extern sptr<Filter<cfloat, float, AMConfig>> demodulateurAM();

/** @brief Discrimination polaire pour la démodulation FM.
 *
 * <h3>Discrimination FM</h3>
 *
 * Implémente un discriminateur polaire (calcul de la fréquence instantanée pour un signal en bande de base),
 * qui peut servir de brique pour la démodulation FM :
 * @f[
 * y(t) = \frac{d\arg x}{dt}(t)
 * @f]
 *
 * L'implémentation est basée sur https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/,
 * qui permet un calcul efficace, sans arctangente :
 * @f[
 * y(t)  = \frac{\mathcal{Re}(x(t)) \frac{d\mathcal{Im}(x(t))}{dt} - \mathcal{Im}(x(t)) \frac{d\mathcal{Re}(x(t))}{dt}}{\left|x(t)\right|^2}
 * @f]
 *
 * Les dérivées étant approximées à l'ordre 1 :
 * @f[
 * \frac{dx}{dt}(k) \sim \frac{x_{k+1}-x_{k-1}}{2}
 * @f]
 *
 *
 * @returns   Filtre cfloat (signal bande de base, complexe) vers float (fréquence instantanée, sous forme de pulsation normalisée (entre @f$-\pi@f$ et @f$\pi@f$)).
 * @warning Du fait de l'aproximation utilisée pour le calcul de la dérivée, ce bloc génére un retard de 1 échantillon.
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_discriminateur_fm
 * @image html ex-discriminateur-fm.png width=800px
 */
extern sptr<FilterGen<cfloat,float>> discriminateur_fm();

/** @brief TODO */
extern sptr<Filter<cfloat, cfloat, FMDemodConfig>> demodulateurFM();

/** @} */

/** @addtogroup telecom-codes-synchro
 *  @{
 */



/** @brief Maximal length sequence generation.
 *
 *  <h3>Maximal length sequence generation</h3>
 *
 *  @param n Generating polynomial degree (must be between 1 and 16).
 *
 *  This function generate a binary code of length @f$m=2^n-1@f$, using a shift register
 *  and a primitive polynomial.
 *
 *
 *  @par Example
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_code_mls
 *  @image html ex-code-mls.png width=800px
 *
 *  @sa code_Barker()
 */
extern BitStream code_mls(int n);

/** @brief Génération d'un code de Barker.
 *
 *  <h3>Génération d'un code de Barker</h3>
 *
 *  @param n Longueur du code (2, 3, 4, 5, 7, 11 ou 13)
 *
 *  @par Exemple
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_code_Barker
 *  @image html ex-code-barker.png width=800px
 *
 *  @sa code_mls()
 */
inline BitStream code_Barker(int n)
{
  return nfr::code_Barker(n);
}

/** @cond
 *  Renvoie un polynôme primitif de degré reglen.
 *  Le polynôme renvoyé est stocké "à l'envers", le LSB étant le coefficient
 *  de X^{n-1}, et le MSB celui de X^0. */
inline uint32_t primitive_polynom_binary(int reglen)
{
  return nfr::polynome_primitif_binaire(reglen);
}
/** @endcond */

/** @brief Calcul d'un polynôme primitif.
 *
 *  <h3>Calcul d'un polynôme primitif.</h3>
 *
 *  Cette fonction est utilisée pour la génération de codes à longueur maximale.
 *  Un polynôme primitif @f$p(x)@f$ de degré @f$n@f$ est :
 *   - Irréductible,
 *   - Le plus petit @f$k@f$ tel que @f$p(x)@f$ divise @f$x^k-1@f$ est @f$k=2^n-1@f$.
 *
 *  @param n Degré du polynôme (doit être compris entre 1 et 16).
 *
 *  @sa code_mls()
 */
inline tsd::Poly<int> primitive_polynomial(int n)
{
  return nfr::polynome_primitif(n);
}


/** @} */


/** @brief Vecteur binaire (alternative à la classe BitStream). */
using ArrayHd = ArrayXb;

/** @brief Vecteur de LLR (flottantes) */
using ArrayLLR = ArrayXf;

/** @brief Vecteur de LLR (codées sur 8 bits) */
using ArrayLLRi = Eigen::Array<char, Eigen::Dynamic, 1>;


/** @addtogroup telecom-codes
 *  @{
 */


/** @brief Interface abstraite vers un code correcteur d'erreur. */
struct Code
{
  /** @brief Taille de bloc */
  int n;
  /** @brief Nb bits utiles */
  int k;
  /** @brief Nom du code */
  std::string nom;

  /** @brief Taux de transmission (nb bits utiles / nb bits transmis) */
  float taux() const{return (1.0f * k) / n;}

  /** @brief Fonction d'encodage */
  virtual BitStream encode(const BitStream &u) = 0;

  /** @brief Fonction de décodage */
  virtual BitStream decode(const ArrayLLRi &llri) = 0;
};


/** @} */



}




