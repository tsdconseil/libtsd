#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <vector>
#include <random>



namespace tsd::vue {
class Figures;
}

namespace tsd::telecom {


/** @addtogroup telecom-ps
 *  @{
 */


/** @brief Spécification d'un filtre de mise en forme */
struct SpecFiltreMiseEnForme
{
  /** @brief Pas de filtre, des impulsions brutes sont transmises */
  static SpecFiltreMiseEnForme aucun();

  /** @brief Filtre gaussien + moyenne glissante. */
  static SpecFiltreMiseEnForme gaussien(float BT);

  /** @brief Filtre "NRZ" (moyenne glissante). */
  static SpecFiltreMiseEnForme nrz();

  /** @brief Filtre RCS (racine de cosinus sur-élevé). */
  static SpecFiltreMiseEnForme srrc(float β);

  /** @brief Type de filtre. */
  enum Type
  {
    /** @brief Filtrage "NRZ" (moyenne glissante) */
    NRZ,
    /** @brief Pas de filtre, émisssion d'un train d'impulsions brut. */
    AUCUN,
    /** @brief Filtrage gaussien + moyenne glissante. */
    GAUSSIEN,
    /** @brief Filtre RCS (racine de cosinus sur-élevé). */
    SRRC
  };

  /** @brief Type de filtre. */
  Type type = Type::SRRC;

  /** @brief Produit B * T (filtre Gaussien) */
  float BT            = 0.8;

  /** @brief Facteur de dépassement pour le filtre RCS (voir @ref design_rif_cs()) */
  float β = 0.2;



  /** @brief Calcul des coefficients d'un filtre de mise en forme
   *
   *  <h3>%Filtre de mise en forme - coefficients</h3>
   *
   *  @param ncoefs Nombre de coefficients souhaités
   *  @param osf    Facteur de sur-échantillonnage
   *  @returns      Vecteur des coefficients
   */
  ArrayXf get_coefs(int ncoefs, int osf) const;


  /** @brief Création d'un filtre de mise en forme avec sur-échantillonnage intégré
   *
   *  <h3>%Filtre de mise en forme</h3>
   *
   *  Création d'un filtre de mise en forme avec sur-échantillonnage intégré.
   *  Ce filtre prends accepte donc en entrée directement les symboles à encoder,
   *  et génère un signal sur-échantillonné (facteur R = @f$f_e/f_{symb}@f$ passé en paramètre) et filtré.
   *
   *  Les coefficients du filtre sont normalisés en énergie :
   *  @f[
   *  \sum h_k^2 = R
   *  @f]
   *  de manière à ce que la puissance moyenne en entrée et en sortie du filtre soit identique.
   *
   *  @param ncoefs Nombre de coefficients
   *  @param R      Facteur de sur-échantillonnage
   *
   *  @sa filtre_adapte()
   */
  sptr<FiltreGen<cfloat>> filtre_mise_en_forme(int ncoefs, int R) const;


  /** @brief Idem filtre de mise en forme, mais sans le sur-échantillonnage
   *
   *  <h3>%Filtre adapté</h3>
   *
   *  @param ncoefs Nombre de coefficients
   *  @param osf Facteur de sur-échantillonnage
   *
   *  @sa filtre_mise_en_forme()
   */
  sptr<FiltreGen<cfloat>> filtre_adapte(int ncoefs, int osf) const;

  /** @brief Filtrage adapté et sous-échantillonnage à la fréquence symbole intégré
   *
   *  <h3>%Filtre adapté avec sous-échantillonnage à la fréquence symbole</h3>
   *
   *  @param ncoefs Nombre de coefficients
   *  @param osf Facteur de sur-échantillonnage en entrée
   *
   *  @sa filtre_mise_en_forme(), filtre_adapte()
   */
  sptr<FiltreGen<cfloat>> filtre_adapte_decimation(int ncoefs, int osf) const;

  struct Analyse
  {
    tsd::vue::Figures fig;
  };

  Analyse analyse(int ncoefs, int osf) const;

};

extern std::ostream& operator<<(std::ostream &ss, const SpecFiltreMiseEnForme &t);



/** @}
 *  @addtogroup telecom-mods-wf
 *  @{
 */





/** @brief Spécification d'une forme d'onde */
struct FormeOnde
{
  /** @brief Génération des symboles I/Q à partir d'un flux binaire. */
  virtual ArrayXcf génère_symboles(const BitStream &bs);

  /** @brief Génération des échantillons I/Q à partir d'un flux binaire (y compris filtre de mise en forme).
   *
   * <h3>Génération des échantillons I/Q à partir d'un flux binaire</h3>
   *
   * Cette fonction génère, à partir d'un train binaire, des échantillons I/Q, comprenant le filtre de mise en forme et le sur-échantillonnage.
   *
   * @param bs            Train binaire à encoder.
   * @param ncoefs        Nombre de coefficient pour l'implémentation du filtre.
   * @param osf           Facteur de sur-échantillonnage.
   * @param[out] retard   Retard, en nombre d'échantillon, entre le début du flux de sortie, et le milieu du premier symbole (le retard est du au filtre de mise en forme).
   *
   *
   */
  ArrayXcf génère_échantillons(const BitStream &bs, int ncoefs, int osf, float &retard);


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

  virtual sptr<CtxGen> get_contexte_tx(int ncoefs, int osf);

  /** @brief Décodage des symboles I/Q (par seuillage) et génération d'un train binaire. */
  virtual void decode_symboles(BitStream &bs, const ArrayXcf &x);

  /** @brief Renvoie le ieme symbole de la constellation. */
  virtual cfloat lis_symbole(unsigned int i) const = 0;

  /** @brief Symbole le plus proche parmi les points de la constellation. */
  virtual int symbole_plus_proche(const cfloat &point) const;

  /** @brief Taux d'erreur binaire théorique (pour cette forme d'onde) en fonction du SNR normalisé. */
  virtual float ber(float EbN0_dB) = 0;

  /** @brief Taux d'erreur binaire théorique en fonction du SNR normalisé. */
  ArrayXf ber(const ArrayXf &EbN0_dB);

  /** @brief Renvoie les points de la constellation. */
  virtual ArrayXcf constellation() const = 0;

  /** @brief Excursion fréquentielle, en multiple de la fréquence symbole. */
  virtual float excursion() const;

  /** @brief Renvoie une description de la modulation (courte chaine de caractères). */
  virtual std::string desc() const = 0;

  /** @brief Renvoie une description de la modulation (courte chaine de caractères). */
  virtual std::string desc_courte() const {return desc();}



  struct Infos
  {
    /** @brief Vrai pour les modulations PSK, ASK, QAM. */
    bool est_lineaire = true;

    /** @brief Indique une modulation de phase */
    bool est_psk = false;

    /** @brief Indique une modulation d'amplitude */
    bool est_ask = false;

    /** @brief Indique une modulation de fréquence */
    bool est_fsk = false;

    /** @brief Indique une modulation en quadrature */
    bool est_qam = false;

    /** @brief Indice de modulation (pour les modulations FSK) */
    float index = 1.0f;

    /** @brief Nombre de symboles possibles */
    int M;

    /** @brief Nombre de bits par symbole (@f$\log_2(M)@f$) */
    int k;
  };

  Infos infos;

  int cnt = 0;

  /** @brief Spécification du filtre de mise en forme */
  SpecFiltreMiseEnForme filtre;

};

extern std::ostream& operator<<(std::ostream &ss, const FormeOnde &t);

/** @brief Création d'une forme d'onde de type modulation de phase.
 *
 * <h3>Modulation de phase</h3>
 *
 * Création d'une forme d'onde M-PSK. Le résultat peut être utilisé pour créer un modulateur
 * (@ref modulateur_création()) ou un démodulateur (@ref démodulateur_création()).
 *
 * @param M Nombre de bits / symboles.
 * @param filtre Filtre de mise en forme (par défaut : NRZ)
 * @return Pointeur vers une forme d'onde abstraite (@ref FormeOnde).
 *
 * @warning @f$M@f$ doit être une puissance de 2.
 *
 * @par Exemple 1 : tracé de quelques constellations PSK
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_psk
 * @image html waveform-psk.png "Exemples de formes d'ondes PSK : BPSK, QPSK, 8PSK, 16PSK" width=800px
 *
 * @par Exemple 2 : Calcul de taux d'erreur binaires
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_psk2
 * @image html waveform-psk2.png "Taux d'erreur binaire" width=800px
 *
 * @sa forme_onde_qam(), forme_onde_qpsk()
 */
extern sptr<FormeOnde> forme_onde_psk(unsigned int M, const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde BPSK.
 *
 * <h3>Création d'une forme d'onde BPSK</h3>
 *
 * Cette fonction est un raccourci vers @ref forme_onde_psk() pour M = 2.
 *
 * @sa forme_onde_psk(), forme_onde_qam(), forme_onde_fsk()
 */
extern sptr<FormeOnde> forme_onde_bpsk(const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde M-ASK.
 *
 * <h3>Création d'une forme d'onde M-ASK</h3>
 *
 * @f[
 * x_n = K_1 + \frac{s_n}{M-1}\cdot K_2
 * @f]
 *
 *
 * @sa forme_onde_psk(), forme_onde_bpsk(), forme_onde_qam(), forme_onde_fsk()
 */
extern sptr<FormeOnde> forme_onde_ask(int M = 2, float K1 = -1, float K2 = 2, const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde QPSK.
 *
 * <h3>Création d'une forme d'onde QPSK</h3>
 *
 * Cette fonction est un raccourci vers @ref forme_onde_psk() pour M = 4.
 *
 * @sa forme_onde_psk(), forme_onde_qam(), forme_onde_fsk()
 */
extern sptr<FormeOnde> forme_onde_qpsk(const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde π/4 - QPSK.
 *
 * <h3>Création d'une forme d'onde π/4 - QPSK</h3>
 *
 *
 * @sa forme_onde_qpsk()
 */
extern sptr<FormeOnde> forme_onde_π4_qpsk(const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde QAM
 *
 * <h3>Forme d'onde QAM</h3>
 *
 * Création d'une forme d'onde en modulation d'amplitude en quadrature :
 * les points de la constellation forment une grille régulière.
 *
 * @par Exemple : tracé de constellations QAM
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_qam
 * @image html waveform-qam.png "Constellations QAM16, QAM64, QAM256" width=800px
 *
 * @sa forme_onde_psk(), forme_onde_fsk()
 */
extern sptr<FormeOnde> forme_onde_qam(unsigned int M, const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());

/** @brief Création d'une forme d'onde FSK.
 *
 * <h3>Forme d'onde FSK</h3>
 *
 * Création d'une forme d'onde en modulation de fréquence (FSK, pour <i>Frequency Shift Keying</i>).
 *
 * Cette modulation est caractérisée par <b>l'indice de modulation</b>, qui est le rapport entre l'excursion maximale (deux fois la déviation) et le débit symbole :
 * @f[
 * h = \frac{2 \Delta f}{f_{symb}}
 * @f]
 *
 * la fréquence instantanée variant entre @f$f_c - \Delta f@f$ et @f$f_c + \Delta f@f$.
 *
 * @param M Nombre de valeurs possibles par symbole,
 * @param index Indice de modulation
 * @param filtre  Type de filtrage
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_waveform_fsk
 * @image html waveform-fsk.png width=800px
 *
 * @sa forme_onde_psk(), forme_onde_qam()
 */
extern sptr<FormeOnde> forme_onde_fsk(unsigned int M = 2, float index = 0.4, const SpecFiltreMiseEnForme &filtre = SpecFiltreMiseEnForme::nrz());



/** @} */






/** @addtogroup telecom-mods
 *  @{
 */


// Attention, démodulateur de protocole ! (à renommer)
template<typename TC, typename TR>
struct ProtocoleDemodulateur
{
  virtual ~ProtocoleDemodulateur(){};
  virtual int configure(const TC &config) = 0;
  virtual std::vector<TR> step(const ArrayXcf &x) = 0;
};

/** @} */


/** @addtogroup telecom-ps
 *  @{
 */


/** @brief Sample and hold oversampling
 *
 *  <h3>Sample and hold oversampling</h3>
 *
 *  Each sample of the input signal is repeated @p R times.
 *  For instance, if R = 2, and @f$x = x_0, x_1, \dots@f$, then @f$y = x_0, x_0, x_1, x_1, \dots@f$
 *
 *  This function can be used to generate NRZ data stream.
 *
 *  @param x input vector
 *  @param R number of samples to hold.
 *  @returns Oversampled signal
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
 *  @par Example 2: generating a random NRZ data stream
 *  @code
 *    int nsymbs = 5;  // Number of symbol to generate
 *    int osf    = 10; // Over-sampling ratio
 *    ArrayXf y = sah(randb(nsymbs), osf);
 *    // Will generate nsymbs * osf samples (each symbol is repeated osf times)
 *  @endcode
 */
template<typename T>
Vecteur<T> sah(const Vecteur<T> &x, int R)
{
  int n = x.rows();
  Vecteur<T> y(n * R);
  for(auto i = 0; i < n; i++)
    y.segment(i*R, R).setConstant(x(i));
  return y;
}

/** @brief Conversion train binaire @f$\to@f$ index.
 *
 * <h3>Conversion train binaire @f$\to@f$ index</h3>
 *
 * Cette fonction convertit un vecteur binaire (valeurs : 0 ou 1),
 * en un vecteur de symboles, avec @f$k@f$ bits / symboles,
 * suivant l'encodage binaire standard :
 * @f[
 * y_i = \sum_{j=0}^{k-1} x_{ki+j} 2^j,\quad i = 0,\dots,(n+k-1)/k
 * @f]
 *
 * @param x Signal d'entrée (valeurs : 0 ou 1)
 * @param k Nombre de bits par symbole
 * @return Vecteurs de @f$(n+k-1)/k@f$ échantillons.
 *
 * @note Si le nombre d'échantillons @f$n@f$ du signal d'entrée n'est pas un multiple de @f$k@f$, des zéros sont insérés à la fin de manière à produire le dernier échantillon.
 *
 * @sa symdemap_binaire()
 *
 */
extern ArrayXi symmap_binaire(const BitStream &x, int k);

/** @brief Conversion index @f$\to@f$ train binaire
 *
 * <h3>Conversion index @f$\to@f$ train binaire</h3>
 *
 * Cette fonction réalise l'inverse de @ref symmap_binaire(), c'est-à-dire qu'à partir d'une suite de symboles
 * entiers compris entre 0 et @f$2^{k}-1@f$, elle produit une chaine de symboles binaires (0 ou 1).
 *
 * @param x Signal d'entrée
 * @param k Nombre de bits par symbole
 * @param[out] bs Chaine binaire (valeurs : 0 ou 1)
 *
 */
extern void symdemap_binaire(BitStream &bs, const ArrayXi &x, int k);

/** @brief Differential encoder (polynomial = @f$1/(1+X)@f$), MSB first.
 *
 * <h3>Encodeur différentiel</h3>
 *
 * Cette fonction génére un train binaire encodé de manière différentielle :
 * @f[
 * y_n = x_n \oplus y_{n-1}
 * @f]
 *
 * Soit la fonction de transfert :
 * @f[
 * P = \frac{1}{1+X}
 * @f]
 *
 * Typiquement utilisé en DBPSK : dans ce cas, la phase est inchangée pour @f$x_n=0@f$,
 * et décalée de 180° pour @f$x_n=1@f$.
 *
 * @param x Train binaire d'entrée (@f$n@f$ bits)
 * @param y Train binaire de sortie (@f$n@f$ bits)
 *
 * @sa diff_decode()
 *
 */
extern void diff_encode(BitStream &y, const BitStream &x);

/** @brief Differential decoder (polynomial = 1+X), MSB first.
 *
 * <h3>Décodeur différentiel</h3>
 *
 * Restauration du signal original à partir d'un signal encodé en différentiel :
 * @f[
 * y_n = x_n \oplus x_{n-1}
 * @f]
 *
 * @note Le bit précédent le premier (@f$x_{-1}@f$) n'est pas connu, aussi il n'est pas possible de calculer
 * @f$y_0@f$. Par conséquent, le train binaire de sortie contiendra 1 bit de moins que le train d'entrée.
 *
 * @param x Train binaire d'entrée (@f$n@f$ bits)
 * @param y Train binaire de sortie (@f$n-1@f$ bits)
 *
 * @sa diff_encode()
 *
 */
extern void diff_decode(BitStream &y, const BitStream &x);

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
extern void decode_hard(BitStream &y, const ArrayXf &llr);

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

/** @brief Ajoute un bruit blanc gaussien complexe.
 *
 * <h3>Simulation d'un canal complexe AWGN (Additive White Gaussian Noise)</h3>
 *
 * @f[
 * y_k = x_k + b_k^{(r)} + \mathbf{i}\cdot b_k^{(i)}, \quad b^{(r)}, b^{(i)} : \mathcal{N}\left(0,\sigma^2\right)
 * @f]
 *
 * @param x Signal d'entrée (complexe)
 * @param σ Ecart-type du bruit
 *
 * @warning Le signal étant complexe, la puissance du bruit ajoutée est de @f$2\sigma^2@f$.
 *
 * @sa bruit_thermique()
 */
extern ArrayXcf bruit_awgn(IArrayXcf &x, float σ);

/** @brief Ajoute un bruit blanc gaussien réel.
 *
 * <h3>Simulation d'un canal réel AWGN (Additive White Gaussian Noise)</h3>
 *
 * @f[
 * y_k = x_k + b_k , \quad b_k : \mathcal{N}\left(0,\sigma^2\right)
 * @f]
 *
 * @param x Signal d'entrée (réel)
 * @param σ Ecart-type du bruit
 * @sa bruit_thermique()
 *
 */
extern ArrayXf bruit_awgn(IArrayXf &x, float σ);


/** @brief Type de canal (avec ou sans trajet dominant) */
enum TypeCanal
{
  /** @brief Sans trajet dominant */
  RAYLEIGH = 0,
  /** @brief Avec trajet dominant */
  RICE
};

/** @brief Configuration pour un canal dispersif */
struct CanalDispersifConfig
{
  /** Type de canal (avec ou sans trajet dominant) */
  TypeCanal type;

  /** @brief Fréquence Doppler max */
  float fd;

  /** @brief Fréquence d'échantillonnage */
  float fe;

  /** @brief Facteur Ricien. */
  float K;
};


/** @brief Création d'un simulateur de canal dispersif.
 *
 * <h3>Simulateur de canal dispersif</h3>
 *
 * Cet objet permet de simuler un canal de Rayleigh (sans trajet dominant)
 * ou de Rice (avec trajet dominant), en bande de base.
 *
 * @param config Configuration (type de canal, Doppler max et facteur Ricien).
 * @return %Filtre signal bande de base (cfloat) @f$\to@f$ Signal bande de base, après atténuation.
 *
 *
 * @par Exemple
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_canal_dispersif
 *  @image html canal-dispersif.png width=800px
 */
extern sptr<Filtre<cfloat, cfloat, CanalDispersifConfig>> canal_dispersif(const CanalDispersifConfig &config);

/** @} */



/** @addtogroup telecom-canalisation
 *  @{
 */

/** @brief Frequency Hopping Spread Sequence configuration */
struct FHSSConfig
{
  /** @brief Séquence de sauts de fréquence */
  Eigen::ArrayXi seq;
  /** Facteur de sur-échantillonnage en entrée */
  int osf_in;
  /** Facteur de sur-échantillonnage en sortie */
  int osf_out;
  /** Fréquence individuelle de chaque saut, normalisée par rapport à la fréquence d'éch d'entrée */
  float df;
  /** Durée de chaque slot, en nombre de symboles de sortie */
  int duree_slot = 0;
};

/** @brief Instanciation of a FHSS (Frequency Hopping Spread Sequence) spreader
 *
 *  <h3>FHSS (Frequency Hopping Spread Sequence)</h3>
 *
 *  @param config  Configuration structure
 *  @return Filtre cfloat @f$\to@f$ cfloat
 */
extern sptr<Filtre<cfloat,cfloat,FHSSConfig>> fhss_modulation(const FHSSConfig &config);


/** @brief DSSS configuration */
struct DSSSConfig
{
  ArrayXf chips;
  // Facteur de sur-échantillonnage en entrée
  int osf_in;
};


/** @brief Instanciation of a DSSS (Direct Sequence Spread Sequence) spreader
 *
 *  <h3>DSSS (Direct Sequence Spread Sequence)</h3>
 *
 *  @param config  Configuration structure
 *  @return Filtre cfloat @f$\to@f$ cfloat
 */

extern sptr<Filtre<cfloat,cfloat,DSSSConfig>> dsss_modulation(const DSSSConfig &config);


/** @brief Configuration d'une transposition en bande de base */
struct TranspoBBConfig
{
  /** @brief Fréquence intermédiaire (normalisée, entre 0 et 0,5). */
  float fi = 0;

  /** @brief Adaptation de rythme demandée (1 = aucune, 0.5 = 1/2, etc) */
  float ratio_ar = 1;
};


/** @brief Transposition de fréquence à partir d'un signal réel ou complexe.
 *
 * <h3>Transposition de fréquence à partir d'un signal réel ou complexe</h3>
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
  sptr<Filtre<T,cfloat,TranspoBBConfig>> transpo_bb(const TranspoBBConfig &config);



/** @} */


/** @addtogroup telecom-crec
 *  @{
 */

struct Ted
{
  // npts : nombre de points nécessaires
  // osf  : Index de sur-échantillonage
  unsigned int npts, osf;
  virtual float calcule(cfloat x0, cfloat x1, cfloat x2) = 0;//const ArrayXcf &x) = 0;
};


/** @brief Clock recovery configuration structure */
struct ClockRecConfig
{
  /** @brief Input signal oversampling factor (e.g. ratio of input signal frequency vs symbol frequency) */
  int osf = 8;

  /** @brief Timing error detector object (default is Gardner detector). A ted can be created with the @ref ted_init() function. */
  sptr<Ted> ted;

  /** @brief Interpolator object (default is cardinal cubic spline interpolator).
   *  A interpolator can be created for instance with the @ref itrp_sinc() function.*/
  sptr<tsd::filtrage::Interpolateur<cfloat>> itrp;

  /** @brief Time constant of the loop, in symbols (default is 5 symbols) */
  float tc = 5;

  /** @brief Enable debug mode (generation of plots) */
  bool debug_actif = false;

  // Coefficients du filtre adapté
  ArrayXf h_fa;
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
extern sptr<FiltreGen<cfloat>> clock_rec_init(const ClockRecConfig &config);

extern sptr<FiltreGen<cfloat>> clock_rec2_init(const ClockRecConfig &config);

/** @brief Les différents types de détecteur d'erreur d'horloge */
enum class TedType
{
  GARDNER = 0,
  MM,
  EARLY_LATE
};

/** @brief Création d'un détecteur d'erreur d'horloge (ted / timing error detector) */
extern sptr<Ted> ted_init(TedType type);


/** @brief Interface pour un détecteur de phase */
using Ped = std::function<float (cfloat x)>;

/*struct Ped
{
  virtual float calcule(const std::complex<float> &x) = 0;
  std::string nom;
  unsigned int M = 2;
  bool require_agc = false;
  float agc_tc = 3.0f;
};*/

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
extern Ped ped_init(PedType type, sptr<FormeOnde> wf);


extern Ped ped_costa(int M);
extern Ped ped_ploop(int M);
extern Ped ped_tloop(int M);
extern Ped ped_decision(sptr<FormeOnde> wf);

/** @brief Interface abstraite pour un filtre de boucle. */
struct FiltreBoucle
{
  /** @brief Calcul du prochain déphasage à appliquer à partir de l'erreur de phase courante */
  virtual float step(float err_phase) = 0;

  /** @brief Redémarrage de la boucle */
  virtual void reset() = 0;
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
 *  @sa filtre_boucle_ordre_2()
 */
extern sptr<FiltreBoucle> filtre_boucle_ordre_1(float τ);

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
 *  @par Bibliographie
 *  <i>DVBS2 : Carrier phase synchronization techniques for broadband satellite transmissions, ESA, 2003</i>
 *
 *  @sa filtre_boucle_ordre_1()
 */
extern sptr<FiltreBoucle> filtre_boucle_ordre_2(float BL, float η);


/** @} */


/** @addtogroup telecom-mods
 *  @{
 */

/** @brief Paramétrage d'un modulateur numérique */
struct ModConfig
{
  /** @brief Spécifications de la forme d'onde */
  sptr<FormeOnde> wf;

  /** @brief Fréquence d'échantillonnage (Hz) */
  float fe = 1;

  /** @brief Fréquence intermédiaire (Hz) */
  float fi = 0;

  /** @brief Fréquence symbole (Hz) */
  float fsymb = 1;

  bool sortie_reelle = true;
  bool debug_actif   = false;


  int ncoefs_filtre_mise_en_forme = 0;
};


/** @brief Interface abstraite vers un modulateur */
struct Modulateur
{
  /** @brief Modulation.
   *
   *  <h3>Modulation</h3>
   *
   *  @param       bs  Train binaire
   *  @return      x   Flot d'échantillons I/Q
   */
  virtual ArrayXcf step(const BitStream &bs) = 0;


  /** Compléte l'émission avec des échantillons à zéros, filtrés proprement */
  virtual ArrayXcf flush(int nech) = 0;

  /** @brief Délais, en nombre d'échantillons.
   *
   * <h3>Délais, en nombre d'échantillons</h3>
   *
   * Nombre d'échantillons entre le premier sorti et le début du premier symbole transmis.
   */
  virtual float delais() const = 0;


  /** @brief Modifie la forme d'onde en cours de route.
   *
   *  <h3>Modifie la forme d'onde en cours de route</h3>
   *
   *
   *  Intérêt : si le filtre de mise en forme est partagé par 2 modulateurs,
   *  ce qui peut arriver par exemple si la modulation est différente
   *  pour l'en-tête et les données (au moment du changement de forme d'onde,
   *  l'état du filtre de mise en forme est préservé).
   */
  virtual void def_forme_onde(sptr<FormeOnde> fo) = 0;
};

/** @brief Interface abstraite vers un démodulateur */
struct Démodulateur
{
  /** @brief Démodulation.
   *
   *  <h3>Démodulation</h3>
   *
   *  @param      x   Flot I/Q à démoduler
   *  @param[out] bs  Train binaire (hard decision)
   */
  virtual void step(const ArrayXcf &x, BitStream &bs){ArrayXXf llr; step(x, bs, llr);}

  /** @brief Démodulation, avec calcul des LLR.
   *
   *  <h3>Démodulation (avec LLR)</h3>
   *
   *  @param      x   Flot I/Q à démoduler
   *  @param[out] bs  Train binaire (hard decision)
   *  @param[out] llr Log-vraisemblances de chaque symbole (une ligne par symbole possible, une colonne par échantillon)
   */
  virtual void step(const ArrayXcf &x, BitStream &bs, ArrayXXf &llr) = 0;


  /** @brief Délais, en nombre d'échantillons.
   *
   * <h3>Délais, en nombre d'échantillons</h3>
   *
   * Nombre d'échantillons entre le premier sorti et le début du premier symbole transmis.
   */
  virtual float delais() = 0;

  /** Régle le décalage d'horloge, avec un délais compris entre -1 et 1
   */
  virtual void regle_horloge(float delais){}

  virtual void reset(int cnt = 0) = 0;
};


/** @brief Création d'un modulateur numérique.
 *  <h3>Création d'un modulateur numérique</h3>
 *
 *
 * Le bloc modulateur permet de convertir un train binaire en un
 * signal bande de base (ou déjà transposé à une fréquence intermédiaire),
 * mis en forme et sur-échantillonné (de manière à être prêt à être transmis à un ADC).
 *
 *
 * La structure de paramètrage (@ref ModConfig) spécifie la forme d'onde
 * (type de modulation, filtre de mise en forme),
 * ainsi que les paramètres de fréquence, notamment :
 * - La fréquence symbole @f$f_{symb}@f$ (<code>fsymb</code>),
 * - La fréquence d'échantillonnage souhaitée en sortie @f$f_e@f$ (<code>fe</code>),
 * - La fréquence intermédiaire souhaitée en sortie @f$f_i@f$ (<code>fi</code>).
 *
 *
 * Le facteur de sur-échantillonnage global vaut donc :
 * @f[
 * R = \frac{f_e}{f_{symb}}
 * @f]
 *
 * Afin de réduire la complexité de calcul,
 * ce facteur global est décomposé en deux
 * parties :
 * @f[
 * R = R_1 \cdot R_2
 * @f]
 *
 * Le facteur @f$R_1@f$ sera appliqué au moment de l'application
 * du filtre de mise en forme, et le facteur @f$R_2@f$ par un interpolateur final.
 *
 * Les différentes étapes sont les suivantes :
 *   - <b>Génération des symboles :</b> les bits d'entrées sont regroupés par
 *   groupes de @f$k@f$ bits (@f$k@f$ étant le nombre de bits par symbole de la forme d'onde configurée),
 *   et chaque groupe de @f$k@f$ bits est transformé en un symbole de la constellation.
 *   - <b>%Filtre de mise en forme et sur-échantillonnage :</b>
 *   le filtre est implémenté sous forme polyphase et permet de passer
 *   de 1 échantillon / symbole à @f$R_1@f$ échantillons / symbole :
 *    @f[
 *     x_n = \sum_k h_k \cdot u_{(n-k)/R_1}
 *    @f]
 *   - <b>Adaptation de rythme :</b>
 *   cet étage d'interpolation permet de passer
 *   de @f$R_1@f$ à @f$R=R_1\cdot R_2@f$ échantillons / symbole.
 *   - <b>Transposition :</b>
 *   cet étage (actif uniquement si le paramètre @f$f_i@f$ est non nul)
 *   permet de transposer le signal bande de base vers
 *   la fréquence intermédiaire @f$f_i@f$ :
 *   @f[
 *     y_k = e^{2\pi\mathbf{i} k f_i / f_e} \cdot x_k
 *    @f]
 *      Par ailleurs, si le paramètre <code>sortie_réelle</code> est actif,
 *      la partie imaginaire du résultat est mise à zéro :
 * @f[
 *  z_k = \mathcal{R}(y_k)
 * @f]
 *
 * @warning Attention :
 *   L'étage d'adaptation de rythme est pour l'instant désactivé, si bien que le changement de rythme
 *   est entièrement à la charge du filtre de mise en forme (implémentation RIF polyphase),
 *   ce qui peut être assez coûteux en charge de calcul si @f$R@f$ (facteur global de rééchantillonnage)
 *   est élevé.
 *
 * @par Schéma-bloc
 * <img src="modulateur.png" align="left" width="300px"/>
 * <div style="clear: both"></div>
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
extern sptr<Modulateur> modulateur_création(const ModConfig &config);


//  * @image html modulateur.png "Architecture modulateur" width=300px

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
extern sptr<Démodulateur> démodulateur_création(const ModConfig &modconfig,
                                                const DemodConfig &demodconfig = DemodConfig());


/** @brief Définition du format d'une trame */
struct TrameFormat
{
  /** Paramètres de modulation */
  ModConfig modulation;

  /** @brief Définition du motif de synchronisation. */
  BitStream entete;

  /** @brief Optionnel : forme d'onde spécifique pour l'en-tête */
  sptr<FormeOnde> fo_entete;

  /** @brief Dimension des trames (nombre de bits utiles, après le motif de synchronisation). */
  int nbits = 0;
};



/** @brief Structure de configuration d'un récepteur générique */
struct RécepteurConfig
{
  // TODO : gestion de plusieurs en-têtes...

  /** @brief Format des trames */
  TrameFormat format;

  // std::vector<TrameFormat> formats;

  // TODO: redondant avec format.modulation
  //float fe = 0, fsymb = 0, fi = 0;

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
};

/** @brief Trame décodée par un récepteur */
struct RécepteurTrame
{
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
  // Dim des blocs d'entrée
  int Ne;
};

/** @brief Interface abstraite vers un récepteur de trames.
 *
 *  Documentation détaillée : @ref récepteur_création() */
struct Récepteur
{
  /** @brief Fonction de configuration */
  virtual int configure(const RécepteurConfig &config) = 0;

  /** @brief Traitement d'un buffer de données. */
  virtual std::vector<RécepteurTrame> step(const ArrayXcf &x) = 0;

  /** @brief Lecture des moniteurs CPU. */
  virtual MoniteursStats moniteurs() = 0;

  virtual RécepteurEtat get_etat() = 0;
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
extern sptr<Récepteur> récepteur_création(const RécepteurConfig &rc);



/** @brief Structure de configuration d'un récepteur_création générique */
struct ÉmetteurConfig
{
  /** @brief Format des trames */
  TrameFormat format;

  /** @brief Activation ou non des plots de mise au point */
  bool debug_actif = false;
};


/** @brief Interface abstraite vers un générateur de trames. */
struct Émetteur
{
  /** @brief Fonction de configuration */
  virtual int configure(const ÉmetteurConfig &config) = 0;

  /** @brief Traitement d'un buffer de données. */
  virtual ArrayXcf step(const BitStream &x) = 0;

  /** @brief Lecture des moniteurs CPU. */
  virtual MoniteursStats moniteurs() = 0;

  /** @brief Retard, en nombre d'échantillons */
  virtual float retard() const = 0;
};

/** @brief Création d'un générateur de trames.
 *
 * <h3>Création d'un générateur de trames.</h3>
 *
 * La structure de configuration (@ref ÉmetteurConfig) indique le format de la trame, c'est-à-dire :
 *   - L'en-tête de synchronisation,
 *   - Le nombre de bits utiles,
 *   - Les paramètres de la modulation.
 *
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
 *
 * @par Schéma-bloc
 * <img src="emetteur.png" align="left" width="500px"/>
 * <div style="clear: both"></div>
 *
 * @warning
 * L'application d'un code correcteur n'est pas encore implémentée.
 *
 * @par Exemple
 * @snippet exemples/src/sdr/ex-sdr.cc ex_émetteur
 * @image html ex-emetteur.png width=1000px
 *
 *  @sa récepteur_création(), modulateur_création(), démodulateur_création()
 */
extern sptr<Émetteur> émetteur_création(const ÉmetteurConfig &ec);


//  * @image html emetteur.png "Architecture du bloc émetteur de trames" width=600px

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
extern ArrayXf doppler_distri(ArrayXd f, float fd, double fc);


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
extern float bruit_thermique(float bp, float T = 25);

/** @brief Paramétrage d'un émulateur de canal de propagation */
struct ECPConfig
{
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
extern sptr<Filtre<cfloat, cfloat, ECPConfig>> ecp_création(const ECPConfig &config);




/** @} */


/** @addtogroup telecom-eq
 *  @{
 */

/** @brief Création d'un égaliseur basé sur un ou des filtres RIF ajustés itérativement.
 *
 * <h3>Egaliseur RIF</h3>
 *
 * Création d'un égaliseur, échantillonné soit à la fréquence syumbole (@f$K = 1@f$), ou
 * échantillonné avec une période fractionnaire (@f$K > 1@f$).
 * Les structures d'égalisation suivantes sont possibles :
 *  - <b>Feed Forward Equalization (FFE)</b> <code>(structure = "dde")</code> :
 *       Un filtre RIF est réglé (à chaque période symbole) afin de minimiser le carré de l'erreur de sortie.
 *  - <b>Decision Feedback Equalization (DFE)</b> <code>(strucuture = "dfe")</code> :
 *      A la fois un filtre RIF direct (fonctionnant à @f$K\cdot f_{symb}@f$) et
 *      un filtre RIF de rétro-action (fonctionnant à la fréquence symbole) sur les décisions sont utilisés.
 *
 * Les fonctions d'erreur suivantes sont possibles :
 *  - <b>Basé sur la décision symbole</b> <code>(errf = "slicer")</code> : @f$E=(d-y)^2@f$. Avec cette fonction d'erreur, l'algorithme est aussi appelé LMS (Least Mean Square).
 *  - <b>Amplitude constante</b> (CMA / Constant Modulus Algorithm) <code>(errf = "cma")</code>
 *       @f$E=\left(R-|y|^2\right)^2@f$
 *
 * @image html ffe.png "Egalisation FFE" width=600px
 *
 * @param forme_onde        Forme d'onde (utilisée seulement si la fonction d'erreur est basée sur la décision symbole).
 * @param structure         Structure d'égalisation (directe : "ffe", ou avec rétro-actions : "dfe").
 * @param fonction_erreur   Fonction d'erreur (basé sur la décision symbole : "dec", ou constant modulus algorithm : "cma").
 * @param K                 Facteur de sur-échantillonnage.
 * @param α                 Taux de mise à jour pour l'algorithme LMS.
 * @param N1                Nombre de coefficients du filtre RIF d'égalisation.
 * @param N2                Nombre de coefficients du filtre RIF pour les rétro-actions (seulement si structure DFE).
 *
 *
 * @sa égaliseur_zfe()
 */
extern sptr<FiltreGen<cfloat>> égaliseur_rif_création(sptr<FormeOnde> forme_onde,
    const std::string &structure, const std::string &fonction_erreur,
    int K, float α, int N1, int N2);

/** @brief Calcul du filtre inverse par zéro-forçage.
 *
 * <h3>Calcul du filtre inverse par zéro-forçage</h3>
 *
 * Etant donné la réponse du canal @f$h@f$, cette fonction calcule les coefficients
 * d'un filtre RIF @f$g@f$, en essayant d'approximer
 * @f[
 * g\star h = \delta_d
 * @f],
 *
 * @f$d@f$ étant un délais global. Autrement dit, @f$g@f$ est un filtre inverse (au délais près) de @f$h@f$.
 *
 * @image html zfe.png "Egalisation ZFE" width=600px
 *
 * @note Cette fonction requiert de pouvoir mesurer la réponse du canal (par exemple en envoyant un signal de type impulsionnel côté émetteur).
 *
 * @warning
 *  - L'inversion n'est qu'approximative, le filtre inverse exact ayant une réponse impulsionnelle de support non borné.
 *  - Si la réponse du canal présente des zéros (ou des magnitudes faibles) dans le domaine fréquentielle, ce type d'égalisation n'est pas recommandée (amplification du bruit).
 *
 * @param h Réponse impulsionnelle du canal,
 * @param n Nombre de coefficients souhaités pour le filtre inverse.
 * @returns %Filtre RIF inverse (coefficients).
 *
 *  @par Exemple
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_eg_zfe
 *  @image html zfe-0.png "Réponses impulsionnelles (canal et du filtre d'égalisation)" width=800px
 *  <br/>
 *  @image html zfe-1.png "Réponses fréquentielles  (canal et du filtre d'égalisation)" width=800px
 *  <br/>
 *  @image html zfe-2.png "Exemple d'égalisation sur une flux NRZ" width=800px
 *
 * @sa égaliseur_création()
 */
extern ArrayXf égaliseur_zfe(IArrayXf h, int n);



extern Eigen::MatrixXf égaliseur_zfe_matrice(IArrayXf h, int n);


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
extern float capacite_canal_awgn(float snr, float B = 1);

/** @} */





/** @addtogroup telecom-crec
 *  @{
 */

/** @brief Structure de configuration pour une PLL */
struct PLLConfig
{
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
extern sptr<Filtre<float, float, RPLLConfig>> rpll_création(const RPLLConfig &config);

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
extern sptr<Filtre<cfloat, cfloat, PLLConfig>> cpll_création(const PLLConfig &config);

/** @} */


/** @addtogroup telecom-ber
 *  @{
 */


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
extern CmpBitsRes cmp_bits(const BitStream &b0, const BitStream &b1);

/** @brief Idem @ref cmp_bits(), avec gestion des ambiguité de phase M-PSK
 *
 *  <h3>Comparaison de chaines binaires (PSK)</h3>
 *
 */
extern CmpBitsRes cmp_bits_psk(const BitStream &b0, const BitStream &b1, int k);


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
extern void plot_eye(tsd::vue::Figure &f, const ArrayXf &x, float T);

/** @} */

/** @addtogroup telecom-snr
 *  @{
 */

/** @brief Interface abstraite pour un estimateur de SNR. */
struct EstimateurSNR
{
  /** @brief Calcul deux vecteurs (S et N) correspondant resp. aux énergies du signal et du bruit, à partir d'un signal bruité x. */
  virtual void step(const ArrayXcf &x, ArrayXf &S, ArrayXf &N) = 0;
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
extern sptr<EstimateurSNR> snr_Matzner(float γ = 0.1);

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
extern sptr<Filtre<float, float, AMConfig>> modulateurAM();

/** @brief TODO */
extern sptr<Filtre<cfloat, float, AMConfig>> demodulateurAM();

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
extern sptr<FiltreGen<cfloat,float>> discriminateur_fm();

/** @brief TODO */
extern sptr<Filtre<cfloat, cfloat, FMDemodConfig>> demodulateurFM();

/** @} */

/** @addtogroup telecom-codes-synchro
 *  @{
 */



/** @brief Génération d'un code à séquence maximale.
 *
 *  <h3>Génération d'un code à séquence maximale</h3>
 *
 *  @param n Longueur du registre à décalage (doit être compris entre 1 et 16).
 *
 *  Cette fonction génère un code binaire de longueur @f$m=2^n-1@f$, grâce à un registre à décalage
 *  et un polynôme primitif.
 *
 *  @par Exemple
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
extern BitStream code_Barker(int n);

/** @cond
 *  Renvoie un polynôme primitif de degré reglen.
 *  Le polynôme renvoyé est stocké "à l'envers", le LSB étant le coefficient
 *  de X^{n-1}, et le MSB celui de X^0. */
extern uint32_t polynome_primitif_binaire(int reglen);
/** @endcond */

/** @brief Calcul d'un polynôme primitif.
 *
 *  <h3>Calcul d'un polynôme primitif.</h3>
 *
 *  Cette fonction est utilisée pour la génération de codes à longueur maximale.
 *
 *  @param n Degré du polynôme (doit être compris entre 1 et 16).
 *
 *  @sa code_mls()
 */
extern Poly<int> polynome_primitif(int n);


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




