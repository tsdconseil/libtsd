#pragma once

#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include <vector>

namespace tsd::vue {
class Figures;
}

namespace tsd::telecom {




class BitStream;

/** @addtogroup telecom-ps
 *  @{
 */


/** @brief Spécification d'un filtre de mise en forme */
struct SpecFiltreMiseEnForme
{
  /** @brief Pas de filtre, des impulsions brutes sont transmises */
  static SpecFiltreMiseEnForme aucun();
  /** @brief Filtrage gaussien */
  static SpecFiltreMiseEnForme gaussien(float BT);
  /** @brief Filtrage "NRZ" (moyenne glissante) */
  static SpecFiltreMiseEnForme nrz();
  /** @brief Filtrage SRRC (racine de cosinus sur-élevé) */
  static SpecFiltreMiseEnForme srrc(float β);

  /** @brief Type de filtre. */
  enum Type
  {
    /** @brief Filtrage NRZ (simple répétition des symboles) */
    NRZ,
    /** @brief Emisssion d'un train d'impulsions brut */
    AUCUN,
    /** @brief Filtrage Gaussien */
    GAUSSIEN,
    /** @brief Racine de cosinus surélevé */
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

  /** @brief Spécification du filtre de mise en forme */
  SpecFiltreMiseEnForme filtre;

  /** @brief Nombre de symboles possibles */
  int M;

  /** @brief Nombre de bits par symbole (@f$\log_2(M)@f$) */
  int k;

  int cnt = 0;


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

}


