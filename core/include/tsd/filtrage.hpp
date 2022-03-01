#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"
#include "tsd/figure.hpp"

namespace tsd::vue {
  struct Figure;
  struct Figures;
}

namespace tsd::filtrage {


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
  extern void verifie_frequence_normalisee(float f, const std::string &msg = "");


  /** @} */

  /** @addtogroup filtrage-fenetres
    *  @{ */



  /** @cond private  */
  enum class Fenetre
  {
    AUCUNE = 0, // e.g. rectangulaire
    HANN,
    TRIANGLE,
    HAMMING,
    BLACKMAN,
    CHEBYCHEV
  };

  extern std::ostream& operator<<(std::ostream& ss, const Fenetre &t);


  extern std::string Fenetre2string(Fenetre f);

  extern ArrayXf fenetre(Fenetre type, int n, bool symetrique = true);
  /** @endcond */

  /** @brief Création d'une fenêtre sans paramètre (rectangulaire, Hann, Hamming, triangulaire ou Blackman)
   *
   *  <h3>Fenêtre sans paramètre</h3>
   *
   *  Cette fonction permet de créer une fenêtre simple, sans paramètre (en dehors de la dimension et du fait qu'elle soit symétrique ou non).
   *
   *  @param type Choix de la fenêtre : "re" (rectangulaire), "hn" (Hann), "hm" (Hamming), "tr" (triangulaire), "bm" (Blackman).
   *  @param n Nombre de points
   *  @param symetrique Si vrai, réalisation d'une fenêtre symétrique autour du point central (ce qui est adapté pour la conception
   *  d'un filtre), et sinon réalisation d'une fenêtre périodique (ce qui est adapté pour l'analyse spectrale).
   *  @return Vecteurs des coefficients de la fenêtre (vecteur de dimension @p n)
   *  @sa fenetre_kaiser(), fenetre_chebychev()
   *
   *  @par Exemple : création d'une fenêtre de Von Hann
   *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_fenetre
   *  @image html filtrage-fenetre.png width=600px
   */
  extern ArrayXf fenetre(const std::string &type, int n, bool symetrique = true);

  /** @brief Création d'une fenêtre de Chebychev.
   *
   *  <h3>Fenêtre de Chebychev</h3>
   *
   *  La fenêtre de Chebychev a la propriété d'avoir une ondulation d'amplitude constante (voir exemple ci-dessous).
   *  Le design est aussi très pratique, car on choisi l'ordre (le nombre de coefficients), et l'atténuation souhaitée,
   *  et c'est la largeur du lobe principal qui sert de variable d'ajustement.
   *
   *  @param n            Nombre de coefficients.
   *  @param atten_db     Atténuation en dB sur la bande coupée.
   *  @param symetrique  Si vrai, réalisation d'une fenêtre symétrique (adaptée pour la conception
   *  d'un filtre), sinon réalise une fenêtre périodique (adaptée pour l'analyse spectrale).
   *  @return Vecteurs des coefficients de la fenêtre (vecteur de dimension @p n)
   *  @par Exemple : création d'une fenêtre avec 60 dB d'atténuation
   *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_fenetre_cheby
   *  @image html filtrage-fenetre-cheby.png width=600px
   *
   *
   *  @sa design_rir_fen()
   */
  extern ArrayXf fenêtre_chebychev(int n, float atten_db, bool symetrique = true);


  /** @cond undoc */
  extern MatrixXf slepian_evec(int N, float B);
  /** @endcond */

  /** @brief Création d'une fenêtre de Slepian
   *
   *
   */
  extern ArrayXf fenêtre_slepian(int N, float B);


  // TODO : DOC
  struct FenInfos
  {
    /** @brief Atténuation pire lobe secondaire */
    float atten_ls;

    /** @brief Largeur lobe principal */
    float largeur_lp;

    /** @brief Vrai si la fenêtre est symétrique */
    bool symetrique;

    /** @brief Vrai si non symétrique */
    bool periodique;

    /** @brief Atténuation premier lobe secondaire */
    float atten_pls;

    tsd::vue::Figures fig;
  };

  // TODO : DOC
  extern FenInfos fenetre_analyse(const std::string &nom, const ArrayXf &x, bool do_plot = true);


  /** @brief Calcul du paramètre @f$\beta@f$ et de l'ordre d'un filtre de Kaiser.
   *
   *  <h3>Calcul des paramètres d'un filtre de Kaiser</h3>
   *
   *  @param atten_db Atténuation en dB dans la bande coupée (nombre positif)
   *  @param δf Largeur de la bande de transition (normalisée par rapport à la fréquence d'échantillonnage)
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
   *  @sa fenetre_kaiser(), fenetre_kaiser1()
   *
   *  @par Exemple
   *  @code
   *  // Paramètres pour une atténuation de 60 dB,
   *  // et une bande de transition de 1 dixième de la fréquence d'échantillonnage.
   *  auto [β, n] = kaiser_param(60, 0.1);
   *  @endcode
   */
  extern std::tuple<float, int> kaiser_param(float atten_db, float δf);

  /** @brief Création d'une fenêtre de Kaiser.
   *
   *  <h3>Création d'une fenêtre de Kaiser</h3>
   *
   *  @param atten_db Atténuation en dB sur la bande coupée
   *  @param δf       Largeur du lobe principal (en fréquence normalisée)
   *  @param symetrique Si vrai, réalisation d'une fenêtre symétrique (adaptée pour la conception
   *  d'un filtre), sinon réalise une fenêtre périodique (adaptée pour l'analyse spectrale).
   *  @return Vecteurs des coefficients de la fenêtre
   *  @sa fenetre_kaiser1(), kaiser_param()
   *
   *  @par Exemple : création d'une fenêtre avec 60 dB d'atténuation
   *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_fenetre_kaiser
   *  @image html filtrage-fenetre-kaiser.png width=600px
   */
  extern ArrayXf fenêtre_kaiser(float atten_db, float δf, bool symetrique = true);

  /** @brief Création d'une fenêtre de Kaiser (d'après paramètre de forme @f$\beta@f$).
   *
   * <h3>Création d'une fenêtre de Kaiser</h3>
   *
   * Cette fonction permet de créer la fenêtre d'après paramètre de forme @f$\beta@f$.
   *
   * @param n     Nombre de coefficients.
   * @param β     Paramètre de forme de la fenêtre.
   * @param symetrique Si vrai, réalisation d'une fenêtre symétrique (adaptée pour la conception
   *  d'un filtre), sinon réalise une fenêtre périodique (adaptée pour l'analyse spectrale).
   * @return Vecteurs des coefficients de la fenêtre
   *
   *  @sa fenetre_kaiser(), kaiser_param() */
  extern ArrayXf fenêtre_kaiser1(int n, float β, bool symetrique = true);



  /** @} */



  /** @addtogroup filtrage-analyse
   *  @{ */

/** @brief Magnitude d'une filtre RIF ou RII.
 *
 *  <h3>Magnitude d'une filtre RIF ou RII</h3>
 *
 *  Cette fonction calcule de la magnitude (valeur absolue) de la réponse fréquentielle d'un filtre (RIF ou RII) :
 *  @f[
 *  y_k = \left|H(f_k)\right|
 *  @f]
 *
 *  @param h Fonction de transfert à analyser
 *  @param npts Résolution fréquentielle
 *  @return Un tuple de deux vecteurs : le vecteur de fréquences @f$f_k@f$
 *  (normalisées, entre 0 et 0,5), et la magnitude de la réponse @f$y_k@f$.
 *
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_frmag
 *  @image html frmag.png width=600px
 *
 *  @sa frgroup(), repfreq(), frphase()
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frmag(const FRat<T> &h, int npts = 1024);



/** @brief Réponse fréquentielle d'un filtre RIF ou RII.
 *
 * <h3>Réponse fréquentielle d'un filtre RIF ou RII</h3>
 *
 *  Cette fonction calcule de la réponse fréquentielle (complexe) :
 *  @f[
 *  y_k = H(f_k)
 *  @f]
 *
 *  @param h  Fonction de transfert à analyser
 *  @param fr Vecteur de fréquences (normalisées)
 *  @return   Réponse fréquentielle @f$H@f$
 *
 *  @sa frmag(), frphase(), frgroup()
 */
template<typename T>
  ArrayXcf repfreq(const FRat<T> &h, const ArrayXf &fr);



/** @brief Réponse impulsionnelle. */
template<typename T>
  ArrayXf repimp(const FRat<T> &h, int npts = -1);


/** @brief Phase d'un filtre RIF ou RII.
 *
 *  <h3>Phase d'un filtre RIF ou RII</h3>
 *
 *  Cette fonction calcule de la phase de la réponse fréquentielle d'un filtre (RIF ou RII) :
 *  @f[
 *  y_k = \arg H(f_k)
 *  @f]
 *
 *  @param h Fonction de transfert à analyser
 *  @param npts Résolution fréquentielle
 *  @return Un tuple de deux vecteurs : le vecteur de fréquences (normalisées, entre 0 et 0,5), et la phase de la réponse (en radians).
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frphase(const FRat<T> &h, int npts = 1024);

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
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_frgroup
 *  @image html frgroup.png width=600px
 *
 */
template<typename T> std::tuple<ArrayXf, ArrayXf> frgroup(const FRat<T> &h, int npts = 1024);




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
 * @snippet exemples/src/filtrage/ex-filtrage.cc exemple_analyse
 * @image html filtrage-analyse.png width=1000px
 *
 * @sa affiche_filtre() (pour un affichage plus simple : seulement réponse fréquentielle et impulsionnelle)
 *
 */
template<typename T>
  tsd::vue::Figures analyse_filtre(const FRat<T> &h, float fe = 1.0f);


/** @brief Affichage des réponses impulsionnelles et fréquentielles d'un filtre.
 *
 *  <h3>Affichage des réponses impulsionnelles et fréquentielles d'un filtre</h3>
 *
 * Cette fonction créée une nouvelle figure, et trace les réponses fréquentielles
 * et temporelles du filtre.
 *
 *  @param h      Fonction de transfert à analyser (ou vecteur de coefficients).
 *  @param fe     Fréquence d'échantillonnage (optionnel).
 *  @return       La figure créée.
 *
 * @par Exemple
 * @snippet exemples/src/filtrage/ex-filtrage.cc exemple_affiche
 * @image html filtrage-affiche.png width=1000px
 *
 * @sa analyse_filtre() (pour un affichage plus complet)
 */
template<typename T>
  tsd::vue::Figures affiche_filtre(const FRat<T> &h, float fe = 1.0f);


extern FenInfos filtre_pb_analyse(const ArrayXf &h);
extern FenInfos filtre_pb_analyse(int ncoefs, const ArrayXf &fr, const ArrayXf &mag, bool do_plot = true);

/** @cond undoc */
template<typename T> ArrayXf repimp(const Vecteur<T> &h, int npts = -1);
extern ArrayXcf repfreq(const ArrayXf &h, const ArrayXf &fr);
template<typename T> ArrayXf repfreq(const Vecteur<T> &h, int npts = 1024);
template<typename T> std::tuple<ArrayXf, ArrayXf> frmag(const Vecteur<T> &h, int npts = 1024);
template<typename T> std::tuple<ArrayXf, ArrayXf> frphase(const Vecteur<T> &h, int npts = 1024);
template<typename T> std::tuple<ArrayXf, ArrayXf> frgroup(const Vecteur<T> &h, int npts = 1024);
extern tsd::vue::Figures analyse_filtre(const ArrayXf &h, float fe = 1.0f);
extern tsd::vue::Figures affiche_filtre(const ArrayXf &h, float fe = 1.0f);
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
 *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_plz
 *  @image html filtrage-plz.png width=600px
 *  */
template<typename T>
  void plot_plz(tsd::vue::Figure &fig, const FRat<T> &h);



/** @brief Réponse en amplitude d'un filtre RIF symétrique ou anti-symétrique (phase linéaire).
 *
 *  <h3>Réponse en amplitude d'un filtre RIF</h3>
 *
 *  Cette fonction calcule la réponse en amplitude @f$A(\omega)@f$, pour un filtre RIF
 *  réel à phase linéaire (les coefficients doivent être symétriques ou anti-symétriques autour du point central).
 *
 *  Pour un filtre symétrique, la réponse en amplitude est :
 *  @f[
 *  A(\omega) = h_M + 2 \sum_{n=0}^{M-1} h_n \cos((M-n)\omega)
 *  @f]
 *
 *  @param L Résolution fréquentielle.
 *  @param h Coefficients du filtre.
 *  @param symetrique Vrai si les coefficients sont symétriques autour de @f$(N-1)/2@f$ (filtre RIF de type I ou II).  Sinon un filtre de type III ou IV est supposé.
 *  @returns Un tuple de deux vecteurs : le vecteur de fréquences (normalisées, entre 0 et 0,5),
 *  et la réponse en amplitude.
 *
 *  @warning La réelle symétrie (ou anti-symétrie) des coefficients n'est pas vérifiée !
 *
 *  @sa frmag(), frphase()
 */
extern std::tuple<ArrayXf, ArrayXf> rifamp(const Eigen::ArrayXf &h, int L = 1024, bool symetrique = true);


/** @brief Calcul du retard d'une filtre RIF à phase linéaire.
 *
 *  <h3>Calcul du retard d'une filtre RIF à phase linéaire</h3>
 *
 *  Cette fonction renvoie le retard, en nombre d'échantillons, introduit par un filtre
 *  RIF réel à phase linéaire, c'est-à-dire dont les coefficients sont symétriques (ou anti-symétriques)
 *  autour du point central :
 *  @f[
 *  \tau = \frac{N-1}{2}
 *  @f]
 *
 *  @param N Nombre d'échantillons.
 *  @returns Délais @f$\tau@f$ du filtre.
 *
 */
extern float rif_delais(int N);

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
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_hilbert
 *  @image html design-rif-hilbert.png width=600px
 *
 *  @sa hilbert(), hilbert_transformeur()
 */
extern ArrayXf design_rif_hilbert(int n, const std::string &fenetre = "hn");



/** @brief Spécification d'un Biquad */
struct BiquadSpec
{
  /** @brief Type de filtre biquad (voir @ref design_biquad()) */
  enum Type
  {
    /** @brief Passe-bas */
    PASSE_BAS,
    /** @brief Passe-haut */
    PASSE_HAUT,
    /** @brief Passe-bande */
    PASSE_BANDE,
    /** @brief Coupe-bande */
    COUPE_BANDE,
    /** @brief Résonnance */
    RESONATEUR,
    /** @brief Amplification des basses fréquences (amplification constante). */
    PLATEAU_BF,
    /** @brief Amplification des hautes fréquences (amplification constante). */
    PLATEAU_HF
  } type;

  /** @brief Fréquence de coupure (ou centrale) normalisée (entre 0 et 0,5). */
  float f;

  /** @brief Facteur de qualité (note : pour @f$Q>1/\sqrt(2)\sim 0{,71}@f$, il y aura une résonnance).
   *  Non utilisé pour low et high shelf. */
  float Q;

  /** @brief Utilisé uniquement pour les filtres à résonnance, et shelf. */
  float gain_dB = 1;
};

extern std::ostream& operator<<(std::ostream &ss, const BiquadSpec &t);

/** @brief Design filtre biquad.
 *
 * <h3>Design filtre biquad</h3>
 *
 * Pour la descriptionc complète, voir @ref design_biquad().
 *
 * @param spec Spécification (type, fréquence de coupure, facteur de qualité, etc.)
 *
 * @sa design_biquad()
 */
extern FRat<float> design_biquad(const BiquadSpec &spec);

/** @brief Design filtre biquad.
 *
 * <h3>Design filtre biquad</h3>
 *
 * Ces filtres RII du second ordre sont adaptés de prototypes analogiques via la transformée bilinéaire.
 *
 * Les prototypes analogiques sont les suivant (pour une pulsation de coupure de 1 radian/s, et @f$Q@f$ étant le facteur de qualité) :
 *
 * - %Filtre passe-bas :
 * @f[
 * H(s) = \frac{1}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - %Filtre passe-haut :
 * @f[
 * H(s) = \frac{s^2}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - %Filtre passe-bande :
 * @f[
 * H(s) = \frac{s/Q}{s^2+\frac{1}{Q}s+1}
 * @f]
 * - %Filtre coupe-bande :
 * @f[
 * H(s) = \frac{s^2+1}{s^2+\frac{1}{Q}s+1}
 * @f]
 *
 *
 * @param type    Type de filtre ("lp", "hp", "bp", "sb", ...).
 * @param f       Fréquence de coupure (ou centrale pour les filtres passe ou stoppe bande) normalisée, entre 0 et 0,5.
 * @param Q       Facteur de qualité (note : pour @f$Q>1/\sqrt(2)\sim 0{,71}@f$, il y aura une résonnance).
 * @param gain_dB Gain, en dB, pour les filtres de type résonnance ou plateau.
 *
 * @par Exemple : filtres passe-bas, avec différentes valeurs pour le facteur de qualité
 * @snippet exemples/src/filtrage/ex-filtrage.cc ex_biquad_lp
 * @image html ex-biquad-pb.png width=800px
 *
 * @par Bibliographie
 * - <i>Cookbook formulae for audio equalizer biquad filter coefficients,</i> Robert Bristow-Johnson,
 *      https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html,
 * - <i>F0 and Q in filters, Mini tutorial,</i> Analog Devices,
 */
extern FRat<float> design_biquad(const std::string type, float f, float Q, float gain_dB = 0);

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
 *  Le filtre est d'abort conçu dans le domaine analogique (transformée de Laplace),
 *  puis converti en filtre digital (transformée en z) grâce à la
 *  transformée bilinéaire.
 *
 * @param n           Ordre du filtre.
 * @param type        Type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...)
 * @param prototype   "butt", "cheb1", "cheb2" ou "ellip"
 * @param fc          Fréquence de coupure normalisée (entre 0 et 0,5)
 * @param δ_bp        Ondulation maximale en decibels dans la bande passante (utilisé seulement pour un filtre de Chebychev de type I ou un filtre elliptique).
 * @param δ_bc        Atténuation minimale en decibels dans la bande coupée (utilisé seulement pour un filtre de Chebychev type II ou un filtre elliptique).
 * @return h          Fonction de transfert (digitale)
 *
 * @par Exemple
 * @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_riia
 * @image html design-riia.png width=800px
 *
 * @sa trf_bilineaire(), @ref filtre_sois(), filtre_rii()
 */
extern FRat<cfloat> design_riia(int n, const std::string &type,
    const std::string &prototype, float fc, float δ_bp = 0.1f, float δ_bc = 60);



/** @brief Design par échantillonnage fréquentiel.
 *
 *  <h3>Design par échantillonnage fréquentiel</h3>
 *
 * Cette technique permet d'approximer avec un filtre RIF de @f$n@f$ coefficients (@f$n@f$ étant impair)
 * une réponse fréquentielle arbitraire, passée en paramètre.
 * La réponse fréquentielle doit être donnée sous la forme d'un tableau @f$d@f$ de @f$m=\frac{n+1}{2}@f$ éléments réels,
 * de type :
 * @f[
 * d_k = H(f_k),\ f_k = k \cdot \frac{1}{2m-1},\ k = 0,\dots, m-1
 * @f]
 *
 *
 * Les @f$f_k@f$ peuvent être obtenus par la fonction @ref design_rif_freq_freqs().
 *
 * @note
 * Notez que @f$m@f$ éléments de réponse fréquentielle permettent de spécifier de manière unique un filtre
 * réel avec @f$n=2m-1@f$ coefficients (en effet, chaque élement de la réponse fréquentielle, sauf le premier,
 * est utilisé une deuxième fois, pour les fréquences négatives).
 * Par conséquent, si @f$n\neq 2m-1@f$, la réponse fréquentielle souhaitée est, avant calcul du filtre, ré-échantillonnée (par interpolation linéaire) avec @f$m'@f$ valeurs de tel sorte que
 * @f$n=2m'-1@f$ (dans tous les cas, le nombre de coefficients doit être impair).
 *
 *
 *  @param n     Ordre du filtre (doit être impair).
 *  @param d     Vecteur définissant la réponse fréquentielle souhaitée (sur les fréquences positives).
 *  @returns     Vecteur des coefficients du filtre (dimension = @f$n@f$).
 *
 * @par Exemple
 * @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_freq
 * @image html design-rif-freq.png width=800px
 *
 * @sa design_rif_freq_freqs()
 */
extern ArrayXf design_rif_freq(int n, const ArrayXf &d);

/** @brief Calcul des @f$m@f$ valeurs de fréquences utilisée pour le design par échantillonnage fréquentiel.
 *
 * <h3>Calcul des @f$m@f$ valeurs de fréquences utilisée pour le design par échantillonnage fréquentiel</h3>
 *
 * Cette fonction retourne le vecteur de fréquences :
 * @f[
 *  f_k = k \cdot \frac{1}{2m-1},\ k = 0,\dots, m-1
 * @f]
 * où @f$m=(n+1)/2@f$ est le nombre de points à renseigner.
 *
 * @param n Ordre du filtre (doit être impair).
 *
 * @sa design_rif_freq()
 */
extern ArrayXf design_rif_freq_freqs(int n);

/** @cond private */
extern tsd::vue::Figures design_rif_freq_analyse(int n, const ArrayXf &d);
/** @endcond */

/** @brief Design équiripple / Remez.
 *
 *  <h3>Design équiripple / Remez</h3>
 *
 *  @param n      Ordre du filtre
 *  @param d      Réponse souhaitée
 *  @param w      Poids de pondération (doit être de même longueur que d)
 *  @param debug  Si vrai, affiche des plots des calculs intermédiaires.
 *  @returns      Vecteur des coefficients du filtre (vecteur de dimension n)
 */
extern ArrayXf design_rif_eq(int n, IArrayXf d, IArrayXf w, bool debug = false);

extern ArrayXf design_rif_eq2(int n, IArrayXf D, IArrayXf W, bool debug = false);

extern ArrayXf design_rif_eq(int n, const std::vector<SpecFreqIntervalle> &spec, bool debug = false, bool cheby_mode = true);

/** @brief Sinus cardinal normalisé avec fréquence de coupure paramétrable.
 *
 *  <h3>Sinus cardinal normalisé avec fréquence de coupure paramétrable</h3>
 *
 *  @f[
 *  y(t) = \frac{\sin(2 \pi t f_c)}{\pi  t}
 *  @f]
 *
 *  La transformée de Fourier de cette fonction est une porte de largeur @f$\pm f_c@f$,
 *  ce qui en fait donc le protype idéal (mais non réalisable) pour un filtre passe-bas.
 *
 *  @param t Point d'échantillonnage temporel (en nombre d'échantillons)
 *  @param fc Fréquence de coupure
 *  @returns Valeur de @f$y@f$
 *
 *  @sa sinc(), design_rif_fen() */
extern float sinc2(float t, float fc);

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
extern float sinc(float t);

/** @brief Design RIF par sinus-cardinal fenêtré.
 *
 * <h3>Design RIF par sinus-cardinal fenêtré</h3>
 *
 * @param n       Ordre du filtre,
 * @param type    Type de filtre ("lp" pour passe-bas, "hp" pour passe-haut, ...),
 * @param fc      Fréquence de coupure normalisée,
 * @param fen     Type de fenêtre ("hn", "hm", "tr" ou "re"), voir fonction @ref fenetre(),
 * @param fc2     Deuxième fréquence de coupure (uniquement pour les filtres passe-bande ou stoppe-bande),
 * @returns       Vecteur des coefficients du filtre
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
 *  Ce filtre présente l'intérêt de ne générer aucune interférence inter-symboles
 *  (si les symboles sont séparés de @f$T= \frac{1}{2 f_c}@f$).
 *
 *  @param n  ordre du filtre
 *  @param β  facteur de dépassement
 *  @param fc fréquence de coupure normalisée
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_cs
 *  @image html design-rif-cs.png width=800px
 *
 *  @sa design_rif_rcs() */
extern ArrayXf design_rif_cs(int n, float β, float fc);

/** @brief Design d'un filtre en racine de cosinus sur-élevé.
 *
 *  <h3>Design d'un filtre en racine de cosinus sur-élevé</h3>
 *
 *  @param n  ordre du filtre
 *  @param β  facteur de dépassement
 *  @param fc fréquence de coupure normalisée
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_rcs
 *  @image html design-rif-rcs.png width=800px
 *
 *  @sa design_rif_cs(), design_rif_rcs1()
 */
extern ArrayXf design_rif_rcs(int n, float β, float fc);

/** @brief Design d'un filtre en racine de cosinus sur-élevé (1)
 *
 *  <h3>Design d'un filtre en racine de cosinus sur-élevé (1)</h3>
 *
 *  Cette fonction est équivalente à @ref design_rif_rcs(), si ce n'est
 *  qu'au lieu de la fréquence de coupure @f$f_c@f$, c'est le facteur de sur-échantillonnage
 *  qui est passé en paramètre (@f$\textrm{OSF} = \frac{1}{2 f_c}@f$).
 *
 *  @param n   Ordre du filtre.
 *  @param β   Facteur de dépassement.
 *  @param osf Rapport de sur-échantillonnage (@f$\frac{f_e}{f_{symb}}@f$).
 *  @param nrm Type de normalisation des coefficients ('s' pour somme = 1, 'e' pour somme des carrés égale à 1).
 *
 *  @sa design_rif_rcs()
 */
extern ArrayXf design_rif_rcs1(int n, float β, float osf, char nrm = 's');

/** @brief Coefficients d'une approximation RIF d'un filtre gaussien (non fenêtré).
 *
 * <h3>%Filtre Gaussien (approximation RIF)</h3>
 *
 *  @param n   Nombre de coefficients du filtre.
 *  @param σ   Ecart-type (en nombre d'échantillons).
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_gaussien
 *  @image html design-rif-gaussien.png width=800px
 *
 *  @sa design_rif_gaussien_telecom()
 */
extern ArrayXf design_rif_gaussien(int n, float σ);


/** @brief Coefficients d'une approximation RIF d'un filtre gaussien, pour une modulation GFSK.
 *
 *  <h3>%Filtre Gaussien avec porte</h3>
 *
 *  Calcul de la convolution d'un filtre Gaussien et d'un filtre en moyenne glissante, de largeur égale
 *  au facteur de sur-échantillonnage.
 *  L'écart-type du filtre gaussien est calculé en fonction du produit BT.
 *
 *  Les paramètres de cette fonction sont adaptés à la réalisation
 *  d'un filtre de mise en forme pour une modulation GFSK.
 *
 *  @param n   Nombre de coefficients du filtre
 *  @param BT  Produit Bande-passante - Temps
 *  @param osf Facteur de sur-échantillonnage (rapport entre la fréquence d'échantillonnage et la fréquence symbole)
 *  @return    Vecteur des coefficients du filtre
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rif_gaussien_telecom
 *  @image html design-rif-gaussien-telecom.png width=800px
 *
 *  @sa design_rif_gaussien()
 */
extern ArrayXf design_rif_gaussien_telecom(int n, float BT, int osf);


/** @brief Calcul l'écart-type (relatif à la période symbole) en fonction du BT */
extern float design_rif_gaussien_telecom_BT_vers_sigma(float BT);





/** @brief Calcul d'un filtre RIF issu d'une concaténation série de deux filtres RIF.
 *
 * <h3>Mise en série de deux filtres RIF</h3>
 *
 * Cette fonction calcule les coefficients d'un filtre RIF qui a la même réponse
 * que la mise en série de deux filtres RIF spécifiés en paramètres :
 * @f[
 * h = h_1 \star h_2
 * @f]
 *
 * Le nombre de coefficients est déterminé de manière à être le minimum requis,
 * soit @f$n = n_1 + n_2 - 1@f$.
 *
 * @param h1 Coefficients du premier filtre (@f$n_1@f$ coefficients)
 * @param h2 Coefficients du deuxième filtre (@f$n_2@f$ coefficients)
 * @return   Coefficients du filtre concaténé  (@f$n_1+n_2-1@f$ coefficients)
 *
 */
extern ArrayXf design_rif_prod(const ArrayXf &h1, const ArrayXf &h2);

/** @brief Paramètres principaux d'un filtre CIC */
struct CICConfig
{
  /** @brief Decimation ratio */
  int R = 1;
  /** @brief Number of integrators / differentiators */
  int N = 1;
  /** @brief Design parameter (typically M=1) */
  int M = 1;
};

/** @brief Fonction de transfert théorique d'un filtre CIC

    <h3>Fonction de transfert théorique d'un filtre CIC</h3>

    @param config Paramètres principaux (voir @ref CICConfig)
    @returns CIC transfert function (not taking into account the decimation)
    This function computes the theorical transfert function of a CIC filter, when one does not look at the decimation effect. The CIC filter responses is defined as:
           @f[H(z) = \frac{1}{R^N}\left(1+z^{-1}+\cdots+z^{-(R-1)}\right)^N@f]
     (e.g. a cascade of @f$N@f$ moving average filters, each of identical length @f$R@f$.
 **/
extern FRat<float> design_cic(const CICConfig &config);


struct CICComp
{
  ArrayXf h;
  tsd::vue::Figure
    fig_reponse_ideale,
    fig_comp_rimp,
    fig_spectre_global,
    fig_spectre_bf;
};

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
extern CICComp design_cic_comp(const CICConfig &config, float Fin, int R2, float fc, int ncoefs);


/** @brief Fonction de transfert pour un bloqueur DC.
 *
 *  <h3>Fonction de transfert bloqueur DC</h3>
 *
 *  Calcul de la fonction de transfert suivante :
 *  @f[
 *  H(z) = \frac{1 - z^{-1}}{1 - \alpha \cdot z^{-1}}
 *  @f]
 *  Avec
 *  @f[
 *  \alpha = \frac{\sqrt{3} - 2 \sin(\pi f_c)}{\sin(\pi f_c) + \sqrt{3} \cos(\pi f_c)}
 *  @f]
 *
 *  @param fc Fréquence de coupure normalisée (0 - 0,5)
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_design_bloqueur_dc
 *  @image html bloqueur-dc-resp.png width=1000px
 *
 *  @par Bibliographie
 *  <i>The DC Blocking Filter,</i> J.M. de Freitas, 2007
 *
 *  @sa filtre_dc()
 */
extern FRat<float> design_bloqueur_dc(float fc);

/** @brief Fonction de transfert d'un filtre exponentiel
 *
 *  <h3>Fonction de transfert d'un filtre exponentiel</h3>
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
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_design_rii1
 *  @image html design-rii1.png width=1000px
 *
 *  @sa rii1_coef()
 */
extern FRat<float> design_rii1(float fc);

/** @brief Design filtre RII d'ordre 1 (d'après la fréquence de coupure).
 *
 *  <h3>Design filtre RII d'ordre 1</h3>
 *
 *  Cette fonction calcule le coefficient (facteur d'oubli) d'un filtre RII du premier ordre
 *  (filtre exponentiel), en fonction de la fréquence de coupure souhaitée.
 *
 *  @f[
 *    y_{n+1} = (1-\gamma) y_n + \gamma x_n
 *            = y_n + \gamma (x_n - y_n)
 *  @f]
 *
 *  @f[
 *    \gamma = 1 - e^{-2\pi f_c}
 *  @f]
 *
 *  @param fc Normalized cut-off frequency, in [0, 0.5]
 *  @returns Coeficient @f$\gamma@f$ for the IIR filter
 *
 *  @sa design_rii1(), rii1_tc_vers_coef()
 *
 */
extern float rii1_coef(float fc);

/** @brief Design filtre RII d'ordre 1 (d'après la constante de temps souhaitée).
 *
 *  <h3>Design filtre RII d'ordre 1</h3>
 *
 *  Cette fonction calcule le coefficient (facteur d'oubli) d'un filtre RII du premier ordre
 *  (filtre exponentiel), en fonction de la constante de temps souhaitée.
 *
 *  @f[
 *    y_{n+1} = (1-\gamma) \cdot y_n + \gamma \cdot x_n
 *            = y_n + \gamma \cdot (x_n - y_n)
 *  @f]
 *
 *  @f[
 *    \gamma = 1 - e^{-1/τ}
 *  @f]
 *
 *  @param τ Constante de temps, en nombre de symboles.
 *  @returns Coefficient d'oubli @f$\gamma@f$ du filtre RII
 *
 */
extern float rii1_tc_vers_coef(float τ);

extern float rii1_coef_vers_tc(float γ);

/** @brief Compute cut-off frequency, from forget factor of first order IIR filter.
 *
 *  @param γ Facteur d'oubli
 *  @returns Fréquence de coupure
 *
 *  @sa rii1_coef() */
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
  Tableau<T> forme_polyphase(const Eigen::Ref<const Vecteur<T>> x, unsigned int M);

/** @brief Compute the standard form from the polyphase representation
 *
 *  <h3>Calcul de la forme standard à partir de la forme polyphase</h3>
 *
 *  @param X Forme polyphase (tableau 2d)
 *  @returns Forme normale (tableau 1d)
 *
 *  @sa forme_polyphase()
 */
template<typename T>
  Vecteur<T> iforme_polyphase(const Eigen::Ref<const Tableau<T>> X);

/** @brief Transformée bilinéaire : conversion transformée de Laplace vers transformée en z.
 *
 * <h3>Transformée bilinéaire</h3>
 *
 *  La transformée bilinéaire permet d'approximer une fonction de transfert analogique (transformée
 *  de Laplace) grâce à une fonction de tranfert digitale (transformée en z).
 *
 *  Le calcul de la transformée en z est fait en approximant la
 *  variable @f$s@f$ de la transformée de Laplace par :
 *  @f[
 *  s \mapsto 2 f_e \cdot \frac{1 - z^{-1}}{1 + z^{-1}}
 *  @f]
 *
 *  @param ha Transformée de Laplace du système (fraction rationnelle)
 *  @param fe Fréquence d'échantillonnage
 *  @returns  Transformée en z (fraction rationnelle)
 *
 *  @sa fd_vers_fa(), fa_vers_fd()
 */
extern FRat<cfloat> trf_bilineaire(const FRat<cfloat> &ha, float fe);


extern float ωa_vers_ωd(float ωa, float fe);
extern float ωd_vers_ωa(float ωd, float fe);

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


/** @brief Retarde le signal d'entrée d'un nombre entier d'échantillons.
 *
 *  <h3>Ligne à retard</h3>
 *
 *  Ce filtre produit autant d'échantillons en sortie qu'il y en a en entrée.
 *  Les échantillons précédant le premier échantillon d'entrée sont supposés nuls.
 *
 *  @param n Délais entier (@f$n \geq 0@f$)
 *
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_ligne_a_retard
 *  @image html filtrage-ligne-a-retard.png width=800px
 *  */
template<typename T>
  sptr<FiltreGen<T>> ligne_a_retard(int n);

struct HilbertTransformeurConfig
{
  int ntaps = 63;
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
 * @return          %Filtre float -> cfloat
 *
 * @sa design_rif_hilbert(), hilbert(), hilbert_tfd()
 */
extern sptr<Filtre<float, cfloat, HilbertTransformeurConfig>>
  hilbert_transformeur(int n = 31, const std::string &fenetre = "hn");



/** @brief Implémentation directe d'une filtre RIF (Réponse Impulsionnelle Finie).
 *
 * <h3>Implémentation directe d'une filtre RIF (Réponse Impulsionnelle Finie)</h3>
 *
 *  Implémenté suivant l'équation :
 *  @f[
      y_n = \sum h_k x_{n-k}
    @f]

 *  @param h Vecteur des coefficients du filtre
 *  @tparam T Type des données à traiter (float, cfloat, etc.)
 *  @tparam Tc Type des coefficients
 *  @return %Filtre T -> T
 *
 *  @note Si le filtre a un nombre important de coefficients, utilisez plutôt @ref filtre_rif_fft(), qui sera plus efficace
 *  (filtrage dans le domaine fréquentiel).
 *
 *  @sa filtre_rif_fft()
 */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rif(const Eigen::Ref<const Vecteur<Tc>> h);



/** @brief Filtre identité.
 *
 *  <h3>Filtre identité</h3>
 *
 *  Ce filtre laisse le signal inchangé.
 */
template<typename T>
  sptr<FiltreGen<T>> filtre_id();


/** @brief Décimateur 1:R
 *
 *  <h3>Décimateur 1:R</h3>
 *
 *  Ce "filtre" supprime @f$R-1@f$ échantillons tous les @f$R@f$.
 */
template<typename T>
  sptr<FiltreGen<T>> decimateur(int R);


/** @brief Implémentation efficace d'un filtre RIF par FFT
 *
 *  <h3>Implémentation fréquentielle d'un filtre RIF</h3>
 *
 *  Ce bloc réalise un filtrage RIF via la technique OLA (Ovelap-And-Add).
 *  La complexité est donc bien moindre que l'implémentation standard si le nombre de coefficients @f$M@f$
 *  est important (de l'ordre de @f$\log M@f$ opérations par échantillon au lieu de @f$M@f$).
 *
 *  @param h Vecteur des coefficients du filtre
 *  @tparam T Type de données à traiter
 *  @return %Filtre T -> T
 *
 *  @note Notez que cette technique introduit un délais un peu plus important que l'implémentation temporelle.
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
 *  @param h    Fonction de transfert (fraction rationnelle).
 *  @tparam T   Type de données à traiter.
 *  @tparam Tc  Type des coefficients.
 *  @return     %Filtre T -> T.
 *
 *  @sa @ref filtre_sois()
 *
 *  @warning Si l'ordre du filtre est important, du fait des erreurs de troncature, cette implémentation a de fortes chances
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
    @param config Paramètres principaux (voir @ref CICConfig)
 *  @param mode 'd' for decimation or 'u' for upsampling.
 *  @tparam T Type d'entrée / sortie
 *  @tparam Ti Type pour les calculs interne (int, int64_t, ...)
 *
 *  @note Pour le type interne, il faut absolument choisir un type entier, car la façon dont le filtre est
 *  implémenté fait que les calculs ne fonctionneront pas avec un type flottant.
 *
 *
 *  @par Exemple pour l'interpolation
 *  @snippet exemples/src/filtrage/ex-cic.cc exemple_cic_upsampling
 *  @image html filtrage-cic-interpolation.png width=800px
 *  Notez les repliements du signal utile sur le spectre.
 *  <br>
 *  @par Exemple pour la décimation
 *  @snippet exemples/src/filtrage/ex-cic.cc exemple_cic_decimation
 *  @image html filtrage-cic-decimation.png width=800px
 */
template<typename T, typename Ti>
  sptr<FiltreGen<T>> filtre_cic(const CICConfig &config, char mode = 'd');


/** @brief Résultat de l'analyse d'un filtre CIC */
struct CICAnalyse
{
  CICConfig config;
  /** Input sample frequency */
  float Fin;
  /** Fréquence où mesurer l'atténuation */
  float Fint;
  FRat<float> h;
  int nbits;
  ArrayXf fr, mag;
  tsd::vue::Figure fig_repliement, fig_spectre_global, fig_spectre_bf;
};

/** @brief This function computes and shows the frequency response of a CIC filter and
 *  then analyse the aliasing that occurs after decimation.
 *
 * @param config Paramètres principaux (voir @ref CICConfig)
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
extern CICAnalyse cic_analyse(const CICConfig &config, float Fin, float Fint = 0);



/** @brief Frequency response of a CIC filter
 *
 * @param config Paramètres principaux (voir @ref CICConfig)
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
extern ArrayXf cic_freq(const CICConfig &config, const ArrayXf &f);




/** @brief Type de structure pour l'implémentation d'un filtre RII */
typedef enum
{
  /** Numérateur, puis dénominateur */
  FormeDirecte1,
  FormeDirecte1Transposee,
  /** Dénominateur, puis numérateur */
  FormeDirecte2,
  FormeDirecte2Transposee
} RIIStructure;

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
 *  @return %Filtre réel vers réel
 *
 *  La forme directe 2 est légérement plus efficace que la forme 1. Cependant, si les coefficients du filtre
 *  sont amenés à changer en cours de fonctionnement, il vaut mieux utiliser la première forme, car la forme 2
 *  risque de générer des discontinuités.
 *
 *  @warning Il est recommandé que la fonction de transfert passée en entrée soit sous la forme pôles / zéros,
 *  et non pas sous la forme de la liste des coefficients, car avec cette dernière représentation, la position
 *  des pôles et zéros (et donc la réponse fréquentielle) est très instable (très sensible aux erreurs de troncature
 *  lors de la quantification des coefficients). Notez que la fonction @ref design_riia() renvoie bien une fonction
 *  de transfert sous la forme pôles / zéros.
 *
 *  @sa filtre_rii(), design_riia()
 */
template<typename T>
  sptr<FiltreGen<T>> filtre_sois(const FRat<cfloat> &h, RIIStructure structure = FormeDirecte2);


template<typename T>
  sptr<FiltreGen<T>> filtre_sois(const FRat<float> &h, RIIStructure structure = FormeDirecte2);

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
 */
template<typename T>
  sptr<FiltreGen<T>> filtre_rii1(float γ);


/** @brief %Filtre pour la suppression du DC (basses fréquences).
 *
 * <h3>%Filtre pour la suppression du DC</h3>
 *
 *  Implémente le filtre :
 *  @f[
 *  y_k = x_k - x_{k-1} + \alpha y_{k-1}
 *  @f]
 *
 *  @param fc Fréquence de coupure normalisée (entre 0 et 0.5)
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc exemple_filtre_dc
 *  @image html bloqueur-dc-ex.png width=1000px
 *
 *  @par Bibliographie
 *  <i>The DC Blocking Filter,</i> J.M. de Freitas, 2007
 *
 *  @sa design_bloqueur_dc()
 */
template<typename T>
  sptr<FiltreGen<T>> filtre_dc(float fc);


/** @brief Moyenne glissante.
 *
 *  <h3>Moyenne glissante</h3>
 *
 *
 *  Ce filtre est une implémentation optimisée d'une moyenne glissante :
 *  @f[
 *  y_n = \frac{1}{N}\cdot\sum_{k=0}^{K-1} x_{n-k}
 *  @f]
 *  @remark Ce filtre est implémenté efficacement grâce à un intégrateur et à un peigne de profondeur @f$K@f$ :
 *  @f[
 *   y_n = y_{n-1} + \frac{x_n - x_{n-K}}{K}
 *  @f]
 *
 *  @param K profondeur de la moyenne
 *  @tparam T Type de données à traiter
 *  @tparam Tacc Type à utiliser pour l'accumulateur
 *
 */
template<typename T, typename Tacc>
  sptr<FiltreGen<T>> filtre_mg(int K);



/** @} */



/** @addtogroup filtrage-fini
  *  @{ */

/** @brief Filtrage d'un signal de durée finie par un filtre défini par sa fonction de transfert.
 *
 *  <h3>Filtrage d'un signal par un filtre défini par sa fonction de transfert</h3>
 *
 *  Cette  fonction ne fonctionne que sur un signal de durée finie. Pour filtrer des données
 *  reçues au fil de l'eau, il faut utiliser une des structure avec contexte (voir
 *  @ref filtre_rif(), @ref filtre_sois(), etc.).
 *
 *  @param h Fonction de transfert (RIF ou RII) ou coefficients (RIF seulement) du filtre à appliquer.
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
 *  @sa filtfilt()
 */
template<typename T, typename Tc>
Vecteur<T> filtrer(const FRat<Tc> &h, const Vecteur<T> &x)
{
  if(h.est_fir())
  {
    auto f = filtre_rif<Tc, T>(h.numer.coefs.reverse());
    return f->step(x);
  }
  auto f = filtre_sois<T>(h);
  return f->step(x);
}

template<typename T, typename Tc>
Vecteur<typename T::Scalar> filtrer(const Vecteur<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  auto f = filtre_rif<Tc, typename T::Scalar>(h);
  return f->step(x);
}

/** @brief Filtrage zéro-phase (bi-directionnel)
 *
 *  <h3>Filtrage zéro-phase (bi-directionnel)</h3>
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
 *
 *  @sa filtrer()
 */
template<typename T, typename Tc>
  Vecteur<typename T::Scalar> filtfilt(const Vecteur<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  return filtrer(h, filtrer<T,Tc>(h,x).reverse()).reverse();
}


template<typename T, typename Tc>
  Vecteur<T> convol(const Eigen::Ref<const Vecteur<Tc>> h, const Eigen::Ref<const Vecteur<T>> x)
{
  // TODO : via FFT si la dimension de h en vaut la peine
  auto f = filtre_rif<Tc, T>(h);
  return f->step(x);
}

/** @brief Calcul du signal analytique (via un filtrage RIF).
 *
 *  <h3>Calcul du signal analytique (via un filtrage RIF)</h3>
 *
 *  @param x Signal d'entrée
 *  @param ncoefs Nombre de coefficients pour l'approximation RIF
 *  @return Signal analytique (complexe)
 *
 *  @sa design_rif_hilbert(), hilbert_tfd(), hilbert_transformeur() */
extern ArrayXcf hilbert(IArrayXf x, int ncoefs = 127);

/** @brief Calcul du signal analytique (via la TFD)
 *
 *  <h3>Calcul du signal analytique (via la TFD)</h3>
 *
 *  Cette fonction met tout simplement à zéro les fréquences négatives du signal,
 *  dans le domaine fréquentiel.
 *  Contrairement à @ref hilbert(), cette technique n'introduit aucun délais dans le domaine temporel.
 *
 *  @warning
 *  Du fait des hypothèses sous-jacentes à la la TFD, des artefacts peuvent apparaitre
 *  aux bords du signal.
 *
 *  @param x Signal d'entrée
 *  @return Signal analytique (complexe)
 *  @sa design_rif_hilbert(), hilbert(), hilbert_transformeur() */
extern ArrayXcf hilbert_tfd(IArrayXf x);

/** @} */


/** @addtogroup filtrage-rythme-itrp
  *  @{ */

/** @brief Interface générique pour un interpolateur */
template<typename T>
struct Interpolateur
{
  /** @brief Longueur de la fenêtre (ligne à retard). */
  int K = 0;

  /** @brief */
  float delais = 0;

  /** @brief Description de l'interpolateur */
  std::string nom;

  virtual ~Interpolateur(){}

  /** @brief Interpolation (méthode abstraite).
   *
   *  <h3>Interpolation</h3>
   *  Cette méthode abstraite réalise l'interpolation, à partir
   *  des échantillons de la ligne à retard @f$x@f$.
   *
   *  @f[
   *    y = x\left(k + K/2 + \tau [K]\right)
   *  @f]
   *
   *  @param x        Fenêtre circulaire contenant les  de @f$K@f$ derniers points d'entrée.
   *  @param k        Index de l'échantillon le plus ancien dans la ligne à retard (@f$k+K-1@f$ est donc l'index de l'échantillon le plus récent).
   *  @param τ        Retard fractionnaire, entre 0 et 1.
   *  @returns        Valeur interpoléee (@f$y@f$).
   *
   */
  virtual T step(const Vecteur<T> &x, int k, float τ) = 0;
};

/** @brief %Interpolateur générique à base filtre RIF.
 *
 *  <h3>%Interpolateur générique à base filtre RIF</h3>
 *
 *  @f[
 *    y =\sum_{i=0}^{K-1} h_i \cdot x_{i+k[K]}
 *  @f]
 *
 */
template<typename T>
struct InterpolateurRIF: Interpolateur<T>
{
  virtual ~InterpolateurRIF(){}

  /** @brief Calcul des coefficients d'un filtre RIF en fonction du délais souhaité.
   *
   *  <h3>Calcul des coefficients d'un filtre RIF</h3>
   *
   *  Cette méthode abstraite doit être redéfinie par la classe dérivée,
   *  de manière à calculer les coefficients de filtrage pour un délais fractionnaire @f$\tau@f$ donné (entre 0 et 1).
   *  Les coefficients renvoyés doivent être ceux d'un filtre avec un délais exactement égal à :
   *  @f[
   *  d = K/2 + \tau
   *  @f]
   *
   *  @f$K@f$ étant le nombre de coefficients du filtre.
   *
   *  @sa Interpolateur, itrp_sinc()
   */
  virtual ArrayXf coefs(float τ) = 0;

  // k : index de l'échantillon le plus ancien
  T step(const Vecteur<T> &x, int k, float τ)
  {
    ArrayXf h = coefs(τ);
    tsd_assert(h.rows() == this->K);
    T res = 0;
    for(auto i = 0; i < this->K; i++)
      res += h(i) * x((i + k) % this->K);
    return res;
  }
};


/** @cond undoc
 *  (voir à supprimer) */


/* Interpolateur demi-bande */
template<typename T> sptr<FiltreGen<T>> filtre_allpass_ups();

/* Décimateur demi-bande */
template<typename T> sptr<FiltreGen<T>> filtre_allpass_decim();

/** @endcond */

/** @brief %Interpolation par splines cardinales.
 *
 *  <h3>%Interpolation par splines cardinales</h3>
 *
 *  @sa InterpolateurRIF, itrp_lineaire(), itrp_lagrange(), itrp_sinc()
 */
template<typename T>
  sptr<InterpolateurRIF<T>> itrp_cspline();

/** @brief %Interpolation linéaire (équivalente à Lagrange degré 1).
 *
 *  <h3>%Interpolation linéaire (équivalente à Lagrange degré 1)</h3>
 *
 *  @sa InterpolateurRIF, itrp_cspline(), itrp_lagrange(), itrp_sinc()
 */
template<typename T>
  sptr<InterpolateurRIF<T>> itrp_lineaire();

/** @brief %Interpolateur de Lagrange (la fonction sinus cardinal est interpolée par un polynôme).
 *
 *  <h3>%Interpolateur de Lagrange</h3>
 *  L'interpolation polynomiale (Lagrange) consiste à interpoler le signal grâce à un polynôme d'ordre donné,
 *  ce qui revient, en terme de filtre RIF, à interpoler la fonction sinus cardinal par un polynôme.
 *
 *
 *  @sa InterpolateurRIF, itrp_cspline(), itrp_lineaire(), itrp_sinc()
 */
template<typename T>
  sptr<InterpolateurRIF<T>> itrp_lagrange(unsigned int degré);


/** @brief Structure de configuration pour un interpolateur Sinc. */
struct InterpolateurSincConfig
{
  /** @brief Nombre de coefficients (dimension de la ligne à retard). */
  int ncoefs = 31;

  /** @brief Nombre de pas pour la LUT. */
  int nphases = 256;

  /** @brief Fréquence de coupure normalisée. */
  float fcut = 0.5;

  /** @brief Type de fenêtre */
  std::string fenetre = "hn";
};

/** @brief %Interpolateur à sinus cardinal fenêtré.
 *
 *  <h3>%Interpolateur à sinus cardinal fenêtré</h3>
 *
 *  %Interpolateur RIF dont les coefficients sont ceux d'un sinus cardinal fenêtré,
 *  de fréquence de coupure @f$f_c@f$ et avec un retard de @f$K/2 + \tau@f$ échantillons :
 *
 *  @f[
 *  h_k = \textrm{sinc}\left(k - K/2 - \tau, f_c\right) \cdot w_{k-\tau}
 *  @f]
 *
 *  La fonction sinc étant définie ici comme la TF inverse d'une porte fréquentielle de largeur @f$\pm f_c@f$.
 *
 *  @sa InterpolateurRIF, itrp_cspline(), itrp_lineaire(), itrp_lagrange()
 */
template<typename T>
  sptr<InterpolateurRIF<T>> itrp_sinc(const InterpolateurSincConfig &config = InterpolateurSincConfig());




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
 * @tparam T  Type d'entrée / sortie
 *
 *
 * @sa filtre_rif_ups()
 */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rif_decim(const Eigen::Ref<const Vecteur<Tc>> h, int R);


/** @brief A documenter : filtre RIF optimisé pour filtrage demi-bande */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rif_demi_bande(const Eigen::Ref<const Vecteur<Tc>> c);

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
 *  @return un filtre (type T vers T), le flux de sortie étant upsamplé d'un facteur @f$R@f$ par rapport au flux d'entrée.
 *
 *
 *  @sa filtre_rif_ups_délais(), filtre_rif_decim(), filtre_rif()
 */
template<typename Tc, typename T = Tc>
  sptr<FiltreGen<T>> filtre_rif_ups(const Eigen::Ref<const Vecteur<Tc>> h, int R);


/** @brief Calcul du retard d'une filtre RIF avec upsampling.
 *
 * <h3>Calcul du retard d'une filtre RIF avec upsampling</h3>
 *
 * Si le nombre de coefficients, @f$K@f$ est un multiple de @f$R@f$ :
 * @f[
 *  \tau = \frac{K - 1}{2}
 * @f]
 *
 * Sinon :
 * @f[
 *  \tau = \frac{K - 1}{2} + R - (K [R])
 * @f]
 *
 *
 * @sa filtre_rif_ups()
 */
extern float filtre_rif_ups_délais(int nc, int R);



/** @brief Adatateur de rythme, ratio arbitraire.
 *
 *  <h3>Adatateur de rythme (ratio arbitraire)</h3>
 *
 *  Cette fonction renvoi un bloc de ré-échantillonnage, permettant de changer le rythme d'un signal reçu au fil de l'eau.
 *  L'implémentation est basée sur une cascade de décimateurs (si le ratio est inférieur à 1) ou d'interpolateurs (si le ratio est supérieur à 1) demi-bandes,
 *  suivis d'un interpolateur de ratio arbitraire.
 *
 *  @param ratio Ratio de décimation / interpolation (rapport entre les fréquences d'échantillonnage de sortie et d'entrée).
 *  @sa resample(), filtre_rif_ups(), filtre_rif_decim()
 **/
template<typename T>
  sptr<FiltreGen<T>> filtre_reechan(float ratio);

/** @brief Interpolation d'un ratio arbitraire (calcul au fil de l'eau)
 *
 * <h3>Interpolation (ratio arbitraire)</h3>
 *
 *  @param ratio Ratio de décimation / interpolation (rapport entre les fréquences d'échantillonnage de sortie et d'entrée).
 *  @param itrp  Pointeur vers un interpolateur générique
 *
 *  @sa itrp_cspline(), itrp_lineaire(), itrp_lagrange(), itrp_sinc()
 */
template<typename T> sptr<FiltreGen<T>> filtre_itrp(float ratio, sptr<Interpolateur<T>> itrp = itrp_cspline<T>());





/** @brief Type d'interpolation pour un signal d'entrée échantillonné de manière aléatoire. */
enum class InterpOption
{
  LINEAIRE = 0,
  CSPLINE
};


/** @brief Interpolation d'un signal échantilloné de manière aléatoire.
 *
 * <h3>Interpolation d'un signal échantilloné de manière aléatoire</h3>
 *
 * Cette fonction permet de rééchantilloner un signal de manière réguliére, alors
 * que l'on ne connait les valeurs qu'en un ensemble de points, pas forcément équidistants.
 * Plus précisément, étant donné un ensemble de @f$N@f$ points @f$(x_k,f(y_k))@f$,
 * on calcule un ensemble de valeurs @f$f(x'_k)@f$, les @f$x'_k@f$ étant équirépartis entre @f$x_0@f$ et @f$x_{N-1}@f$.
 *
 * @param x    Vecteur des points d'échantillonnage en entrée (doit être un vecteur croissant).
 * @param y    Vecteur des valeurs connues de la fonction à interpoler.
 * @param x2   Vecteur des points d'échantillonnage souhaités (doit être un vecteur croissant).
 * @param mode Type d'interpolation souhaitée (linéaire ou spline naturelles).
 *
 *  @par Exemple
 *  @snippet exemples/src/filtrage/ex-filtrage.cc ex_itrp_irreg
 *  @image html itrp-irreg.png width=800px
 *
 */
template<typename T>
Vecteur<T> interp(const ArrayXf &x, const Vecteur<T> &y, const ArrayXf &x2, InterpOption mode = InterpOption::LINEAIRE);




/** @} */




}






