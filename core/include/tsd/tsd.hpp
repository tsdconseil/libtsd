#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */


#include <complex>
#include <Eigen/Core>
#include <memory>
#include <cstdarg>
#include <vector>
#include <numbers>
#include <random>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/color.h>
#include <functional>

#ifndef FMT_RUNTIME
#if (FMT_VERSION >= 80000)
# define FMT_RUNTIME(AA) fmt::runtime(AA)
#else
# define FMT_RUNTIME(AA) AA
#endif
#endif


namespace cutils::model
{
  class Node;
}

namespace tsd {
extern void msg_impl2(const char *fn, const int ligne, int niveau, const char *fonction, const std::string &str);


typedef std::function<void(const char *fn, const int ligne, int niveau, const char *fonction, const std::string &str)> logger_t;

template<typename ... Ts>
  void msg_impl(const char *fn, const int ligne, int niveau, const char *fonction, const std::string &format_str, Ts &&... args)
{
  auto s = fmt::format(FMT_RUNTIME(format_str), args...);
  msg_impl2(fn, ligne, niveau, fonction, s);
}

extern void set_logger(logger_t logger);

}

#define msg(...)          tsd::msg_impl(__FILE__, __LINE__, 1, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_majeur(...)   tsd::msg_impl(__FILE__, __LINE__, 2, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_avert(...)    tsd::msg_impl(__FILE__, __LINE__, 3, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_erreur(...)   tsd::msg_impl(__FILE__, __LINE__, 4, __PRETTY_FUNCTION__, __VA_ARGS__)
#define echec(...)        tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, __VA_ARGS__)

/** @brief Display an error message in the log before asserting */
#define tsd_assert_msg(AA, ...)  if(!(AA)) {tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, __VA_ARGS__);}
#define tsd_assert(AA)           if(!(AA)) {tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, "{}", "Echec assertion : " #AA ".");}

#define dsp_assertion(AA, ...) {printf("\033[1;37;41mEchec assertion : " #AA ".\nFichier: %s, ligne: %d.\033[0m", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); fflush(0); {*((char *) 0) = 5;} assert(AA);}

/** @brief Returns -1 if the condition is not true (and add an error message in the log) */
#define tsd_check(AA, ...) if(!(AA)){dsp_assertion(AA, __VA_ARGS__); return -1;}

#ifdef DEBUG_MODE
# define dsp_assert_safe tsd_assert
# define dsp_check_safe tsd_check
#else
# define dsp_assert_safe(...)
# define dsp_check_safe(...)
#endif


namespace tsd
{

/** @addtogroup tsd
 *  @{ */


  /** @brief Valeur de @f$\pi@f$, en pr??cision flottante 64 bits. */
  static const auto pi = std::numbers::pi;

  /** @brief Valeur de @f$\pi@f$, en pr??cision flottante 64 bits. */
  static const auto ?? = pi;

  /** @brief Valeur de @f$\pi@f$, en pr??cision flottante 32 bits. */
  static const auto ??_f = std::numbers::pi_v<float>;


  using namespace std::complex_literals;

  /** @brief Raccourci pour nombre complexe flottant 32 bits */
  typedef std::complex<float> cfloat;

  /** @brief Raccourci pour nombre complexe flottant 64 bits */
  typedef std::complex<double> cdouble;

  /** @cond private */
  using Eigen::ArrayXf;
  using Eigen::ArrayXd;
  using Eigen::ArrayXXf;
  using Eigen::ArrayXXcf;
  using Eigen::ArrayXi;
  using Eigen::MatrixXf;
  using Eigen::MatrixXcf;
  using Eigen::VectorXf;
  using Eigen::VectorXcf;


  using fmt::format;

  /** @brief Vecteur de complexes 32 bits */
  using Eigen::ArrayXcf;

  /** @endcond */

  template<typename t>
    using Tableau = Eigen::Array<t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

  /** @brief Raccourcis pour un vecteur colonne de donn??es de type g??n??rique (param??tre de template t) */
  template<typename t>
    using Vecteur = Eigen::Array<t, Eigen::Dynamic, 1, Eigen::ColMajor>;

  /** @brief Vecteur flottant en entr??e (const) */
  using IArrayXf   = const Eigen::Ref<const ArrayXf>;

  /** @brief Tableau complexe flottant en entr??e (const) */
  using IArrayXXcf = const Eigen::Ref<const ArrayXXcf>;

  /** @brief Vecteur complexe flottant en entr??e (const) */
  using IArrayXcf  = const Eigen::Ref<const ArrayXcf>;

  /** @brief Vecteur flottant en sortie */
  using OArrayXf   = Eigen::Ref<ArrayXf>;
  using IArrayXXf  = const Eigen::Ref<const ArrayXXf>;
  using OArrayXXf  = Eigen::Ref<ArrayXXf>;

  /** @brief Tableau de boul??ens */
  using ArrayXXb = Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic>;

  /** @brief Vecteur colonne de boul??ens */
  using ArrayXb  = Eigen::Array<bool,Eigen::Dynamic,1>;

  using IArrayXXb = const Eigen::Ref<const ArrayXXb>;
  using IArrayXb  = const Eigen::Ref<const ArrayXb>;

  /** @brief Raccourci pour d??finir un pointeur partag?? */
  template<typename T>
    using sptr = std::shared_ptr<T>;


  /** @brief Concat??nation verticale de deux vecteurs.
   *
   * <h3>Concat??nation verticale</h3>
   *
   * Cette fonction r??alise l'assemblage de deux vecteurs colonnes :
   * @f[
   * c = \left(
   * \begin{array}{c}
   * a_0\\
   * a_1\\
   * \vdots\\
   * b_{n-1}\\
   * b_0\\
   * b_1\\
   * \vdots\\
   * b_{m-1}\\
   * \end{array}
   * \right)
   * @f]
   *
   * (en supposant que a et b aient respectivement @f$n@f$ et @f$m@f$ ??l??ments).
   * L'op??rateur | a ??t?? red??fini pour r??aliser la m??me op??ration.
   *
   * @param a Premier vecteur
   * @param b Deuxi??me vecteur
   * @returns Concat??nation verticale
   *
   * @par Example
   * @code
   *  ArrayXf a = ArrayXf::Ones(1);
   *  ArrayXf b = ArrayXf::Zero(5);
   *  ArrayXf c = vconcat(a, b); // (ou ArrayXf c = a | b)
   * @endcode
   *
   */
  template<typename T1, typename T2>
    auto vconcat(const Eigen::ArrayBase<T1> &a, const Eigen::ArrayBase<T2> &b)
  {
    Vecteur<typename T1::Scalar> c(a.rows()+b.rows());
    if(a.rows() > 0)
      c.head(a.rows()) = a;
    if(b.rows() > 0)
      c.tail(b.rows()) = b;
    return c;
  }

  template<typename T1, typename T2>
    auto operator |(const Eigen::ArrayBase<T1> &a, const Eigen::ArrayBase<T2> &b)
  {
    return vconcat(a, b);
  }

  /** @cond undoc */
  /* Calcule la longueur d'un vecteur colonne */
  template<typename T>
    auto length(const Eigen::ArrayBase<T> &x)
  {
    return x.rows();
  }
  /** @endcond */

  /** @brief Rotation d'un vecteur
   *
   *  <h3>Rotation d'un vecteur</h3>
   *
   *  @param x Vecteur colonne d'entr??e
   *  @param d Nombre entier de pas (positif ou n??gatif)
   *  @returns Vecteur des donn??es d??cal??e de @f$d@f$ pas (modulo la dimension du vecteur) :
   *  @f[
   *   y_k = x_{k + d\ [N]}
   *  @f]
   *  @f$N@f$ ??tant la dimension du vecteur.
   *
   *  @par Exemple
   *  @code
   *  ArrayXf x = linspace(0, 5, 6);
   *  ArrayXf y = rotation_vec(x, 2);
   *  ArrayXf yref(6);
   *  yref << 2, 3, 4, 5, 0, 1;
   *  tsd_assert(y.isApprox(yref));
   *  @endcode
   */
  template<typename T>
  Vecteur<T> rotation_vec(const Vecteur<T> &x, int d)
  {
    auto n = x.rows();
    Vecteur<T> y(n);

    for(auto i = 0; i < n; i++)
      y(i) = x((i + d + n) % n);

    return y;
  }

  /** @brief Diff??rence entre 2 ??lements successifs d'un vecteur
   *
   *  <h3>Diff??rence entre 2 ??lements successifs d'un vecteur</h3>
   *
   *  @returns Un vecteur de longueur @f$n-1@f$ :
   *  @f[
   *  y_k = x_{k+1} - x_k,\quad k = 0,\dots, n-2
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  ArrayXf x = randn(10);
   *  ArrayXf y = diff(x);
   *  ArrayXf yref = x.tail(9) - x.head(9);
   *  tsd_assert(y.isApprox(yref));
   *  @endcode
   *
   *  @sa cumsum()
   */
  template<typename D>
    auto diff(const Eigen::ArrayBase<D> &x)
  {
    auto n = x.rows();
    return x.tail(n-1) - x.head(n-1);
  }

  /** @brief Somme cumul??e d'un vecteur
   *  <h3>Somme cumul??e d'un vecteur</h3>
   *
   *  Calcule un vecteur de dimension identique et ??gale ?? l'acculation des
   *  ??chantillons du vecteurs d'entr??e :
   *  @f[
   *  y_n = \sum_{k=0}^n x_k,\quad k = 0\dots N-1
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  ArrayXf x     = linspace(0, 99, 100);
   *  ArrayXf y     = cumsum(x);
   *  ArrayXf yref  = x * (x + 1) / 2;
   *  tsd_assert(y.isApprox(yref));
   *  @endcode
   *
   *  @sa diff()
   *
   */
  template<typename D>
      auto cumsum(const Eigen::ArrayBase<D> &x)
  {
    auto n = x.rows();
    Vecteur<typename D::Scalar> y(n);
    typename D::Scalar cs = 0;
    for(auto i = 0; i < n; i++)
    {
      cs += x(i);
      y(i) = cs;
    }
    return y;
  }

  /** @brief Corrige les sauts de phases
   *
   *  <h3>Correction des sauts de phase</h3>
   *
   *  @param x Vecteur contenant typiquement des angles
   *  @param r D??finit la classe d'??quivalence sur les valeurs de @f$x@f$
   *  @returns Un vecteur @f$y@f$ tel que @f$y_k = x_k + k \cdot r @f$, et avec le moins de discontinuit??s possibles (plus pr??cis??ment, aucune discontinuit?? sup??rieure ?? @f$r/2@f$ en valeur absolue).
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc exemple_unwrap
   *  @image html unwrap.png width=800px
   *
   *  @sa modulo_pm_??(), modulo_2??()
   *
   */
  extern Eigen::ArrayXf d??plie_phase(const ArrayXf &x, float r = 2*??);

  /** @brief Recherche d'??l??ments vrais dans un tableau boul??en.
   *
   *  <h3>Recherche d'??l??ments vrais dans un tableau boul??en</h3>
   *
   *  @param x Tableau de valeurs boul??ennes
   *  @returns Vecteur des index pour lequels les ??l??ments du tableau sont vrais.
   *
   *  @par Exemple
   *  @code
   *  ArrayXf x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
   *  auto idx = trouve(x >= 0); // idx = {2, 3, 4}
   *  @endcode
   *
   *  @sa trouve_premier()
   */
  extern std::vector<int> trouve(IArrayXb x);

  /** @brief Recherche de l'indice du premier ??l??ment vrai dans un tableau de boul??ens.
   *
   *  <h3>Recherche de l'indice du premier ??l??ment vrai dans un tableau de boul??ens.</h3>
   *
   *
   *  @par Exemple
   *  @code
   *  ArrayXf x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
   *  auto idx = trouve_premier(x >= 0); // idx = 2
   *  @endcode
   *
   *  @sa trouve()
   */
  extern int trouve_premier(IArrayXb x);

  extern ArrayXcf polar(const ArrayXf &??);
  extern ArrayXcf polar(const ArrayXf &??, const ArrayXf &??);


  // Equivalent ?? Scilab x(dec:pas:dec+pas*n)
  // Fonction ?? supprimer ?
  inline auto subarray1d(IArrayXf &x, unsigned int dec, unsigned int n, unsigned int pas)
  {
    auto p = (float *) x.data();
    return Eigen::Map<ArrayXf, 0, Eigen::InnerStride<>>(p + dec, n, Eigen::InnerStride<>(pas));
  }

  // Equivalent de : y = x(1:pas:$)
  /** @brief Sous-??chantillonnage d'un vecteur colonne
   *
   *  <h3>Sous-??chantillonnage d'un vecteur colonne</h3>
   *
   *  A partir d'un vecteur @f$(x_k), k=0\dots n-1@f$, renvoie
   *  un vecteur sous-??chantillonn?? d'un facteur @f$R@f$ :
   *  @f[
   *  y_k = x_{kR},\quad k =0 \dots \left\lfloor N/R\right\rfloor -1
   *  @f]
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc ex_sousech
   *  @image html sousech.png width=800px
   *
   *  @sa surech()
   *
   */
  template<typename T>
   auto sousech(const Vecteur<T> &x, int R)
   {
    return Eigen::Map<const Vecteur<T>,0,Eigen::InnerStride<>>(
        x.data(), x.rows() / R, Eigen::InnerStride<>(R));
   }

  /** @brief Sur-??chantillonnage d'un vecteur colonne
   *
   *  <h3>Sur-??chantillonnage d'un vecteur colonne</h3>
   *
   *  A partir d'un vecteur @f$(x_k), k=0\dots n-1@f$, renvoie
   *  un vecteur sur-??chantillonn?? d'un facteur @f$R@f$, par l'insertion de z??ros :
   *  @f[
   *  y_k = \begin{cases} x_{k/R} & \textrm{si } k\textrm{ est multiple de } R,\\
   *   0 & \textrm{sinon.}\end{cases}
   *  @f]
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc ex_surech
   *  @image html surech.png width=800px
   *
   *  @sa sousech()
   *
   */
  template<typename D>
      auto surech(const Eigen::ArrayBase<D> &x, int R)
  {
    using T = typename D::Scalar;
    Vecteur<T> y = Vecteur<T>::Zero(x.rows()*R);
    auto map = Eigen::Map<Vecteur<T>,0,Eigen::InnerStride<>>(y.data(), y.rows() / R, Eigen::InnerStride<>(R));
    map = x;
    return y;
  }

  /*template<typename T>
   Vecteur<T> surech(const Vecteur<T> &x, int R)
   {
    Vecteur<T> y = Vecteur<T>::Zero(x.rows()*R);
    auto map = Eigen::Map<Vecteur<T>,0,Eigen::InnerStride<>>(y.data(), y.rows() / R, Eigen::InnerStride<>(R));
    map = x;
    return y;
   }*/

  /** @brief Conversion lin??aire vers decibels
   *
   *  <h3>Conversion lin??aire vers decibels</h3>
   *
   *  @param x Valeur en lin??aire
   *  @returns Valeur en dB :
   *  @f[
   *    y = 10 \log_{10}(x)
   *  @f]
   *
   *  @sa db2pow()
   */
  template<typename T>
    auto pow2db(const T &x)
  {
    if constexpr (std::is_integral_v<T>)
        return 10 * std::log10(x * 1.0);
    else
      return 10 * std::log10(x);
  }


  /** @brief Conversion decibels vers lin??aire
   *
   *  <h3>Conversion decibels vers lin??aire</h3>
   *
   *  @param x Valeur en dB
   *  @returns Valeur en lin??aire :
   *  @f[
   *    y = 10^{x/10}
   *  @f]
   *
   *  @sa pow2db()
   */
  template<typename T>
    auto db2pow(const T &x)
  {
    if constexpr (std::is_integral_v<T>)
      return std::pow(10.0, x/10.0);
    else
      return std::pow((T) 10, x/10);
  }


  /** @brief Retourne la plus petite puissance de 2 sup??rieure ou ??gale ?? i.
   *
   * <h3>Prochaine puissance de 2</h3>
   *
   * Calcule la plus petite puissance de 2 sup??rieure ou ??gale ?? @f$i@f$ :
   * @f[
   *  y = \min_{j\geq i,\ j = 2^k} j
   * @f]
   *
   * @par Exemple
   * @code
   * tsd_assert(prochaine_puissance_de_2(3) == 4);
   * @endcode
   *
   */
  extern int prochaine_puissance_de_2(unsigned int i);




  /** @brief Ajoute des z??ros ?? la fin d'un des deux vecteurs de mani??re ?? ce qu'ils aient la m??me dimension.
   *
   *  <h3>Compl??tion de vecteurs avec des z??ros</h3>
   *
   *  Un des deux vecteurs est compl??t?? avec des z??ros de mani??re ?? ce que les deux vecteurs
   *  r??sultants aient la m??me longueur.
   *
   *  @param x  Premier vecteur
   *  @param y  Deuxi??me vecteur
   *  @param p2 Si vrai, arrondi la dimension des vecteurs de sortie ?? la prochaine puissance de 2.
   *  @return Un tuple de deux vecteurs de m??mes dimensions, et commen??ant respectivement par x et y (des z??ros sont ajout??s ?? la fin).
   *
   *  @par Exemple
   *  @code
   *  auto x = linspace(0,4,5);
   *  // x = {0, 1, 2, 3, 4};
   *  auto y = linspace(0,3,4);
   *  // y = {0, 1, 2, 3};
   *
   *  auto [x2,y2] = pad_zeros(x,y);
   *  // x2 = {0, 1, 2, 3, 4};
   *  // y2 = {0, 1, 3, 4, 0};
   *  @endcode
   *
   */
  template<typename D1, typename D2>
    std::tuple<Vecteur<typename D1::Scalar>, Vecteur<typename D2::Scalar>> pad_zeros(
        const Eigen::ArrayBase<D1> &x,
        const Eigen::ArrayBase<D2> &y,
        bool p2 = false)
  {
    Vecteur<typename D1::Scalar> b1p, b2p;

    auto n1 = x.rows(), n2 = y.rows();

    auto n3 = std::max(n1, n2);

    if(p2)
      n3 = prochaine_puissance_de_2(n3);

    b1p = ArrayXf::Zero(n3);
    b2p = ArrayXf::Zero(n3);
    b1p.head(n1) = x;
    b2p.head(n2) = y;

    return {b1p, b2p};
  }



  /** @brief Calcule le carr?? d'un scalaire */
  template<typename T>
    T carr??(const T &v) {return v * v;}

  /** @brief Structure abstraite pour une classe configurable.
   *  @tparam C Type de configuration */
  template<typename C>
  struct Configurable
  {
    virtual ~Configurable(){}

    /** @brief Configuration de la structure
     *  @param c Valeur de la configuration ?? appliquer
     *  @returns entier non nul en cas d'erreur. */
    int configure(const C &c)
    {
      if(callback_modif)
        callback_modif(c);
      int res = configure_impl(c);
      config = c;
      return res;
    }

    // Retourne 0 si support?? et rend "configure_impl" inutile
    virtual int configure(const cutils::model::Node &n)
    {
      //config.deserialise(n);
      //configure_impl(config);
      return -1;
    }

    /** @brief M??thode abstraite ?? impl??menter par la classe d??riv??e */
    virtual int configure_impl(const C &c) = 0;

    /** @brief Lecture de la configuration actuelle. */
    const C &lis_config() const {return config;}

    /** @brief Callback optionnelle qui sera appel??e ?? chaque changement de la configuration. */
    std::function<void (const C &nv_config)> callback_modif;

  protected:
    C config;
  };

  struct Void{};

  /** @brief Structure abstraite pour un puit de donn??es g??n??rique */
  template<typename Te>
  struct SinkGen
  {
    virtual ~SinkGen(){}

    /** @brief Traitement d'un bloc de donn??es
     *  @param x Signal d'entr??e (r??f??rence vers un vecteur Eigen) */
    virtual void step(const Eigen::Ref<const Vecteur<Te>> x) = 0;
  };

  /** @brief Structure abstraite pour un puit de donn??es g??n??rique et configurable */
  template<typename Te, typename Tc = Void>
  struct Sink : Configurable<Tc>, SinkGen<Te>
  {
    virtual ~Sink(){}
  };

  /** @brief Structure abstraite pour une source de donn??es
   *  @tparam Ts Type de donn??es ?? produire (float, cfloat, ...)
   *  @tparam C  Type de configuration */
  template<typename Ts>
  struct SourceGen
  {
    virtual ~SourceGen(){}
    /** @brief Production de n ??chantillons
     *  @param n Nombre d'??chantillons ?? produire
     *  @returns Vecteur de n ??chantillons s*/
    virtual Vecteur<Ts> step(int n) = 0;
  };

  /** @brief Structure abstraite pour source de donn??es, configurable.
   *  @tparam Ts Type de donn??es ?? produire (float, cfloat, ...)
   *  @tparam C  Type de configuration */
  template<typename Ts, typename C = Void>
  struct Source : SourceGen<Ts>, Configurable<C>
  {
    virtual ~Source(){}
  };

  /** @brief Structure abstraite pour un filtre de donn??es
   *  @tparam Te Input type
   *  @tparam Ts Output type */
  template<typename Te, typename Ts = Te>
  struct FiltreGen
  {
    virtual ~FiltreGen(){}

    /** @brief Traitement d'un bloc de donn??es
     *  @param x  Signal d'entr??e (r??f??rence vers un vecteur Eigen)
     *  @param y  Signal de sortie (ne peut ??tre une r??f??rence car taille variable)
     */
    virtual void step(const Eigen::Ref<const Vecteur<Te>> x, Vecteur<Ts> &y) = 0;


    /** @brief Traitement d'un bloc de donn??es
     *  @param x  Signal d'entr??e (r??f??rence vers un vecteur Eigen)
     *  @returns  Signal de sortie
     */
    Vecteur<Ts> step(const Eigen::Ref<const Vecteur<Te>> x)
    {
      Vecteur<Ts> y;
      step(x, y);
      return y;
    }

    Ts step(Te x)
    {
      Vecteur<Te> vx(1);
      Vecteur<Ts> vy(1);
      vx(0) = x;
      step(vx, vy);
      return vy(0);
    }
  };

  /** @brief Structure abstraite pour un filtre de donn??es configurable
   *  @tparam  Te Type d'entr??e
   *  @tparam  Ts Type de sortie
   *  @tparam  Tc Type de configuration
   */
  template<typename Te, typename Ts = Te, typename Tc = Void>
  struct Filtre : Configurable<Tc>, FiltreGen<Te, Ts>
  {
    virtual ~Filtre(){}
  };

  /** @cond undoc */
  namespace filtrage
  {
  template<typename T>
    sptr<FiltreGen<T>> filtre_reechan(float ratio);
  }
  /** @endcond */


  /** @brief      R??-??chantillonnage d'un signal suivant un facteur arbitraire.
   *
   *  <h3>R??-??chantillonnage</h3>
   *
   *  Cette fonction r??-??chantillonne un signal suivant un facteur de d??cimation (@f$r < 1@f$)
   *  ou d'interpolation (@f$r > 1@f$) arbitraire.
   *
   *  @param x    Signal d'entr??e
   *  @param r    Facteur de r??-??chantillonnage (@f$f_e^{(2)} / f_e^{(1)}@f$)
   *  @returns    Signal r??-??chantillonn??
   *
   *  @note Un filtre (voir une cascade de filtres) est automatiquement ins??r?? (avant la d??cimation ou apr??s le sur-??chantillonnage) de
   *  mani??re ?? ??viter le repliement de spectre.
   *
   *  @sa filtre_reechan() (pour traiter des signaux au fil de l'eau), reechan_freq() (r????chantillonnnage sans retard).
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_resample
   *  @image html resample.png width=800px
   *
   *  @note Du fait des filtres anti-repliement utilis??s, le signal d'entr??e est retard??.
   *  Pour r??-??chantilloner sans introduire de retard, voir aussi la fonction @ref r????chan_freq().
   */
  template<typename Derived>
    auto r????chan(const Eigen::ArrayBase<Derived> &x, float r)
  {
    auto f = tsd::filtrage::filtre_reechan<typename Derived::Scalar>(r);
    return f->step(x);
  }


  /** @brief A partir d'un flux de donn??es entrant, d??coup?? en paquet de dimensions quelconques,
   *  produit un flux de donn??es sous forme de <b>paquets de dimension fixe</b>.
   *
   *  <h3>Pr??paration d'un flux de paquets de donn??es de dimension fixes</h3>
   *
   *  A partir d'un flux de donn??es entrant, d??coup?? en paquet de dimensions quelconques,
   *  ce bloc produit un flux de donn??es sous forme de <b>paquets de dimension fixe</b>.
   *
   *  @image html tampon.png width=800px
   *
   *  @param N dimension des paquets de sortie
   *  @param callback Fonction utilisateur qui sera appel??e ?? chaque fois qu'un nouveau paquet de dimension @f$N@f$ est disponible.
   *  @returns Un puit de donn??e acceptant des vecteurs de type @p T.
   *  @tparam T Type de donn??es d'entr??e / sortie (e.g. float, cfloat, double, ...).
   *
   *
   *  @par Exemple : pr??paration des paquets pour une FFT
   *  @snippet exemples/src/ex-tsd.cc ex_tampon
   *  Notez que 10000 n'??tant pas un multiple de 512,
   *  les derniers ??chantillons ne sont pas trait??s (ils le seront si la m??thode t->step())
   *  est appel??e ?? nouveau avec d'autres ??chantillons).
   *
   *
   */
  template<typename T>
    sptr<Sink<T, int>> tampon_cr??ation(int N, std::function<void (const Vecteur<T> &)> callback);



  template<typename T>
  struct est_complex_t : public std::false_type {};

  template<typename T>
  struct est_complex_t<std::complex<T>> : public std::true_type {};

  template<typename T>
  constexpr bool est_complexe() { return est_complex_t<T>::value; }


  // En attendant que la fonction eigen src(index) soit dispo dans Eigen
  template<typename T>
  T sousvec(const T &src, const std::vector<int> &index)
  {
    T res(index.size());
    for(auto i = 0u; i < index.size(); i++)
    {
      assert(index[i] < src.rows());
      res(i) = src(index[i]);
    }
    return res;
  }

  /** @brief Modulo avec r??sultat compris dans l'intervalle @f$[0,m[@f$.
   *
   * <h3>Modulo avec r??sultat compris dans l'intervalle @f$[0,m[@f$</h3>
   *
   * Contrairement ?? la fonction standard <code>fmod(x, m)</code>, qui renvoie
   * une valeur comprise entre @f$-m@f$ et @f$m@f$, cette fonction
   * renvoie une valeur comprise entre @f$0@f$ et @f$m@f$ :
   *  @f[
   *  y = x + k\cdot m,\quad k\in\mathbb{Z},\ y\in \left[0,m\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  auto y = modulo(1.5, 1.0); // -> y = 0.5
   *  @endcode
   *
   *  @sa modulo_2??(), modulo_pm_2??()
   *
   */
  template<typename T>
  T modulo(T x, T m)
  {
    if(m == 0)
      return x;

    auto res = x - m * ((T) floor(x / m)); // floor(-0.1) = -1
    // Peut arriver, si x = -une valeur tr??s faible
    // Par exemple, x = -??
    // res = -?? + y, ce qui peut ??tre ??gal ?? y par troncature.
    if(res == m)
      res = 0;
    return res;
  }

  /** @brief Calcule @f$a@f$ modulo @f$2\pi@f$, r??sultat dans l'intervalle @f$\left[0,2\pi\right[@f$.
   *
   *
   *  <h3>Modulo @f$2\pi@f$</h3>
   *  Le r??sultat est dans l'intervalle @f$\left[0,2\pi\right[@f$ :
   *  @f[
   *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[0,2\pi\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  tsd_assert(modulo_2??(2*??+1e-5) == 1e-5);
   *  @endcode
   *
   *  @sa modulo_pm_??(), d??plie_phase(), modulo()
   *
   **/
  template<typename T>
  inline T modulo_2??(T x)
  {
    return modulo(x, 2 * std::numbers::pi_v<T>);
  }

  /** @brief Calcule @f$a@f$ modulo @f$2\pi@f$, r??sultat dans l'intervalle @f$\left[-\pi,\pi\right[@f$.
   *
   *  <h3>Modulo @f$2\pi@f$</h3>
   *  Le r??sultat est dans l'intervalle @f$\left[-\pi,\pi\right[@f$ :
   *  @f[
   *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[-\pi,\pi\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  tsd_assert(modulo_pm_??(2*??-1e-5) == -1e-5);
   *  @endcode
   *
   *  @sa modulo_2??(), d??plie_phase(), modulo()
   *
   **/
  template<typename T>
  T modulo_pm_??(T x)
  {
    auto ?? = std::numbers::pi_v<T>;
    return modulo_2??(x + ??) - ??;
  }

  /** @brief Conversion degr??s vers radians.
   *
   *  <h3>Conversion degr??s vers radians</h3>
   *
   *  @f[
   *  y = \frac{\pi\cdot x}{180}
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  tsd_assert(abs(deg2rad(45) - ??/4) < 1e-15);
   *  @endcode
   *
   *  @sa rad2deg()
   */
  template<typename T>
  auto deg2rad(T degrees)
  {
    if constexpr(std::is_integral_v<T>)
      return degrees * ?? / 180;
    else
      return degrees * std::numbers::pi_v<T> / 180.0;
  }

  /** @brief Conversion radians vers degr??s.
   *
   *  <h3>Conversion degr??s vers radians</h3>
   *
   *  @f[
   *  y = \frac{180\cdot x}{\pi}
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  tsd_assert(abs(rad2deg(??/4) - 45) < 1e-15);
   *  @endcode
   *
   *  @sa deg2rad()
   */
  template<typename T>
  T rad2deg(T radians)
  {
    return radians * 180.0 / std::numbers::pi_v<T>;
  }




  /** @} */

  /** @addtogroup tsd-gen
   *  @{ */

  /** @brief Intervalle de points ??quidistants.
   *
   *  <h3>Intervalle de points ??quidistants</h3>
   *
   *  Calcule @f$n@f$ points ??quidistants entre @f$a@f$ et @f$b@f$.
   *
   *  @param a point initial (premi??re valeur)
   *  @param b point final (derni??re valeur)
   *  @param n nombre de points
   *  @returns Intervalle d??fini par :
   *  @f[
   *  x_k = a + k \cdot \frac{b-a}{n-1},\quad k = 0,\dots,n-1
   *  @f]
   *
   *  @par Exemple 1 : vecteur des n premiers entiers
   *  @code
   *  ArrayXf t = linspace(0, n-1, n)
   *  // t = 0, 1, 2, ..., n-1
   *  @endcode
   *  @image html linspace.png width=600px
   *  @par Exemple 2 : vecteur de n valeurs de temps, ?? la fr??quence d'??chantillonnage de fe Hz :
   *  @code
   *  ArrayXf t = linspace(0, (n-1)/fe, n)
   *  @endcode
   *
   *  @sa logspace()
   */
  static inline auto linspace(float a, float b, int n)
  {
    return ArrayXf::LinSpaced(n, a, b);
  }


  /** @brief Intervalle de points logarithmiquement ??quidistants (suite g??om??trique).
   *
   *  <h3>Intervalle de points logarithmiquement ??quidistants</h3>
   *
   *  Calcule @f$n@f$ points logarithmiquement ??quidistants entre @f$a@f$ et @f$b@f$,
   *  c'est-??-dire la suite g??m??trique suivante :
   *
   *  @f[
   *  x_k = 10^{a + k  \frac{b-a}{n-1}},\quad k = 0,\dots,n-1
   *  @f]
   *
   *
   *  @param a Logarithme en base 10 du point initial
   *  @param b Logarithme en base 10 du point final
   *  @param n Nombre de points
   *  @returns Suite g??om??trique
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_logspace
   *  @image html logspace.png width=600px
   *
   *  @sa linspace()
   */
  static inline auto logspace(float a, float b, int n)
  {
    return (10 * ArrayXf::Ones(n)).pow(linspace(a, b, n));
  }

  /** @brief Intervalle entier
   *
   *  <h3>Intervalle entier</h3>
   *
   *  @return Vecteur des @f$b-a+1@f$ entiers entre @f$a@f$ et @f$b@f$ :
   *  @f[
   *   x_k = a + k,\quad k = 0,\dots,(b-a+1)
   *  @f]
   *
   *  @sa linspace(), trange()
   */
  static inline ArrayXi intervalle_entier(int a, int b)
  {
    return ArrayXi::LinSpaced((b - a) + 1, a, b);
  }

  /** @brief Intervalle temporel, en fonction de la fr??quence d'??chantillonnage
   *
   *  <h3>Intervalle temporel</h3>
   *
   *  @param n nombre de points
   *  @param fe fr??quence d'??chantillonnage
   *  @return Ensemble de valeurs temporelles ??quir??parties et ??chantillon??es ?? la fr??quence @f$f_e@f$ :
   *  @f[
   *   x_k = \frac{k}{f_e},\quad k = 0,\dots,n-1
   *  @f]
   *
   *  @sa linspace(), irange()
   */
  static inline ArrayXf intervalle_temporel(unsigned int n, float fe)
  {
    return linspace(0, (n-1) / fe, n);
  }

  /** @brief Loi normale (vecteur colonne).
   *
   *  <h3>Loi normale</h3>
   *  G??n??ration d'un vecteur (colonne) d'??chantillons d'une loi normale (@f$\mathcal{N}(0,1)@f$).
   *
   *  @param n nombre de points ?? g??n??rer
   *  @returns Un vecteur de valeurs al??atoires, ??chantillon??es suivant une loi normale (variance unitaire).
   *
   *  @note Notez que pour g??n??rer une loi Gaussienne plus g??n??rale (variance diff??rente de 1, et moyenne diff??rente de z??ro),
   *  il suffit d'effectuer une mise ?? l'??chelle.
   *  Par exemple, pour une moyenne de 5, et une variance de @f$1/2@f$ :
   *  @code
   *  ArrayXf x = 5 + 0.5 * randn(n);
   *  @endcode
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randn
   *  @image html randn.png width=600px
   *
   *  @sa randu(), randb()
   */
  extern ArrayXf randn(int n);

  /** @brief G??n??rateur al??atoire utilis?? en interne */
  extern std::default_random_engine generateur_aleatoire;

  /** @cond private */

  /** @brief Loi normale (tableau 2d)
   *  @sa randu() */
  extern ArrayXXf randn_2d(unsigned int n, unsigned int m);

  /** @brief Loi uniforme, intervalle [0,1] (tableau 2d) */
  extern ArrayXXf randu_2d(unsigned int n, unsigned int m);

  /** @endcond */

  /** @brief Loi uniforme (vecteur colonne).
   *
   *  <h3>Loi uniforme</h3>
   *
   *  @param n Nombre de points ?? g??n??rer
   *  @returns Un vecteur de valeurs al??atoires, ??chantillon??es suivant une loi uniforme entre 0 et 1.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randu
   *  @image html randu.png width=600px
   *
   *  @sa randn(), randi(), randb()
   */
  extern ArrayXf randu(int n);

  /** @brief Loi al??atoire binaire.
   *
   *  <h3>Loi binaire</h3>
   *
   *  @param n nombre de points ?? g??n??rer
   *  @returns Un vecteur de valeurs al??atoire de ??gales ?? 0 ou 1 uniquement.
   *
   *  Cette fonction peut ??tre utile pour g??n??rer un train de bits al??atoires.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randb
   *  @image html randb.png width=600px
   *
   *
   *  @sa randn(), randu(), randi()
   */
  extern ArrayXb randb(int n);

  /** @brief Loi al??atoire cat??gorielle
   *
   *  <h3>Loi cat??gorielle</h3>
   *
   *  @param M Nombre de cat??gories
   *  @param n Nombre de valeurs ?? g??n??rer
   *  @returns Un vecteur de @f$n@f$ entiers entre @f$0@f$ et @f$M-1@f$.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randi
   *  @image html randi.png width=600px
   *
   *  @sa randn(), randu(), randb()
   */
  extern ArrayXi randi(int M, int n);

  /** @brief Calcul efficace d'une exponentielle complexe
   *
   *  <h3>G??n??ration d'une exponentielle complexe</h3>
   *
   *  Cette fonction g??n??re un signal exponentiel, le calcul ??tant bas?? sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonom??trique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une d??rive, tant en gain qu'en phase.
   *
   *  @param f Fr??quence normalis??e (entre -0,5 et 0,5)
   *  @param n Nombre de points ?? g??n??rer
   *  @returns Vecteur complexe d??fini par :
   *  @f[
   *   y_k = e^{2\pi\mathbf{i} k f},\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigexp
   *  @image html sigexp.png width=600px
   *
   *  @sa sigcos(), sigsin(), sigcar(), sigtri()
   */
  extern ArrayXcf sigexp(float f, int n);

  /** @brief Calcul efficace d'une sinusoide.
   *
   *  <h3>G??n??ration d'une sinusoide</h3>
   *
   *  Cette fonction g??n??re un signal sinusoidal, le calcul ??tant bas?? sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonom??trique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une d??rive, tant en gain qu'en phase.
   *
   *  @param f Fr??quence normalis??e (entre -0,5 et 0,5)
   *  @param n Nombre de points ?? g??n??rer
   *  @returns Vecteur r??el d??fini par :
   *  @f[
   *   y_k = \sin\left(2\pi k f\right),\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigsin
   *  @image html sigsin.png width=600px
   *
   *  @sa sigcos(), sigexp(), sigcar(), sigtri()
   */
  extern ArrayXf sigsin(float f, int n);

  /** @brief Calcul efficace d'un cosinus
   *
   *  <h3>G??n??ration d'un cosinus</h3>
   *
   *  Cette fonction g??n??re un signal sinusoidal, le calcul ??tant bas?? sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonom??trique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une d??rive, tant en gain qu'en phase.
   *
   *  @param f Fr??quence normalis??e (entre -0,5 et 0,5)
   *  @param n Nombre de points ?? g??n??rer
   *  @returns Vecteur r??el d??fini par :
   *  @f[
   *   y_k = \cos\left(2\pi k f\right),\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigcos
   *  @image html sigcos.png width=600px
   *
   *  @sa sigsin(), sigexp()
   */
  extern ArrayXf sigcos(float f, int n);


  /** @brief Signal triangulaire (p??riodique).
   *
   *  <h3>Signal triangulaire (p??riodique)</h3>
   *
   *  @param p Demi-p??riode, en nombre d'??chantillons.
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigtri
   *  @image html sigtri.png width=600px
   *
   *  @sa sigcar(), sigsin(), sigcos(), sigexp()
   */
  extern ArrayXf sigtri(int p, int n);

  /** @brief Signal carr?? (p??riodique)
   *
   *  <h3>Signal carr?? (p??riodique)</h3>
   *
   *  @param p Demi-p??riode, en nombre d'??chantillons.
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigcar
   *  @image html sigcar.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern ArrayXf sigcar(int p, int n);

  /** @brief Impulsion discr??te
   *
   *  <h3>Impulsion discr??te</h3>
   *
   *  @param n Dimension du vecteur
   *  @param p Position de l'impulsion
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigimp
   *  @image html sigimp.png width=600px
   *
   */
  extern ArrayXf sigimp(int n, int p = 0);


  /** @brief Signal en dent de scie
   *
   *  <h3>Signal en dent de scie</h3>
   *
   *  @param p P??riode, en nombre d'??chantillons.
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigscie
   *  @image html sigscie.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern ArrayXf sigscie(int p, int n);

  /** @brief Sinuso??de modul??e par une Gaussienne
   *
   *  <h3>Sinuso??de modul??e par une Gaussienne</h3>
   *
   *  @f[
   *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2} \cdot \sin \left(2\pi fk\right)
   *  @f]
   *
   *
   *  @param f Fr??quence (normalis??e).
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *  @param a Param??tre de d??croissance
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_siggsin
   *  @image html siggsin.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern ArrayXf siggsin(float f, int n, float a = 10);


  /** @brief Impulsion filtr??e par une Gaussienne
     *
     *  <h3>Impulsion filtr??e par une Gaussienne</h3>
     *
     *  @f[
     *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2}
     *  @f]
     *
     *
     *  @param n Nombre d'??chantillons ?? g??n??rer.
     *  @param a Param??tre de d??croissance
     *
     *  @par Exemple
     *  @snippet exemples/src/ex-tsd.cc ex_siggaus
     *  @image html sigsiggaus.png width=600px
     *
     *  @sa sigtri(), sigsin(), sigcos(), sigexp(), siggsin()
     */
    extern ArrayXf siggauss(int n, float a = 10);

  /** @brief Chirp lin??aire
   *
   *  <h3>Chirp lin??aire</h3>
   *
   *  @f[
   *  x_k = \cos \phi_k, \quad \phi_k = 2 \pi \sum_{i=0}^k f_k
   *  @f]
   *
   *  Les @f$f_k@f$ ??tant r??partis lin??airement entre @f$f_0@f$ et @f$f_1@f$.
   *
   *  @param f0 Fr??quence initiale (normalis??e).
   *  @param f1 Fr??quence finale (normalis??e).
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigchirp
   *  @image html sigchirp.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp(), sigchirp2()
   */
  extern ArrayXf sigchirp(float f0, float f1, int n);

  /** @brief Chirp quadratique
   *
   *  <h3>Chirp quadratique</h3>
   *
   *  @f[
   *  x_k = \cos \phi_k, \quad \phi_k = 2 \pi \sum_{i=0}^k f_k
   *  @f]
   *
   *  Les @f$f_k@f$ ??tant r??partis quadratiquement entre @f$f_0@f$ et @f$f_1@f$ :
   *  @f[
   *  f_k = f_0 + (f_1 - f_0) \left(\frac{k}{n-1}\right)^2
   *  @f]
   *
   *  @param f0 Fr??quence initiale (normalis??e).
   *  @param f1 Fr??quence finale (normalis??e).
   *  @param n Nombre d'??chantillons ?? g??n??rer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigchirp2
   *  @image html sigchirp2.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp(), sigchirp()
   */
  extern ArrayXf sigchirp2(float f0, float f1, int n);


  struct OHConfig
  {
    float freq = 0;

    // Si vrai, applique un shift de la fr??quence indiqu??e
    bool shift = false;
    float df = 0;
  };

  /** @brief G??n??ration d'un signal exponentiel via un oscillateur harmonique
   *
   *  <h3>Oscillateur harmonique (sortie complexe)</h3>
   *
   *  Cette fonction renvoie une source de donn??es, qui peut ??tre appell??e plusieurs fois (g??n??ration d'un flux continu d'??chantillons, ?? la diff??rence de @ref sigexp(), qui ne peut g??n??rer qu'un nombre fini et pr??determin?? d'??chantillons).
   *
   *  L'oscillateur en quadrature est bas?? sur un simple filtre r??cursif (complexe) d'ordre 1 :
   *  @f[
   *  z_k = z_{k-1} \cdot e^{2\pi\mathbf{i}f}
   *  @f]
   *  qui permet donc de g??n??rer un signal exponentiel :
   *  @f[
   *  z_k = e^{2\pi\mathbf{i}kf}
   *  @f]
   *
   *  @note Pour g??n??rer une exponentielle de <b>fr??quence fixe</b>, cette m??thode sera beaucoup plus efficace
   *  que l'utilisation des fonctions trigonom??triques de la librairie standard.
   *
   *  @param freq Fr??quence normalis??e (entre -0,5 et 0,5).
   *  @return Une source de donn??e complexe.
   *
   *  @par Exemple
   *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohc
   *  @image html ohc.png width=600px
   *
   *  @sa source_ohr(), sigexp()
   */
  extern sptr<Source<cfloat, OHConfig>> source_ohc(float freq);

  /** @brief G??n??ration d'un signal sinuso??dal via un oscillateur harmonique
   *
   *  <h3>Oscillateur harmonique (sortie r??elle)</h3>
   *
   *  Cette fonction renvoie une source de donn??es, qui peut ??tre appell??e plusieurs fois (g??n??ration d'un flux continu d'??chantillons, ?? la diff??rence de @ref sigcos(), qui ne peut g??n??rer qu'un nombre fini et pr??determin?? d'??chantillons).
   *
   *  L'oscillateur est bas?? sur un oscillateur en quadrature standard (voir @ref source_ohc()), dont
   *  on n'utilise que la partie r??elle (cosinus).
   *
   *
   *  @param freq Fr??quence normalis??e (entre -0,5 et 0,5).
   *  @return Une source de donn??e r??elle.
   *
   *  @par Exemple
   *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohr
   *  @image html ohr.png width=600px
   *
   *  @sa source_ohc(), sigcos()
   */
  extern sptr<Source<float, OHConfig>> source_ohr(float freq);

  struct OLUT
  {
    OLUT(int resolution = 1024);
    cfloat step(float ??);
  private:
    struct Impl;
    sptr<Impl> impl;
  };

  /** @} */



}






