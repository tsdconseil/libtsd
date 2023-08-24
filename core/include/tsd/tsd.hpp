#ifndef TSD_TSD_H
#define TSD_TSD_H

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include <numbers>
#include <random>
#include <algorithm>

#include "tsd/tableau.hpp"




namespace cutils::model
{
  class Node;
}


using std::sin;
using std::cos;
using std::tan;
using std::log;
using std::log2;
using std::log10;
using std::exp;
using std::abs;
using std::sqrt;
using std::pow;
using std::polar;
using std::max;
using std::min;
using std::clamp;
using std::real;
using std::imag;
using std::norm;
using std::round;
using std::ceil;
using std::floor;

namespace tsd
{

using std::sin;
using std::cos;
using std::tan;
using std::log;
using std::log10;
using std::exp;
using std::abs;
using std::sqrt;
using std::pow;
using std::polar;
using std::max;
using std::min;
using std::real;
using std::imag;
using std::norm;

/** @addtogroup tsd
 *  @{ */


  /** @brief Valeur de @f$\pi@f$, en précision flottante 64 bits. */
  static const double π = std::numbers::pi;

  /** @brief Valeur de @f$\pi@f$, en précision flottante 32 bits. */
  static const float π_f = std::numbers::pi_v<float>;


  using namespace std::complex_literals;


  /** @brief Retourne -1, 1 ou 0 (en cas d'égalité). */
  template <typename T> entier signe(T val) {
      retourne (T(0) < val) - (val < T(0));
  }

  /** @cond private */

  using fmt::format;



  /** @endcond */

  template<typename t>
    using Tableau = TabT<t, 2>;

  /** @brief Raccourcis pour un vecteur colonne de données de type générique (paramètre de template t) */
  template<typename t>
    using Vecteur = VecT<t>;



  /** @brief Concaténation verticale de deux vecteurs.
   *
   * Cette fonction réalise l'assemblage de deux vecteurs colonnes :
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
   * (en supposant que a et b aient respectivement @f$n@f$ et @f$m@f$ éléments).
   * L'opérateur | a été redéfini pour réaliser la même opération.
   *
   * @param a Premier vecteur
   * @param b Deuxième vecteur
   * @returns Concaténation verticale
   *
   * @par Exemple
   * @code
   *  soit a = Vecf::ones(1),
   *       b = Vecf::zeros(5),
   *       c = vconcat(a, b); // (ou c = a | b)
   * @endcode
   *
   */
  template<typename T>
    auto vconcat(const Vecteur<T> &a, const Vecteur<T> &b)
  {
    Vecteur<T> c(a.rows()+b.rows());
    si(a.rows() > 0)
      c.head(a.rows()) = a;
    si(b.rows() > 0)
      c.tail(b.rows()) = b;
    retourne c;
  }

  /** @brief Concaténation verticale de deux vecteurs (raccourci).
   *  @sa vconcat() */
  template<typename T>
    auto operator |(const Vecteur<T> &a, const Vecteur<T> &b)
  {
    return vconcat(a, b);
  }

  /** @brief Concaténation verticale de deux tableaux. */
  template<typename T>
    auto vconcat(const Tableau<T> &a, const Tableau<T> &b)
  {
    assertion_msg(a.cols() == b.cols(), "Concaténation verticale : le nombre de colonnes devrait être identique.");
    Tableau<T> c(a.cols(), a.rows()+b.rows());
    si(a.rows() > 0)
      c.block(0, 0, a.rows(), a.cols()) = a;
    si(b.rows() > 0)
      c.block(a.rows(), 0, b.rows(), a.cols()) = b;
    retourne c;
  }

  /** @brief Concaténation horizontale de deux tableaux. */
  template<typename T>
    auto hconcat(const Tableau<T> &a, const Tableau<T> &b)
  {
    assertion_msg(a.rows() == b.rows(), "Concaténation horizontale : le nombre de lignes devrait être identique.");
    Tableau<T> c(a.cols() + b.cols(), a.rows());
    si(a.cols() > 0)
      c.block(0, 0, a.rows(), a.cols()) = a;
    si(b.cols() > 0)
      c.block(0, a.cols(), b.rows(), b.cols()) = b;
    retourne c;
  }

  /** @cond undoc */
  /* Calcule la longueur d'un vecteur colonne */
  template<typename T>
    auto length(const Vecteur<T> &x)
  {
    retourne x.rows();
  }
  /** @endcond */

  /** @brief Rotation d'un vecteur
   *
   *  @param x Vecteur colonne d'entrée
   *  @param d Nombre entier de pas (positif ou négatif)
   *  @returns Vecteur des données décalée de @f$d@f$ pas (modulo la dimension du vecteur) :
   *  @f[
   *   y_k = x_{k + d\ [N]}
   *  @f]
   *  @f$N@f$ étant la dimension du vecteur.
   *
   *  @par Exemple
   *  @code
   *  soit x    = linspace(0, 5, 6),
   *       y    = rotation_vec(x, 2),
   *       yref = Vecf::valeurs({2, 3, 4, 5, 0, 1});
   *  assertion(y.isApprox(yref));
   *  @endcode
   */
  template<typename T>
  Vecteur<T> rotation_vec(const Vecteur<T> &x, entier d)
  {
    soit n = x.dim();
    retourne Vecteur<T>::int_expr(n,
        IMAP(x((i + d + n) % n)));
  }

  /** @brief Différence entre 2 élements successifs d'un vecteur
   *
   *  @returns Un vecteur de longueur @f$n-1@f$ :
   *  @f[
   *  y_k = x_{k+1} - x_k,\quad k = 0,\dots, n-2
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  soit x    = randn(10),
   *       y    = diff(x),
   *       yref = x.tail(9) - x.head(9);
   *
   *  assertion(y.isApprox(yref));
   *  @endcode
   *
   *  @sa cumsum()
   */
  template<typename T>
    auto diff(const Vecteur<T> &x)
  {
    soit n = x.rows();
    retourne x.tail(n-1) - x.head(n-1);
  }

  /** @brief Somme cumulée d'un vecteur
   *
   *  Calcule un vecteur de dimension identique et égale à l'acculation des
   *  échantillons du vecteurs d'entrée :
   *  @f[
   *  y_n = \sum_{k=0}^n x_k,\quad k = 0\dots N-1
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  soit x     = linspace(0, 99, 100),
   *       y     = cumsum(x),
   *       yref  = x * (x + 1) / 2;
   *  assertion(y.isApprox(yref));
   *  @endcode
   *
   *  @sa diff()
   *
   */
  template<typename T>
      auto cumsum(const Vecteur<T> &x)
  {
    soit n = x.rows();
    Vecteur<T> y(n);
    T cs = 0;
    Pour(auto i = 0; i < n; i++)
    {
      cs += x(i);
      y(i) = cs;
    }
    retourne y;
  }

  /** @brief Corrige les sauts de phases
   *
   *  @param x Vecteur contenant typiquement des angles
   *  @param r Définit la classe d'équivalence sur les valeurs de @f$x@f$
   *  @returns Un vecteur @f$y@f$ tel que @f$y_k = x_k + k \cdot r @f$, et avec le moins de discontinuités possibles (plus précisément, aucune discontinuité supérieure à @f$r/2@f$ en valeur absolue).
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc exemple_unwrap
   *  @image html unwrap.png
   *
   *  @sa modulo_pm_π(), modulo_2π()
   *
   */
  template<typename T>
    Vecteur<T> déplie_phase(const Vecteur<T> &x, float r = 2*π);

  /** @brief Recherche d'éléments vrais dans un tableau bouléen.
   *
   *  @param x Tableau de valeurs bouléennes
   *  @returns Vecteur des index pour lequels les éléments du tableau sont vrais.
   *
   *  @par Exemple
   *  @code
   *  soit x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
   *  soit idx = trouve(x >= 0); // idx = {2, 3, 4}
   *  @endcode
   *
   *  @sa trouve_premier(), trouve_dernier()
   */
  extern vector<entier> trouve(const Vecb &x);

  /** @brief Recherche de l'indice du premier élément vrai dans un tableau de bouléens.
   *
   *
   *  @par Exemple
   *  @code
   *  soit x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
   *  soit idx = trouve_premier(x >= 0); // idx = 2
   *  @endcode
   *
   *  @sa trouve(), trouve_dernier()
   */
  extern entier trouve_premier(const Vecb &x);

  /** @brief Recherche de l'indice du dernier élément vrai dans un tableau de bouléens.
   *  @sa trouve_premier(), trouve() */
  extern entier trouve_dernier(const Vecb &x);


  /** @brief Recherche du premier maximum local (valeur supérieure à ces deux voisins). */
  template<typename T>
    entier trouve_premier_max_local(const Vecteur<T> &x)
  {
    pour(auto i = 1; i + 1 < x.rows(); i++)
      si((x(i) > x(i-1)) && (x(i) > x(i+1)))
        retourne i;
    retourne -1;
  }


  template<typename T>
   auto polar(const Vecteur<T> &θ)
  {
    retourne Vecteur<std::complex<T>>::int_expr(θ.dim(),
        IMAP(std::polar((T) 1, θ(i))));
  }

  template<typename T>
    auto polar(const Vecteur<T> &ρ, const Vecteur<T> &θ)
  {
    assertion(ρ.dim() == θ.dim());
    retourne Vecteur<std::complex<T>>::int_expr(θ.dim(),
        IMAP(std::polar(ρ(i), θ(i))));
  }

  // Equivalent de : y = x(1:pas:$)
  /** @brief Sous-échantillonnage d'un vecteur colonne
   *
   *  A partir d'un vecteur @f$(x_k), k=0\dots n-1@f$, renvoie
   *  un vecteur sous-échantillonné d'un facteur @f$R@f$ :
   *  @f[
   *  y_k = x_{kR},\quad k =0 \dots \left\lfloor N/R\right\rfloor -1
   *  @f]
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc ex_sousech
   *  @image html sousech.png
   *
   *  @sa surech()
   *
   */
  template<typename T>
   auto sousech(const Vecteur<T> &x, entier R)
   {
    soit n = x.dim();
    Vecteur<T> y(n / R);
    pour(auto i = 0; i < n / R; i++)
      y(i) = x(i * R);
    retourne y;
   }

  /** @brief Sur-échantillonnage d'un vecteur colonne
   *
   *  A partir d'un vecteur @f$(x_k), k=0\dots n-1@f$, renvoie
   *  un vecteur sur-échantillonné d'un facteur @f$R@f$, par l'insertion de zéros :
   *  @f[
   *  y_k = \begin{cases} x_{k/R} & \textrm{si } k\textrm{ est multiple de } R,\\
   *   0 & \textrm{sinon.}\end{cases}
   *  @f]
   *
   *  @par Exemple :
   *  @snippet exemples/src/ex-tsd.cc ex_surech
   *  @image html surech.png
   *
   *  @sa sousech()
   *
   */
  template<typename T>
      auto surech(const Vecteur<T> &x, entier R)
  {
    entier n = x.dim();
    soit y = Vecteur<T>::zeros(n * R);
    pour(auto i = 0; i < n; i++)
      y(i * R) = x(i);
    retourne y;
  }


  /** @brief Conversion linéaire vers decibels
   *
   *  @param x Valeur en linéaire
   *  @returns Valeur en dB :
   *  @f[
   *    y = 10 \log_{10}(x)
   *  @f]
   *
   *  @sa db2pow(), mag2db()
   */
  template<typename T>
    auto pow2db(const T &x)
  {
    si constexpr (std::is_integral_v<T>)
      retourne 10 * log10(x * 1.0);
    sinon
      retourne 10 * log10(x);
  }

  /** @brief Conversion magnitude vers decibels
   *
   *  @param x Valeur en linéaire
   *  @returns Valeur en dB :
   *  @f[
   *    y = 20 \log_{10}(x)
   *  @f]
   *
   *  @sa pow2db(), db2pow()
   */
  auto mag2db(const auto &x)
  {
    retourne 2 * pow2db(x);
  }

  /** @brief Conversion decibels vers linéaire (puissance)
   *
   *  @param x Valeur en dB
   *  @returns Valeur en linéaire :
   *  @f[
   *    y = 10^{x/10}
   *  @f]
   *
   *  @sa pow2db(), db2mag(), mag2db()
   */
  template<typename T>
    auto db2pow(const T &x)
  {
    si constexpr (std::is_integral_v<T>)
      retourne pow(10.0, x/10.0);
    sinon
      retourne pow((T) 10, x/10);
  }

  /** @brief Conversion décibels vers linéaire (magnitude)
   *
   *  @param x Valeur en dB
   *  @returns Valeur en linéaire :
   *  @f[
   *    y = 10^{x/20}
   *  @f]
   *
   *  @sa mag2db(), pow2db(), db2pow()
   */
  auto db2mag(const auto &x)
  {
    retourne db2pow(x/2);
  }


  /** @brief Retourne la plus petite puissance de 2 supérieure ou égale à i.
   *
   * Calcule la plus petite puissance de 2 supérieure ou égale à @f$i@f$ :
   * @f[
   *  y = \min_{j\geq i,\ j = 2^k} j
   * @f]
   *
   * @par Exemple
   * @code
   * assertion(prochaine_puissance_de_2(3) == 4);
   * @endcode
   *
   */
  extern entier prochaine_puissance_de_2(entier i);




  /** @brief Ajoute des zéros à la fin d'un des deux vecteurs de manière à ce qu'ils aient la même dimension.
   *
   *  Un des deux vecteurs est complété avec des zéros de manière à ce que les deux vecteurs
   *  résultants aient la même longueur.
   *
   *  @param x  Premier vecteur
   *  @param y  Deuxième vecteur
   *  @param p2 Si vrai, arrondi la dimension des vecteurs de sortie à la prochaine puissance de 2.
   *  @return Un tuple de deux vecteurs de mêmes dimensions, et commençant respectivement par x et y (des zéros sont ajoutés à la fin).
   *
   *  @par Exemple
   *  @code
   *  soit x = linspace(0,4,5);
   *  // x = {0, 1, 2, 3, 4};
   *  soit y = linspace(0,3,4);
   *  // y = {0, 1, 2, 3};
   *
   *  soit [x2,y2] = pad_zeros(x,y);
   *  // x2 = {0, 1, 2, 3, 4};
   *  // y2 = {0, 1, 3, 4, 0};
   *  @endcode
   *
   */
  template<typename T>
    tuple<Vecteur<T>, Vecteur<T>> pad_zeros(
        const Vecteur<T> &x,
        const Vecteur<T> &y,
        bouléen p2 = non)
  {
    Vecteur<T> b1p, b2p;

    soit n1 = x.dim(), n2 = y.dim();
    soit n3 = max(n1, n2);

    si(p2)
      n3 = prochaine_puissance_de_2(n3);

    b1p = Vecteur<T>::zeros(n3);
    b2p = Vecteur<T>::zeros(n3);
    b1p.head(n1) = x;
    b2p.head(n2) = y;

    retourne {b1p, b2p};
  }



  /** @brief Calcule le carré d'un scalaire */
  template<typename T>
    T carré(const T &v) {retourne v * v;}

  /** @brief Structure abstraite pour une classe configurable.
   *  @tparam C Type de configuration */
  template<typename C>
  struct Configurable
  {
    virtual ~Configurable(){}

    /** @brief Configuration de la structure
     *  @param c Valeur de la configuration à appliquer
     *  @returns entier non nul en cas d'erreur. */
    void configure(const C &c)
    {
      if(callback_modif)
        callback_modif(c);
      configure_impl(c);
      config = c;
    }

    // retourne 0 si supporté et rend "configure_impl" inutile
    virtual entier configure(const cutils::model::Node &n)
    {
      //config.deserialise(n);
      //configure_impl(config);
      retourne -1;
    }

    /** @brief Méthode abstraite à implémenter par la classe dérivée */
    virtual void configure_impl(const C &c) = 0;

    /** @brief Lecture de la configuration actuelle. */
    const C &lis_config() const {return config;}

    /** @brief Callback optionnelle qui sera appelée à chaque changement de la configuration. */
    fonction<void (const C &nv_config)> callback_modif;

  protected:
    C config;
  };

  struct Void{};

  /** @brief Structure abstraite pour un puit de données générique */
  template<typename Te>
  struct SinkGen
  {
    virtual ~SinkGen(){}

    /** @brief Traitement d'un bloc de données
     *  @param x Signal d'entrée (référence vers un vecteur Eigen) */
    virtual void step(const Vecteur<Te> &x) = 0;
  };

  /** @brief Structure abstraite pour un puit de données générique et configurable */
  template<typename Te, typename Tc = Void>
  struct Sink : Configurable<Tc>, SinkGen<Te>
  {
    virtual ~Sink(){}
  };

  /** @brief Structure abstraite pour une source de données
   *  @tparam Ts Type de données à produire (float, cfloat, ...)
   *  @tparam C  Type de configuration */
  template<typename Ts>
  struct SourceGen
  {
    virtual ~SourceGen(){}
    /** @brief Production de n échantillons
     *  @param n Nombre d'échantillons à produire
     *  @returns Vecteur de n échantillons s*/
    virtual Vecteur<Ts> step(entier n) = 0;
  };

  /** @brief Structure abstraite pour source de données, configurable.
   *  @tparam Ts Type de données à produire (float, cfloat, ...)
   *  @tparam C  Type de configuration */
  template<typename Ts, typename C = Void>
  struct Source : SourceGen<Ts>, Configurable<C>
  {
    virtual ~Source(){}
  };

  /** @brief Structure abstraite pour un filtre de données
   *  @tparam Te Input type
   *  @tparam Ts Output type */
  template<typename Te, typename Ts = Te>
  struct FiltreGen
  {
    virtual ~FiltreGen(){}

    /** @brief Traitement d'un bloc de données
     *  @param x  Signal d'entrée (référence vers un vecteur Eigen)
     *  @param y  Signal de sortie (ne peut être une référence car taille variable)
     */
    virtual void step(const Vecteur<Te> &x, Vecteur<Ts> &y) = 0;


    /** @brief Traitement d'un bloc de données
     *  @param x  Signal d'entrée (référence vers un vecteur Eigen)
     *  @returns  Signal de sortie
     */
    Vecteur<Ts> step(const Vecteur<Te> &x)
    {
      Vecteur<Ts> y;
      step(x, y);
      retourne y;
    }

    Ts step(Te x)
    {
      Vecteur<Te> vx(1);
      Vecteur<Ts> vy(1);
      vx(0) = x;
      step(vx, vy);
      retourne vy(0);
    }
  };

  /** @brief Structure abstraite pour un filtre de données configurable
   *  @tparam  Te Type d'entrée
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
    sptr<Filtre<T,T,float>> filtre_reechan(float ratio);
  }
  /** @endcond */


  /** @brief      Ré-échantillonnage d'un signal suivant un facteur arbitraire.
   *
   *  Cette fonction ré-échantillonne un signal suivant un facteur de décimation (@f$r < 1@f$)
   *  ou d'interpolation (@f$r > 1@f$) arbitraire.
   *
   *  @param x    Signal d'entrée
   *  @param r    Facteur de ré-échantillonnage (@f$f_e^{(2)} / f_e^{(1)}@f$)
   *  @returns    Signal ré-échantillonné
   *
   *  @note Un filtre (voir une cascade de filtres) est automatiquement inséré (avant la décimation ou après le sur-échantillonnage) de
   *  manière à éviter le repliement de spectre.
   *
   *  @sa tsd::filtrage::filtre_reechan() (pour traiter des signaux au fil de l'eau), tsd::fourier::rééchan_freq() (rééchantillonnnage sans retard).
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_resample
   *  @image html resample.png width=800px
   *
   *  @note Du fait des filtres anti-repliement utilisés, le signal d'entrée est retardé.
   *  Pour ré-échantilloner sans introduire de retard, voir aussi la fonction @ref rééchan_freq().
   */
  template<typename T>
    auto rééchan(const Vecteur<T> &x, float r)
  {
    soit f = tsd::filtrage::filtre_reechan<T>(r);
    retourne f->step(x);
  }


  /** @brief A partir d'un flux de données entrant, découpé en paquet de dimensions quelconques,
   *  produit un flux de données sous forme de <b>paquets de dimension fixe</b>.
   *
   *  A partir d'un flux de données entrant, découpé en paquet de dimensions quelconques,
   *  ce bloc produit un flux de données sous forme de <b>paquets de dimension fixe</b>.
   *
   *  @image html tampon.png
   *
   *  @param N dimension des paquets de sortie
   *  @param callback Fonction utilisateur qui sera appelée à chaque fois qu'un nouveau paquet de dimension @f$N@f$ est disponible.
   *  @returns Un puit de donnée acceptant des vecteurs de type @p T.
   *  @tparam T Type de données d'entrée / sortie (e.g. float, cfloat, double, ...).
   *
   *
   *  @par Exemple : préparation des paquets pour une FFT
   *  @snippet exemples/src/ex-tsd.cc ex_tampon
   *  Notez que 10000 n'étant pas un multiple de 512,
   *  les derniers échantillons ne sont pas traités (ils le seront si la méthode t->step())
   *  est appelée à nouveau avec d'autres échantillons).
   *
   *
   */
  template<typename T>
    sptr<Sink<T, entier>> tampon_création(entier N, fonction<void (const Vecteur<T> &)> callback);






  // En attendant que la fonction eigen src(index) soit dispo dans Eigen
  template<typename T>
  T sousvec(const T &src, const vector<entier> &index)
  {
    soit n = (entier) index.size();
    T res(n);
    pour(auto i = 0; i < n; i++)
    {
      soit idx = index[i];
      assertion((idx >= 0) && (idx < src.rows()));
      res(i) = src(idx);
    }
    retourne res;
  }

  /** @brief Modulo avec résultat compris dans l'intervalle @f$[0,m[@f$.
   *
   * Contrairement à la fonction standard <code>fmod(x, m)</code>, qui renvoie
   * une valeur comprise entre @f$-m@f$ et @f$m@f$, cette fonction
   * renvoie une valeur comprise entre @f$0@f$ et @f$m@f$ :
   *  @f[
   *  y = x + k\cdot m,\quad k\in\mathbb{Z},\ y\in \left[0,m\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  soit y = modulo(1.5, 1.0); // -> y = 0.5
   *  @endcode
   *
   *  @sa modulo_2π(), modulo_pm_2π()
   *
   */
  template<typename T>
  T modulo(T x, T m)
  {
    si(m == 0)
      retourne x;

    soit res = x - m * ((T) floor(x / m)); // floor(-0.1) = -1
    // Peut arriver, si x = -une valeur très faible
    // Par exemple, x = -ε
    // res = -ε + y, ce qui peut être égal à y par troncature.
    si(res == m)
      res = 0;
    retourne res;
  }

  /** @brief Calcule @f$a@f$ modulo @f$2\pi@f$, résultat dans l'intervalle @f$\left[0,2\pi\right[@f$.
   *
   *  Le résultat est dans l'intervalle @f$\left[0,2\pi\right[@f$ :
   *  @f[
   *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[0,2\pi\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  assertion(modulo_2π(2*π+1e-5) == 1e-5);
   *  @endcode
   *
   *  @sa modulo_pm_π(), déplie_phase(), modulo()
   *
   **/
  template<typename T>
  inline T modulo_2π(T x)
  {
    retourne modulo(x, 2 * std::numbers::pi_v<T>);
  }

  /** @brief Calcule @f$a@f$ modulo @f$2\pi@f$, résultat dans l'intervalle @f$\left[-\pi,\pi\right[@f$.
   *
   *  Le résultat est dans l'intervalle @f$\left[-\pi,\pi\right[@f$ :
   *  @f[
   *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[-\pi,\pi\right[
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  assertion(modulo_pm_π(2*π-1e-5) == -1e-5);
   *  @endcode
   *
   *  @sa modulo_2π(), déplie_phase(), modulo()
   *
   **/
  template<typename T>
  T modulo_pm_π(T x)
  {
    soit π = std::numbers::pi_v<T>;
    retourne modulo_2π(x + π) - π;
  }

  template<typename T>
    bouléen est_pair(T t) NOECLIPSE(requires(std::is_integral_v<T>))
  {
    retourne ((t & 1) == 0);
  }

  template<typename T>
    bouléen est_impair(T t) NOECLIPSE(requires(std::is_integral_v<T>))
  {
    retourne ((t & 1) == 1);
  }


  /** @brief Conversion degrés vers radians.
   *
   *  @f[
   *  y = \frac{\pi\cdot x}{180}
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  assertion(abs(deg2rad(45) - π/4) < 1e-15);
   *  @endcode
   *
   *  @sa rad2deg()
   */
  template<typename T>
  auto deg2rad(T degrees)
  {
    si constexpr(std::is_integral_v<T>)
      retourne degrees * π / 180;
    sinon
      retourne degrees * std::numbers::pi_v<T> / 180;
  }

  /** @brief Conversion radians vers degrés.
   *
   *  @f[
   *  y = \frac{180\cdot x}{\pi}
   *  @f]
   *
   *  @par Exemple
   *  @code
   *  assertion(abs(rad2deg(π/4) - 45) < 1e-15);
   *  @endcode
   *
   *  @sa deg2rad()
   */
  template<typename T>
  T rad2deg(T radians)
  {
    retourne radians * 180.0 / std::numbers::pi_v<T>;
  }




  /** @} */

  /** @addtogroup tsd-gen
   *  @{ */

  /** @brief Intervalle de points équidistants.
   *
   *  Calcule @f$n@f$ points équidistants entre @f$a@f$ et @f$b@f$.
   *
   *  @param a point initial (première valeur)
   *  @param b point final (dernière valeur)
   *  @param n nombre de points
   *  @returns Intervalle défini par :
   *  @f[
   *  x_k = a + k \cdot \frac{b-a}{n-1},\quad k = 0,\dots,n-1
   *  @f]
   *
   *  @par Exemple 1 : vecteur des n premiers entiers
   *  @code
   *  soit t = linspace(0, n-1, n)
   *  // t = 0, 1, 2, ..., n-1
   *  @endcode
   *  @image html linspace.png
   *
   *  @par Exemple 2 : vecteur de n valeurs de temps, à la fréquence d'échantillonnage de fe Hz :
   *  @code
   *  soit t = linspace(0, (n-1)/fe, n)
   *  @endcode
   *
   *  @sa logspace(), intervalle_temporel()
   */
  static inline auto linspace(float a, float b, entier n)
  {
    Vecf x(n);

    si(n > 0)
      x(0) = a;

    si(n > 1)
    {
      soit step = ((double) b-a) / (n-1);
      pour(auto i = 1; i < n; i++)
        x(i) = a + step * i;
    }

    retourne x;
  }


  /** @brief Intervalle de points logarithmiquement équidistants (suite géométrique).
   *
   *  Calcule @f$n@f$ points logarithmiquement équidistants entre @f$a@f$ et @f$b@f$,
   *  c'est-à-dire la suite gémétrique suivante :
   *
   *  @f[
   *  x_k = 10^{a + k  \frac{b-a}{n-1}},\quad k = 0,\dots,n-1
   *  @f]
   *
   *
   *  @param a Logarithme en base 10 du point initial
   *  @param b Logarithme en base 10 du point final
   *  @param n Nombre de points
   *  @returns Suite géométrique
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_logspace
   *  @image html logspace.png
   *
   *  @sa linspace()
   */
  static inline auto logspace(float a, float b, entier n)
  {
    retourne pow(Vecf::ones(n) * 10.0f, linspace(a, b, n));
  }

  /** @brief Intervalle entier
   *
   *  @return Vecteur des @f$b-a+1@f$ entiers entre @f$a@f$ et @f$b@f$ :
   *  @f[
   *   x_k = a + k,\quad k = 0,\dots,(b-a+1)
   *  @f]
   *
   *  @sa linspace(), intervalle_temporel()
   */
  static inline Veci intervalle_entier(entier a, entier b)
  {
    Veci res((b - a) + 1);
    soit j = 0;
    pour(auto i = a; i <= b; i++)
      res(j++) = i;
    retourne res;
  }

  /** @brief Intervalle temporel de @f$n@f$ échantillons, en fonction de la fréquence d'échantillonnage.
   *
   *  @param n nombre de points
   *  @param fe fréquence d'échantillonnage
   *  @return Ensemble de valeurs temporelles équiréparties et échantillonées à la fréquence @f$f_e@f$ :
   *  @f[
   *   x_k = \frac{k}{f_e},\quad k = 0,\dots,n-1
   *  @f]
   *
   *  @sa linspace(), intervalle_entier()
   */
  static inline Vecf intervalle_temporel(entier n, float fe)
  {
    retourne linspace(0, (n-1) / fe, n);
  }

  /** @brief Intervalle temporel de durée @f$T@f$, en fonction de la fréquence d'échantillonnage. */
  static inline Vecf intervalle_temporel(float T, float fe)
  {
    soit n  = (entier) ceil(T*fe);
    retourne linspace(0, (n-1) / fe, n);
  }

  /** @brief Loi normale (vecteur colonne).
   *
   *  Génération d'un vecteur (colonne) d'échantillons d'une loi normale (@f$\mathcal{N}(0,1)@f$).
   *
   *  @param n nombre de points à générer
   *  @returns Un vecteur de valeurs aléatoires, échantillonées suivant une loi normale (variance unitaire).
   *
   *  @note Notez que pour générer une loi Gaussienne plus générale (variance différente de 1, et moyenne différente de zéro),
   *  il suffit d'effectuer une mise à l'échelle.
   *  Par exemple, pour une moyenne de 5, et une variance de @f$1/2@f$ :
   *  @code
   *  soit x = 5 + 0.5 * randn(n);
   *  @endcode
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randn
   *  @image html randn.png
   *
   *  @sa randu(entier, float, float), randb()
   */
  extern Vecf randn(entier n);

  /** @brief Loi normale complexe (vecteur colonne) */
  extern Veccf randcn(entier n);

  /** @cond undoc */
  /** @brief Générateur aléatoire utilisé en interne */
  extern std::default_random_engine generateur_aleatoire;
  /** @endcond */


  /** @brief Loi uniforme (vecteur colonne).
   *
   *  @param n Nombre de points à générer
   *  @param a Valeur minimale
   *  @param b Valeur maximale
   *  @returns Un vecteur de valeurs aléatoires, échantillonées suivant une loi uniforme entre a et b.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randu
   *  @image html randu.png
   *
   *  @sa randn(entier), randi(), randb()
   */
  extern Vecf randu(entier n, float a = -1, float b = 1);

  /** @brief Loi aléatoire binaire.
   *
   *  @param n nombre de points à générer
   *  @returns Un vecteur de valeurs aléatoire de égales à 0 ou 1 uniquement.
   *
   *  Cette fonction peut être utile pour générer un train de bits aléatoires.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randb
   *  @image html randb.png
   *
   *
   *  @sa randn(entier), randu(entier,float,float), randi()
   */
  extern Vecb randb(entier n);

  /** @brief Loi aléatoire catégorielle
   *
   *  @param M Nombre de catégories
   *  @param n Nombre de valeurs à générer
   *  @returns Un vecteur de @f$n@f$ entiers entre @f$0@f$ et @f$M-1@f$.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_randi
   *  @image html randi.png
   *
   *  @sa randn(entier), randu(entier,float,float), randb(), randi(entier)
   */
  extern Veci randi(entier M, entier n);

  /** Loi aléatoire catégorielle (calcul d'un seul élément).
   *
   *  @sa randi(entier, entier)
   */
  extern entier randi(entier M);

  /** @brief Loi uniforme (calcul d'un seul élément)
   *
   *  Retourne un nombre entre -1 et 1.
   *
   *  @sa randu(entier, float, float) */
  extern float randu();

  /** @brief Loi normale (calcul d'un seul élément)
   *  @sa randn(entier) */
  extern float randn();

  /** @brief Calcul efficace d'une exponentielle complexe
   *
   *  Cette fonction génère un signal exponentiel, le calcul étant basé sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonométrique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une dérive, tant en gain qu'en phase.
   *
   *  @param f Fréquence normalisée (entre -0,5 et 0,5)
   *  @param n Nombre de points à générer
   *  @returns Vecteur complexe défini par :
   *  @f[
   *   y_k = e^{2\pi\mathbf{i} k f},\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigexp
   *  @image html sigexp.png
   *
   *  @sa sigcos(), sigsin(), sigcar(), sigtri()
   */
  extern Veccf sigexp(float f, entier n);

  /** @brief Calcul efficace d'une sinusoide.
   *
   *  Cette fonction génère un signal sinusoidal, le calcul étant basé sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonométrique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une dérive, tant en gain qu'en phase.
   *
   *  @param f Fréquence normalisée (entre -0,5 et 0,5)
   *  @param n Nombre de points à générer
   *  @returns Vecteur réel défini par :
   *  @f[
   *   y_k = \sin\left(2\pi k f\right),\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigsin
   *  @image html sigsin.png
   *
   *  @sa sigcos(), sigexp(), sigcar(), sigtri()
   */
  extern Vecf sigsin(float f, entier n);

  /** @brief Calcul efficace d'un cosinus
   *
   *  Cette fonction génère un signal sinusoidal, le calcul étant basé sur un oscillateur harmonique
   *  (ce qui beaucoup plus efficace que via les fonctions trigonométrique de la librairie standard).
   *  Notez cependant que si @f$n@f$ est important, il peut y avoir une dérive, tant en gain qu'en phase.
   *
   *  @param f Fréquence normalisée (entre -0,5 et 0,5)
   *  @param n Nombre de points à générer
   *  @returns Vecteur réel défini par :
   *  @f[
   *   y_k = \cos\left(2\pi k f\right),\quad k = 0\dots n-1
   *  @f]
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigcos
   *  @image html sigcos.png
   *
   *  @sa sigsin(), sigexp()
   */
  extern Vecf sigcos(float f, entier n);


  /** @brief Signal triangulaire (périodique).
   *
   *  @param p Période, en nombre d'échantillons.
   *  @param n Nombre d'échantillons à générer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigtri
   *  @image html sigtri.png
   *
   *  @sa sigcar(), sigsin(), sigcos(), sigexp()
   */
  extern Vecf sigtri(entier p, entier n);

  /** @brief Signal carré (périodique)
   *
   *  @param p Période, en nombre d'échantillons.
   *  @param n Nombre d'échantillons à générer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigcar
   *  @image html sigcar.png
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern Vecf sigcar(entier p, entier n);

  /** @brief Impulsion discrète
   *
   *  @param n Dimension du vecteur
   *  @param p Position de l'impulsion
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigimp
   *  @image html sigimp.png
   *
   */
  extern Vecf sigimp(entier n, entier p = 0);


  /** @brief Signal en dent de scie
   *
   *  @param p Période, en nombre d'échantillons.
   *  @param n Nombre d'échantillons à générer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigscie
   *  @image html sigscie.png
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern Vecf sigscie(entier p, entier n);

  /** @brief Sinusoïde modulée par une Gaussienne
   *
   *  @f[
   *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2} \cdot \sin \left(2\pi fk\right)
   *  @f]
   *
   *
   *  @param f Fréquence (normalisée).
   *  @param n Nombre d'échantillons à générer.
   *  @param a Paramètre de décroissance
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_siggsin
   *  @image html siggsin.png
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern Vecf siggsin(float f, entier n, float a = 10);


  /** @brief Impulsion filtrée par une Gaussienne
     *
     *  @f[
     *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2}
     *  @f]
     *
     *
     *  @param n Nombre d'échantillons à générer.
     *  @param a Paramètre de décroissance
     *
     *  @par Exemple
     *  @snippet exemples/src/ex-tsd.cc ex_siggaus
     *  @image html sigsiggaus.png
     *
     *  @sa sigtri(), sigsin(), sigcos(), sigexp(), siggsin()
     */
    extern Vecf siggauss(entier n, float a = 10);

  /** @brief Chirp linéaire ou quadratique.
   *
   *  @f[
   *  x_k = \cos \phi_k, \quad \phi_k = 2 \pi \sum_{i=0}^k f_k
   *  @f]
   *
   *  Les @f$f_k@f$ étant répartis linéairement ou quadratiquement entre @f$f_0@f$ et @f$f_1@f$.
   *
   *  Pour un chirp linéaire :
   *  @f[
   *  f_k = f_0 + (f_1 - f_0) \left(\frac{k}{n-1}\right)
   *  @f]
   *
   *  Pour un chirp quadratique :
   *  @f[
   *  f_k = f_0 + (f_1 - f_0) \left(\frac{k}{n-1}\right)^2
   *  @f]
   *
   *  @param f0 Fréquence initiale (normalisée).
   *  @param f1 Fréquence finale (normalisée).
   *  @param n Nombre d'échantillons à générer.
   *  @param mode 'l' pour linéaire, 'q' pour quadratique.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_sigchirp
   *  @image html sigchirp.png
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp()
   */
  extern Vecf sigchirp(float f0, float f1, entier n, char mode = 'l');


  /** @brief Signal à la fréquence de Nyquist (-1,1,-1,1,etc.).
   *
   *  Construit une sinusoïde de fréquence @f$f_e/2@f$, soit une séquence
   *  alternée -1,1,-1,1,etc.
   *
   *  @param n Nombre d'échantillons à générer.
   *
   *  @par Exemple
   *  @snippet exemples/src/ex-tsd.cc ex_signyquist
   *  @image html signyquist.png
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp(), sigchirp()
   */
  extern Vecf signyquist(entier n);


  /** @brief Configuration oscillateur harmonique */
  struct OHConfig
  {
    /** @brief Fréquence normalisée */
    float freq = 0;

    // Si vrai, applique un shift de la fréquence indiquée
    //bouléen shift = non;
    //float df = 0;
  };

  /** @brief Génération d'un signal exponentiel via un oscillateur harmonique
   *
   *  Cette fonction renvoie une source de données, qui peut être appellée plusieurs fois (génération d'un flux continu d'échantillons, à la différence de @ref sigexp(), qui ne peut générer qu'un nombre fini et prédeterminé d'échantillons).
   *
   *  L'oscillateur en quadrature est basé sur un simple filtre récursif (complexe) d'ordre 1 :
   *  @f[
   *  z_k = z_{k-1} \cdot e^{2\pi\mathbf{i}f}
   *  @f]
   *  qui permet donc de générer un signal exponentiel :
   *  @f[
   *  z_k = e^{2\pi\mathbf{i}kf}
   *  @f]
   *
   *  @note Pour générer une exponentielle de <b>fréquence fixe</b>, cette méthode sera beaucoup plus efficace
   *  que l'utilisation des fonctions trigonométriques de la librairie standard.
   *
   *  @param freq Fréquence normalisée (entre -0,5 et 0,5).
   *  @return Une source de donnée complexe.
   *
   *  @par Exemple
   *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohc
   *  @image html ohc.png
   *
   *  @sa source_ohr(), sigexp()
   */
  extern sptr<Source<cfloat, OHConfig>> source_ohc(float freq);

  /** @brief Génération d'un signal sinusoïdal via un oscillateur harmonique
   *
   *  Cette fonction renvoie une source de données, qui peut être appellée plusieurs fois (génération d'un flux continu d'échantillons, à la différence de @ref sigcos(), qui ne peut générer qu'un nombre fini et prédeterminé d'échantillons).
   *
   *  L'oscillateur est basé sur un oscillateur en quadrature standard (voir @ref source_ohc()), dont
   *  on n'utilise que la partie réelle (cosinus).
   *
   *
   *  @param freq Fréquence normalisée (entre -0,5 et 0,5).
   *  @return Une source de donnée réelle.
   *
   *  @par Exemple
   *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohr
   *  @image html ohr.png
   *
   *  @sa source_ohc(), sigcos()
   */
  extern sptr<Source<float, OHConfig>> source_ohr(float freq);

  struct OLUT
  {
    OLUT(entier resolution = 1024);
    cfloat step(float θ);

    _PIMPL_
  };

  /** @} */



}



#endif

