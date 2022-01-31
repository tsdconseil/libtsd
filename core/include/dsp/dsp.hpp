#pragma once


#include "tsd/tsd.hpp"

namespace dsp {

/** @brief @f$\pi@f$ value, with 32 bits floating point accuracy. */
auto pi = tsd::pi;

auto π  = tsd::π;

/** @brief Short-cut for the complex 32 bits floating point type. */
using cfloat = tsd::cfloat;

/** @brief Short-cut for the complex 32 bits floating point type. */
using cdouble = tsd::cdouble;

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
  using Tableau = tsd::Tableau<t>;

/** @brief Shortcut for a column vector with generic (templated) data type. */
template<typename t>
  using Vector = tsd::Vecteur<t>;

/** @brief Input floating point array (const) */
using IArrayXf   = tsd::IArrayXf;

/** @brief Input complex floating point array (const) */
using IArrayXXcf = tsd::IArrayXXcf;

/** @brief Input complex floating point column vector (const) */
using IArrayXcf  = tsd::IArrayXcf;

/** @brief Vecteur flottant en sortie */
using OArrayXf   = tsd::OArrayXf;
using IArrayXXf  = tsd::IArrayXXf;
using OArrayXXf  = tsd::OArrayXXf;

/** @brief Boolean array */
using ArrayXXb = tsd::ArrayXXb;

/** @brief Boolean column array */
using ArrayXb  = tsd::ArrayXb;

using IArrayXXb = tsd::IArrayXXb;
using IArrayXb  = tsd::IArrayXb;

/** @brief Short-cut for a shared pointer. */
template<typename T>
  using sptr = std::shared_ptr<T>;


/** @brief Vertical concatenation of 2 vectors.
 *
 * <h3>Vertical concatenation</h3>
 *
 * This function merge 2 column vectors:
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
 * (assuming a and b have respectively @f$n@f$ et @f$m@f$ elements).
 * The | operator has been redefined to make the same operation.
 *
 * @param a First vector
 * @param b Second vector
 * @returns Vertical concatenation
 *
 * @par Example
 * @code
 *  ArrayXf a = ArrayXf::Ones(1);
 *  ArrayXf b = ArrayXf::Zero(5);
 *  ArrayXf c = vconcat(a, b); // (or ArrayXf c = a | b)
 * @endcode
 *
 */
template<typename T1, typename T2>
  auto vconcat(const Eigen::ArrayBase<T1> &a, const Eigen::ArrayBase<T2> &b)
{
  return tsd::vconcat(a, b);
}

template<typename T1, typename T2>
  auto operator |(const Eigen::ArrayBase<T1> &a, const Eigen::ArrayBase<T2> &b)
{
  return vconcat(a, b);
}


/** @brief Vector rotation.
 *
 *  <h3>Vector rotation</h3>
 *
 *  @param x Input column vector
 *  @param d Integer number of steps (positive or negative)
 *  @returns Data vector with rtation of @f$d@f$ steps (modulo the vector dimension):
 *  @f[
 *   y_k = x_{k + d\ [N]}
 *  @f]
 *  @f$N@f$ being the vector length.
 *
 *  @par Example
 *  @code
 *  ArrayXf x = linspace(0, 5, 6);
 *  ArrayXf y = rotation_vec(x, 2);
 *  ArrayXf yref(6);
 *  yref << 2, 3, 4, 5, 0, 1;
 *  tsd_assert(y.isApprox(yref));
 *  @endcode
 */
template<typename T>
Vector<T> rotation_vec(const Vector<T> &x, int d)
{
  return tsd::rotation_vec(x, d);
}


/** @brief Difference between 2 successive elements of a vector.
 *
 *  <h3>Difference between 2 successive elements of a vector</h3>
 *
 *  @returns A column vector of length @f$n-1@f$:
 *  @f[
 *  y_k = x_{k+1} - x_k,\quad k = 0,\dots, n-2
 *  @f]
 *
 *  @par Example
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
  return tsd::diff(x);
}


/** @brief Accumulated sum of a vector.
 *  <h3>Accumulated sum of a vector</h3>
 *
 *  Compute a vector of identical dimension as the input one and
 *  equal to the progressive accumulation of input values:
 *  @f[
 *  y_n = \sum_{k=0}^n x_k,\quad k = 0\dots N-1
 *  @f]
 *
 *  @par Example
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
  return tsd::cumsum(x);
}

/** @brief Fix phase jumps.
 *
 *  <h3>Fix phase jumps</h3>
 *
 *  @param x Vector (typically with angle values).
 *  @param r Defines the equivalence class on @f$x@f$
 *  @returns A vector @f$y@f$ such as @f$y_k = x_k + k \cdot r @f$, and with the less possible number of discontinuities (no discontinuity greater than @f$r/2@f$ in absolute value).
 *
 *  @par Example:
 *  @snippet exemples/src/ex-tsd.cc exemple_unwrap
 *  @image html unwrap.png width=800px
 *
 *  @sa modulo_pm_π(), modulo_2π()
 *
 */
inline Eigen::ArrayXf phase_unwrap(const ArrayXf &x, float r = 2*π)
{
  return tsd::déplie_phase(x, r);
}


/** @brief Find true values indexes in a boolean vector.
 *
 *  <h3>Find true values indexes in a boolean vector</h3>
 *
 *  @param x Boolean vector
 *  @returns Vector containing the indexes for which the input vector elements are true.
 *
 *  @par Example
 *  @code
 *  ArrayXf x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
 *  auto idx = find(x >= 0); // idx = {2, 3, 4}
 *  @endcode
 *
 *  @sa find_first()
 */
inline std::vector<int> find(IArrayXb x)
{
  return tsd::trouve(x);
}

/** @brief Find the first true element index from a boolean vector.
 *
 *  <h3>Find the first true element index from a boolean vector.</h3>
 *
 *
 *  @par Example
 *  @code
 *  ArrayXf x = linspace(-2, 2, 5); // -2, -1, 0, 1, 2
 *  auto idx = find_first(x >= 0); // idx = 2
 *  @endcode
 *
 *  @sa find()
 */
inline int find_first(IArrayXb x)
{
  return tsd::trouve_premier(x);
}



// Equivalent de : y = x(1:pas:$)
/** @brief Decimation of a column vector
 *
 *  <h3>Decimation of a column vector</h3>
 *
 *  From a vector @f$(x_k), k=0\dots n-1@f$, returns
 *  a subset of it, decimated by a factor of @f$R@f$:
 *  @f[
 *  y_k = x_{kR},\quad k =0 \dots \left\lfloor N/R\right\rfloor -1
 *  @f]
 *
 *  @par Example:
 *  @snippet exemples/src/ex-tsd.cc ex_sousech
 *  @image html sousech.png width=800px
 *
 *  @sa upsample()
 *
 */
template<typename T>
 Vector<T> subsample(const Vector<T> &x, int pas)
 {
  return tsd::sousech(x, pas);
 }

/** @brief Oversampling of a column vector.
 *
 *  <h3>Oversampling of a column vector</h3>
 *
 *  From a vector  @f$(x_k), k=0\dots n-1@f$, returns
 *  a upsampled vector of size @f$Rn@f$, by inserting zeroes:
 *  @f[
 *  y_k = \begin{cases} x_{k/R} & \textrm{if } k\textrm{ is a multiple of } R,\\
 *   0 & \textrm{otherwise.}\end{cases}
 *  @f]
 *
 *  @par Example:
 *  @snippet exemples/src/ex-tsd.cc ex_surech
 *  @image html surech.png width=800px
 *
 *  @sa subsample()
 *
 */
template<typename T>
 Vector<T> upsample(const Vector<T> &x, int R)
 {
  return tsd::surech(x, R);
 }


/** @brief Decibels to linear conversion.
 *
 *  <h3>Decibels to linear conversion</h3>
 *
 *  @param x Value in dB
 *  @returns Value in linear units:
 *  @f[
 *    y = 10^{x/10}
 *  @f]
 *
 *  @sa pow2db()
 */
template<typename T>
  auto db2pow(const T &x)
{
  return tsd::db2pow(x);
}

/** @brief Linear to decibel conversion.
 *
 *  <h3>Linear to decibel conversion</h3>
 *
 *  @param x Value in linear units
 *  @returns Valeur en dB:
 *  @f[
 *    y = 10 \log_{10}(x)
 *  @f]
 *
 *  @sa db2pow()
 */
template<typename T>
  auto pow2db(const T &x)
{
  return tsd::pow2db(x);
}



/** @brief Returns the smallest power of 2 greater or equal to i.
 *
 * <h3>Smallest power of 2 greater or equal to some integer</h3>
 *
 * Computes the smallest power of 2 greater or equal to @f$i@f$:
 * @f[
 *  y = \min_{j\geq i,\ j = 2^k} j
 * @f]
 *
 * @par Example
 * @code
 * tsd_assert(next_power_of_2(3) == 4);
 * @endcode
 *
 */
inline int next_power_of_2(unsigned int i)
{
  return tsd::prochaine_puissance_de_2(i);
}


/** @brief Ajoute des zéros à la fin d'un des deux vecteurs de manière à ce qu'ils aient la même dimension.
 *
 *  <h3>Complétion de vecteurs avec des zéros</h3>
 *
 *  Un des deux vecteurs est complétés avec des zéros de manière à ce que les deux vecteurs
 *  résultants aient la même longueur.
 *
 *  @param x  Premier vecteur
 *  @param y  Deuxième vecteur
 *  @param p2 Si vrai, arrondi la dimension des vecteurs de sortie à la prochaine puissance de 2.
 *  @return Un tuple de deux vecteurs de mêmes dimensions, et commençant respectivement par x et y (des zéros sont ajoutés à la fin).
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
  std::tuple<Vector<typename D1::Scalar>, Vector<typename D2::Scalar>> pad_zeros(const Eigen::ArrayBase<D1> &b1, const Eigen::ArrayBase<D2> &b2)
{
  return tsd::pad_zeros(b1, b2);
}

///////////////////////

/** @brief Equidistant point interval.
 *
 *  <h3>Equidistant point interval</h3>
 *
 *
 *  Compute @f$n@f$ equidistant points between @f$a@f$ and @f$b@f$.
 *
 *  @param a Initial point (first value)
 *  @param b Final point (last value)
 *  @param n Number of points
 *  @returns Interval defined as:
 *  @f[
 *  x_k = a + k \cdot \frac{b-a}{n-1},\ k = 0,\dots,n-1
 *  @f]
 *
 *  @par Example 1: vector with the @f$n@f$ first integers:
 *  @code
 *  ArrayXf t = linspace(0, n-1, n)
 *  // t = 0, 1, 2, ..., n-1
 *  @endcode
 *  @par Example 2: vector of @f$n@f$ time points, with a sampling frequency of @f$f_s@f$ Hz:
 *  @code
 *  ArrayXf t = linspace(0, (n-1)/fs, n)
 *  @endcode
 */
static inline auto linspace(float a, float b, unsigned int n)
{
  return tsd::linspace(a,b,n);
}

/** @brief Integer interval
 *
 *  <h3>Integer interval</h3>
 *
 *  @return Vector of the  @f$b-a+1@f$ integers between @f$a@f$ and @f$b@f$:
 *  @f[
 *   x_k = a + k,\quad k = 0,\dots,(b-a+1)
 *  @f]
 *
 *  @sa linspace(), trange()
 */
static inline ArrayXi irange(int a, int b)
{
  return tsd::intervalle_entier(a, b);
}

/** @brief Time interval, with sampling frequency specified.
 *
 *  <h3>Time interval, with sampling frequency specified</h3>
 *
 *  @param n Number of points
 *  @param fs Sampling frequency
 *  @return Set of equidistant points :
 *  @f[
 *   x_k = \frac{k}{f_s},\quad k = 0,\dots,n-1
 *  @f]
 *
 *  @sa linspace(), irange()
 */
static inline ArrayXf trange(unsigned int n, float fs)
{
  return tsd::intervalle_temporel(n, fs);
}



/** @brief Loi normale (vecteur colonne)
 *  @param n nombre de points à générer
 *  @returns Un vecteur de valeurs aléatoires, échantillonées suivant une loi normale (variance unitaire).
 *  @sa randu */
extern ArrayXf randn(unsigned int n){return tsd::randn(n);}

extern ArrayXXf randn_2d(unsigned int n, unsigned int m){return tsd::randn_2d(n,m);}

/** @brief Loi uniforme, intervalle [0,1] (tableau 2d) */
extern ArrayXXf randu(unsigned int n, unsigned int m);

/** @brief Loi uniforme (vecteur colonne)
 *  @param n nombre de points à générer
 *  @returns Un vecteur de valeurs aléatoires, échantillonées suivant une loi uniforme entre 0 et 1.
 *  @sa randn */
extern ArrayXf randu(unsigned int n);

/** @brief Loi aléatoire binaire
 *  @param n nombre de points à générer
 *  @returns Un vecteur de valeurs aléatoire de valeurs 0 ou 1 uniquement.
 *  Cette fonction peut être utile pour générer un train de bits aléatoires.
 */
static inline ArrayXb randb(int n)
{
  return tsd::randb(n);
}

/** @brief Loi aléatoire catégorielle
 *  Entiers entre 0 et M-1 */
extern ArrayXi randi(unsigned int M, unsigned int n);

/** @brief Calcule la longueur d'un vecteur colonne */
static inline unsigned int length(const ArrayXf &x)
{
  return x.rows();
}







extern ArrayXcf polar(const ArrayXf &θ);
extern ArrayXcf polar(const ArrayXf &ρ, const ArrayXf &θ);







inline auto subarray1d(IArrayXf &x, unsigned int dec, unsigned int n, unsigned int pas)
{
  return tsd::subarray1d(x, dec, n, pas);
}

/** @brief Calcule le carré d'un scalaire */
template<typename T>
  T square(const T &v) {return tsd::carré<T>(v);}



/** @brief Structure configurable
 *  @tparam C Type de configuration */
template<typename C>
  using Configurable = tsd::Configurable<C>;

using Void = tsd::Void;

/** @brief Puit (réceptacle de données) générique et non configurable */
template<typename Te>
  using SinkGen = tsd::SinkGen<Te>;


/** @brief Puit (réceptacle de donnée) générique et configurable */
template<typename Te, typename Tc = Void>
  using Sink = tsd::Sink<Te, Tc>;

/** @brief Source de données, configurable
 *  @tparam Ts Type de données à produire (float, cfloat, ...)
 *  @tparam C  Type de configuration */
template<typename Ts, typename C = Void>
  using Source = tsd::Source<Ts, C>;

/** @brief Generic, configuration-less, filter
 *  @tparam Te Input type
 *  @tparam Ts Output type */
template<typename Te, typename Ts = Te>
using FiltreGen = tsd::FiltreGen<Te, Ts>;


/** @brief Generic, run-time configurable, filter
 *  @tparam  Te Input type
 *  @tparam  Ts Output type
 *  @tparam  Tc Configuration type
 */
template<typename Te, typename Ts = Te, typename Tc = Void>
using Filtre = tsd::Filtre<Te, Ts, Tc>;



/** @brief Ré-échantillonnage d'un signal suivant un facteur arbitraire
 *  @param lom Facteur de ré-échantillonnage (@f$f_e^{(2)} / f_e^{(1)}@f$)
 *  @returns Signal ré-échantilloné
 *
 *  Un filtre anti-repliement est automatiquement calculé.
 */
extern ArrayXf resample(IArrayXf x, float lom);
extern ArrayXcf resample(IArrayXcf x, float lom);





/** @brief A partir d'un flux de données entrant, découpé en paquet de dimensions quelconques, produit un flux de données sous forme de <b>paquets de dimension fixe</b>.
    @param N dimension des paquets de sortie
    @param callback Fonction utilisateur qui sera appelée à chaque fois qu'un nouveau paquet de dimension N est disponible.
    @returns Un Sink acceptant des vecteurs de type @p T et reconfigurable (N).
    @tparam T Type de données d'entrée / sortie (e.g. float, cfloat, double, ...). */
template<typename T>
  sptr<Sink<T,int>> creation_tampon(int N,
      std::function<void (const Vector<T> &)> callback);




}



