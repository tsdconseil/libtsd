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
using Eigen::ArrayXXf;
using Eigen::ArrayXcf;
using Eigen::ArrayXXcf;
using Eigen::ArrayXi;
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

/** @brief Vertical concatenation of 2 vectors. */
template<typename T1, typename T2>
  auto vconcat(const Eigen::ArrayBase<T1> &a, const Eigen::ArrayBase<T2> &b)
{
  return tsd::vconcat(a, b);
}


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
static inline ArrayXf irange(int a, int b)
{
  return tsd::irange(a, b);
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
  return tsd::trange(n, fs);
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
extern ArrayXf randb(unsigned int n);

/** @brief Loi aléatoire catégorielle
 *  Entiers entre 0 et M-1 */
extern ArrayXi randi(unsigned int M, unsigned int n);

/** @brief Calcule la longueur d'un vecteur colonne */
static inline unsigned int length(const ArrayXf &x)
{
  return x.rows();
}

/** @brief Rotation d'un vecteur
 *  @param x vecteur colonne d'entrée
 *  @param d nombre entier de pas (peut être négatif)
 *  @returns vecteur des données décalée de d pas (modulo la dimension du vecteur) :
 *  @f[
 *   y_k = x_{k + d\ [N]},\ N \textrm{ étant la dimension du vecteur.}
 *  @f]
 */
template<typename T>
Vector<T> rotation_vec(const Vector<T> &x, int d)
{
  return tsd::rotation_vec(x, d);
}

extern ArrayXf diff(const ArrayXf &x);
extern Eigen::ArrayXi diff(const ArrayXi &x);
extern Eigen::ArrayXi cumsum(const ArrayXi &x);
extern ArrayXf cumsum(const ArrayXf &x);

extern std::vector<int> find(IArrayXb x);

extern ArrayXcf polar(const ArrayXf &θ);
extern ArrayXcf polar(const ArrayXf &ρ, const ArrayXf &θ);



// Equivalent de : y = x(1:pas:$)
template<typename T>
 Vector<T> sousech(const Vector<T> &x, int pas)
 {
  return tsd::sousech(x, pas);
 }

template<typename T>
 Vector<T> surech(const Vector<T> &x, int R)
 {
  return tsd::surech(x, R);
 }

/** Ajoute des zéros à la fin d'un des deux vecteurs de manière à ce qu'ils aient la même dimension */
template<typename D1, typename D2>
  std::tuple<Vector<typename D1::Scalar>, Vector<typename D2::Scalar>> pad_zeros(const Eigen::ArrayBase<D1> &b1, const Eigen::ArrayBase<D2> &b2)
{
  return tsd::pad_zeros(b1, b2);
}


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



