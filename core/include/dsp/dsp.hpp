#pragma once


#include "tsd/tsd.hpp"

namespace dsp {


/** @addtogroup tsd
 *  @{ */


#define dsp_assert      tsd_assert
#define dsp_assert_msg  tsd_assert_msg

/** @brief @f$\pi@f$ value, with 32 bits floating point accuracy. */
static const auto pi = tsd::pi;

static const auto π  = tsd::π;

static const auto π_f  = tsd::π_f;

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
 * This function merges 2 column vectors:
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
 Vector<T> downsample(const Vector<T> &x, int pas)
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
 *  @sa downsample()
 *
 */
template<typename T>
 Vector<T> upsample(const Vector<T> &x, int R)
 {
  return tsd::surech(x, R);
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


/** @brief Add zeros at the end of one of two vectors such as they have the same length.
 *
 *  <h3>Zero padding of vectors</h3>
 *
 *  One of the two vectors is padded with zeros such as the two vectors have the same length.
 *
 *  @param x  First vector
 *  @param y  Second vector
 *  @param p2 If true, make sure the two new vectors length is a power of two.
 *  @return A tuple of two vectors with the same length, and beginnning respectively by x and y (zeros are added at the end).
 *
 *  @par Example
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
  std::tuple<Vector<typename D1::Scalar>, Vector<typename D2::Scalar>>
    pad_zeros(const Eigen::ArrayBase<D1> &x,
              const Eigen::ArrayBase<D2> &y,
              bool p2 = false)
{
  return tsd::pad_zeros(x, y, p2);
}


/** @brief Computes the square of a scalar value. */
template<typename T>
  T square(const T &v) {return tsd::carré<T>(v);}


/** @brief Abstract configurable structure
 *  @tparam C Configuration type */
template<typename C>
  using Configurable = tsd::Configurable<C>;

using Void = tsd::Void;

/** @brief Abstract data sink. */
template<typename Te>
  using SinkGen = tsd::SinkGen<Te>;


/** @brief Abstract, configurable, data sink. */
template<typename Te, typename Tc = Void>
  using Sink = tsd::Sink<Te, Tc>;


/** @brief Abstract data source */
template<typename Ts>
  using SourceGen = tsd::SourceGen<Ts>;

/** @brief Abstract, configurable, data source
 *  @tparam Ts Type de données à produire (float, cfloat, ...)
 *  @tparam C  Type de configuration */
template<typename Ts, typename C = Void>
  using Source = tsd::Source<Ts, C>;

/** @brief Generic data filter
 *  @tparam Te Input type
 *  @tparam Ts Output type */
template<typename Te, typename Ts = Te>
using FilterGen = tsd::FiltreGen<Te, Ts>;


/** @brief Generic, run-time configurable, filter
 *  @tparam  Te Input type
 *  @tparam  Ts Output type
 *  @tparam  Tc Configuration type
 */
template<typename Te, typename Ts = Te, typename Tc = Void>
using Filter = tsd::Filtre<Te, Ts, Tc>;





/** @brief  Finite signal resampling, with arbitrary ratio.
 *
 *  <h3>Finite signal resampling</h3>
 *
 *  This function resample a signal with a decimation (@f$r < 1@f$)
 *  or interpolation (@f$r > 1@f$) ratio.
 *
 *  @param x    Input signal
 *  @param r    Resampling factor (@f$f_e^{(2)} / f_e^{(1)}@f$)
 *  @returns    Resampled signal
 *
 *  @note A filter (or even a cascad of filters) is automatically inserted (before decimation or after interpolation) so
 *  as to avoir spectrum aliasing.
 *
 *  @sa filter_resampling() (to process infinite time signal), resample_freq() (zero-delay resampling).
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_resample
 *  @image html resample.png width=800px
 *
 *  @note Because of the anti-aliasing filters, the output signal will be delayed compared to the input one.
 *  To resample without intruducing delay, one can use @ref resample_freq().
 */
template<typename Derived>
  auto resample(const Eigen::ArrayBase<Derived> &x, float r)
{
  return tsd::rééchan(x, r);
}




/** @brief From an input data stream, formatted in arbitrary packet sizes,
 *  produce a data stream with a <b>fixed packet length</b>.
 *
 *  <h3>Preparation of a fixed packet length data stream</h3>
 *
 *  From an input data stream, formatted in arbitrary packet sizes,
 *  produce a data stream with a <b>fixed packet length</b>.
 *
 *  @image html tampon.png width=800px
 *
 *  @param N        Output packets length
 *  @param callback User function which will be called for each output packet of size @f$N@f$
 *  @returns        A data sink accepting vectors of type @p T and of arbitrary length.
 *  @tparam T       Input / output data type (e.g. float, cfloat, double, ...).
 *
 *
 *  @par Example: preparation of packets for a FFT
 *  @snippet exemples/src/ex-tsd.cc ex_tampon
 *  Note that 10000 not being a multiple of 512,
 *  the last samples are not processed (they would be if the method t->step())
 *  is called with more samples).
 *
 *
 */
template<typename T>
  sptr<Sink<T,int>> buffer_new(int N,
      std::function<void (const Vector<T> &)> callback)
{
  return tsd::tampon_création(N, callback);
}




/** @brief Modulo with result in @f$[0,m[@f$ interval.
 *
 * <h3>Modulo with result in @f$[0,m[@f$ interval</h3>
 *
 * Contrary to the standard function <code>fmod(x, m)</code>, which returns
 * a value between @f$-m@f$ and @f$m@f$, this function
 * returns a value between @f$0@f$ and @f$m@f$:
 *  @f[
 *  y = x + k\cdot m,\quad k\in\mathbb{Z},\ y\in \left[0,m\right[
 *  @f]
 *
 *  @par Example
 *  @code
 *  auto y = modulo(1.5, 1.0); // -> y = 0.5
 *  @endcode
 *
 *  @sa modulo_2π(), modulo_pm_π()
 *
 */
template<typename T>
T modulo(T x, T m)
{
  return tsd::modulo(x, m);
}


/** @brief Computes @f$a@f$ modulo @f$2\pi@f$, the result being in the @f$\left[0,2\pi\right[@f$ interval.
 *
 *
 *  <h3>Modulo @f$2\pi@f$</h3>
 *  The result is in the @f$\left[0,2\pi\right[@f$ interval:
 *  @f[
 *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[0,2\pi\right[
 *  @f]
 *
 *  @par Example
 *  @code
 *  tsd_assert(modulo_2π(2*π+1e-5) == 1e-5);
 *  @endcode
 *
 *  @sa modulo_pm_π(), phase_unwrap(), modulo()
 *
 **/
template<typename T>
inline T modulo_2π(T x)
{
  return tsd::modulo_2π(x);
}

/** @brief Computes @f$a@f$ modulo @f$2\pi@f$, the result being in the @f$\left[-\pi,\pi\right[@f$ interval.
 *
 *  <h3>Modulo @f$2\pi@f$</h3>
 *  The result is in the @f$\left[-\pi,\pi\right[@f$ interval:
 *  @f[
 *  y = x + k\cdot 2\pi,\quad k\in\mathbb{Z},\ y\in \left[-\pi,\pi\right[
 *  @f]
 *
 *  @par Example
 *  @code
 *  tsd_assert(modulo_pm_π(2*π-1e-5) == -1e-5);
 *  @endcode
 *
 *  @sa modulo_2π(), phase_unwrap(), modulo()
 *
 **/
template<typename T>
T modulo_pm_π(T x)
{
  return tsd::modulo_pm_π(x);
}


/** @brief Degrees to radians conversion.
 *
 *  <h3>Degrees to radians conversion</h3>
 *
 *  @f[
 *  y = \frac{\pi\cdot x}{180}
 *  @f]
 *
 *  @par Example
 *  @code
 *  tsd_assert(abs(deg2rad(45) - π/4) < 1e-15);
 *  @endcode
 *
 *  @sa rad2deg()
 */
template<typename T>
auto deg2rad(T degrees)
{
  return tsd::deg2rad(degrees);
}

/** @brief Radians to degrees conversion.
 *
 *  <h3>Radians to degrees conversion</h3>
 *
 *  @f[
 *  y = \frac{180\cdot x}{\pi}
 *  @f]
 *
 *  @par Example
 *  @code
 *  tsd_assert(abs(rad2deg(π/4) - 45) < 1e-15);
 *  @endcode
 *
 *  @sa deg2rad()
 */
template<typename T>
T rad2deg(T radians)
{
  return tsd::rad2deg(radians);
}

/** @} */

/** @addtogroup tsd-gen
 *  @{ */


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
 *
 *  @sa logspace()
 */
static inline auto linspace(float a, float b, unsigned int n)
{
  return tsd::linspace(a,b,n);
}


/** @brief Interval of logarithmicly equidistant points (geometric serie).
 *
 *  <h3>Interval of logarithmicly equidistant points</h3>
 *
 *  Computes @f$n@f$ points logarithmicly equidistants between @f$a@f$ and @f$b@f$,
 *  that is, the following geometric serie:
 *
 *  @f[
 *  x_k = 10^{a + k  \frac{b-a}{n-1}},\quad k = 0,\dots,n-1
 *  @f]
 *
 *
 *  @param a Logarithm (base 10) of the initial point
 *  @param b Logarithm (base 10) of the final point
 *  @param n Number of points
 *  @returns Geometric serie
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_logspace
 *  @image html logspace.png width=600px
 *
 *  @sa linspace()
 */
static inline auto logspace(float a, float b, int n)
{
  return tsd::logspace(a, b, n);
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



/** @brief Normal law (column vector).
 *
 *  <h3>Normal law</h3>
 *  Generation of a samples vector for a normal law (@f$\mathcal{N}(0,1)@f$).
 *
 *  @param n Number of points to generate.
 *  @returns Random vector.
 *
 *  @note Note that so as to generate a more general Gaussian law (variance not unitary, and / or mean not null),
 *  one has just to scale the result of a normal law.
 *  For instance, for a mean of 5, and a variance of @f$1/2@f$ :
 *  @code
 *  ArrayXf x = 5 + 0.5 * randn(n);
 *  @endcode
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_randn
 *  @image html randn.png width=600px
 *
 *  @sa randu(), randb()
 */
static inline ArrayXf randn(unsigned int n)
{
  return tsd::randn(n);
}

/** @cond private */

/** @brief Normal law (2d array)
 *  @sa randu() */
static inline ArrayXXf randn_2d(unsigned int n, unsigned int m){return tsd::randn_2d(n,m);}



/** @brief Loi uniforme, intervalle [0,1] (tableau 2d) */
static inline ArrayXXf randu_2d(unsigned int n, unsigned int m){return tsd::randu_2d(n, m);}


/** @endcond */

/** @brief Uniform law (column vector).
 *
 *  <h3>Uniform law</h3>
 *
 *  @param n Number of points to generate.
 *  @returns A vector of random values, sampling according to a uniform law between 0 and 1.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_randu
 *  @image html randu.png width=600px
 *
 *  @sa randn(), randi(), randb()
 */
static inline ArrayXf randu(unsigned int n){return tsd::randu(n);}

/** @brief Random binary vector.
 *
 *  <h3>Random binary vector</h3>
 *
 *  @param n Number of points to generate.
 *  @returns A vector of random values, 0 or 1.
 *
 *  This function can be used to generate a random bit sequence.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_randb
 *  @image html randb.png width=600px
 *
 *
 *  @sa randn(), randu(), randi()
 */
static inline ArrayXb randb(int n)
{
  return tsd::randb(n);
}

/** @brief Categorial random vector
 *
 *  <h3>Categorial random vector</h3>
 *
 *  @param M Number of categories.
 *  @param n Number of points to generate.
 *  @returns A vector of @f$n@f$ integers between @f$0@f$ et @f$M-1@f$.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_randi
 *  @image html randi.png width=600px
 *
 *  @sa randn(), randu(), randb()
 */
static inline ArrayXi randi(int M, int n)
{
  return tsd::randi(M, n);
}


/** @brief Efficient computing of a complex exponential.
 *
 *  <h3>Generation of a complex exponential</h3>
 *
 *  This function generates a complex exponential signal,
 *  based on an harmonic oscillator
 *  (which is far more efficient than by using the trigonometric functions of the standard library).
 *  Note however that if @f$n@f$ is a big number, there can be some divergence, both in gain and in phase.
 *
 *  @param f Normalized frequency of the exponential (between -0.5 et 0.5)
 *  @param n Number of points to generate.
 *  @returns Complex vector defined by:
 *  @f[
 *   y_k = e^{2\pi\mathbf{i} k f},\quad k = 0\dots n-1
 *  @f]
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigexp
 *  @image html sigexp.png width=600px
 *
 *  @sa sigcos(), sigsin(), sigcar(), sigtri()
 */
static inline ArrayXcf sigexp(float f, int n)
{
  return tsd::sigexp(f, n);
}

/** @brief Efficient computing of a sinusoid.
 *
 *  <h3>Generation of a sinusoid</h3>
 *
 *  This function generates a real sinusoidal signal,
 *  based on an harmonic oscillator
 *  (which is far more efficient than by using the trigonometric functions of the standard library).
 *  Note however that if @f$n@f$ is a big number, there can be some divergence, both in gain and in phase.
 *
 *  @param f Normalized frequency of the sinusoid (between -0.5 et 0.5)
 *  @param n Number of points to generate.
 *  @returns Real vector defined by:
 *  @f[
 *   y_k = \sin\left(2\pi k f\right),\quad k = 0\dots n-1
 *  @f]
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigsin
 *  @image html sigsin.png width=600px
 *
 *  @sa sigcos(), sigexp(), sigcar(), sigtri()
 */
static inline ArrayXf sigsin(float f, int n)
{
  return tsd::sigsin(f, n);
}

/** @brief Efficient computing of a cosinusoid.
 *
 *  <h3>Generation of a cosinusoid</h3>
 *
 *  This function generates a real cosinusoidal signal,
 *  based on an harmonic oscillator
 *  (which is far more efficient than by using the trigonometric functions of the standard library).
 *  Note however that if @f$n@f$ is a big number, there can be some divergence, both in gain and in phase.
 *
 *  @param f Normalized frequency of the cosinus (between -0.5 et 0.5)
 *  @param n Number of points to generate.
 *  @returns Real vector defined by:
 *  @f[
 *   y_k = \cos\left(2\pi k f\right),\quad k = 0\dots n-1
 *  @f]
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigcos
 *  @image html sigcos.png width=600px
 *
 *  @sa sigsin(), sigexp(), sigcar(), sigtri()
 */
static inline ArrayXf sigcos(float f, int n)
{
  return tsd::sigcos(f, n);
}

/** @brief Generation of a triangular, periodic signal.
 *
 *  <h3>Triangular, periodic signal</h3>
 *
 *  @param p Period, in number of samples.
 *  @param n Number of points to generate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigtri
 *  @image html sigtri.png width=600px
 *
 *  @sa sigcar(), sigsin(), sigcos(), sigexp()
 */
static inline ArrayXf sigtri(int p, int n)
{
  return tsd::sigtri(p, n);
}


/** @brief Generation of a square, periodic signal.
 *
 *  <h3>Square, periodic signal</h3>
 *
 *  @param p Period, in number of samples.
 *  @param n Number of points to generate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigcar
 *  @image html sigcar.png width=600px
 *
 *  @sa sigtri(), sigsin(), sigcos(), sigexp()
 */
static inline ArrayXf sigsquare(int p, int n)
{
  return tsd::sigcar(p, n);
}


/** @brief Discret impulsion.
 *
 *  <h3>Discret impulsion</h3>
 *
 *  @param n Number of points to generate.
 *  @param p Impulsion position (in samples).
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigimp
 *  @image html sigimp.png width=600px
 *
 */
static inline ArrayXf sigimp(int n, int p = 0)
{
  return tsd::sigimp(p, n);
}

/** @brief Sawtooth signal.
 *
 *  <h3>Sawtooth signal</h3>
 *
 *  @param p Period, in number of samples.
 *  @param n Number of points to generate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigscie
 *  @image html sigscie.png width=600px
 *
 *  @sa sigtri(), sigsin(), sigcos(), sigexp()
 */
static inline ArrayXf sigsawtooth(int p, int n)
{
  return tsd::sigscie(p, n);
}


/** @brief Sinusoïd modulated by a Gaussian
 *
 *  <h3>Sinusoïd modulated by a Gaussian</h3>
 *
 *  @f[
 *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2} \cdot \sin \left(2\pi fk\right)
 *  @f]
 *
 *
 *  @param f Normalized frequency
 *  @param n Number of points to generate.
 *  @param a Decrease rate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_siggsin
 *  @image html siggsin.png width=600px
 *
 *  @sa sigtri(), sigsin(), sigcos(), sigexp()
 */
static inline ArrayXf siggsin(float f, int n, float a = 10)
{
  return tsd::siggsin(f, n, a);
}


/** @brief Gaussian filtered impulse.
   *
   *  <h3>Gaussian filtered impulse</h3>
   *
   *  @f[
   *  x_k = e^{-a \left(\frac{k-N/2}{N/2}\right)^2}
   *  @f]
   *
   *
   *  @param n Number of points to generate.
   *  @param a Decrease rate.
   *
   *  @par Example
   *  @snippet exemples/src/ex-tsd.cc ex_siggaus
   *  @image html sigsiggaus.png width=600px
   *
   *  @sa sigtri(), sigsin(), sigcos(), sigexp(), siggsin()
   */
static inline ArrayXf siggauss(int n, float a = 10)
{
  return tsd::siggauss(n, a);
}


/** @brief Linear chirp.
 *
 *  <h3>Linear chirp</h3>
 *
 *  @f[
 *  x_k = \cos \phi_k, \quad \phi_k = 2 \pi \sum_{i=0}^k f_k
 *  @f]
 *
 *  The @f$f_k@f$ being lineairly distribued between @f$f_0@f$ and @f$f_1@f$.
 *
 *  @param f0 Initial frequency (normalized).
 *  @param f1 Final frequency (normalized).
 *  @param n  Number of points to generate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigchirp
 *  @image html sigchirp.png width=600px
 *
 *  @sa sigtri(), sigsin(), sigcos(), sigexp(), sigchirp2()
 */
static inline ArrayXf sigchirp(float f0, float f1, int n)
{
  return tsd::sigchirp(f0, f1, n);
}



/** @brief Quadratic chirp.
 *
 *  <h3>Quadratic chirp</h3>
 *
 *  @f[
 *  x_k = \cos \phi_k, \quad \phi_k = 2 \pi \sum_{i=0}^k f_k
 *  @f]
 *
 *  The @f$f_k@f$ being quadraticaly distribued between @f$f_0@f$ and @f$f_1@f$:
 *  @f[
 *  f_k = f_0 + (f_1 - f_0) \left(\frac{k}{n-1}\right)^2
 *  @f]
 *
 *  @param f0 Initial frequency (normalized).
 *  @param f1 Final frequency (normalized).
 *  @param n  Number of points to generate.
 *
 *  @par Example
 *  @snippet exemples/src/ex-tsd.cc ex_sigchirp2
 *  @image html sigchirp2.png width=600px
 *
 *  @sa sigtri(), sigsin(), sigcos(), sigexp(), sigchirp()
 */
extern ArrayXf sigchirp2(float f0, float f1, int n);


using OHConfig = tsd::OHConfig;


/** @brief Generation of an exponential signal through a harmonic oscillator.
 *
 *  <h3>Harmonic oscillator (complex output)</h3>
 *
 *  This function returns a data source, which can be called several times (generation of a continuous flow of samples, contrary to @ref sigexp(), which can only generate a fixed number of samples).
 *
 *  The quadrature oscillator is based on simple first order recursive (complex) filter:
 *  @f[
 *  z_k = z_{k-1} \cdot e^{2\pi\mathbf{i}f}
 *  @f]
 *  which enables to generate an exponential signal:
 *  @f[
 *  z_k = e^{2\pi\mathbf{i}kf}
 *  @f]
 *
 *  @note To generate an exponential of <b>fixed frequency</b>, this method will be more efficient
 *  than using the trigonometric functions of the standard library.
 *
 *  @param freq Normalized frequency of the exponential (between -0.5 et 0.5).
 *  @return Complex data source.
 *
 *  @par Example
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohc
 *  @image html ohc.png width=600px
 *
 *  @sa source_ohr(), sigexp()
 */
extern sptr<Source<cfloat, OHConfig>> source_ohc(float freq);

/** @brief Generation of a sinusoidal signal through a harmonic oscillator.
 *
 *  <h3>Harmonic oscillator (real output)</h3>
 *
 *  This function returns a data source, which can be called several times (generation of a continuous flow of samples, contrary to @ref sigsin(), which can only generate a fixed number of samples).
 *
 *  The oscillator is based on a standard quadrature one (see @ref source_ohc()), from which
 *  we use only the real part (cosinus).
 *
 *
 *  @param freq Normalized frequency of the exponential (between -0.5 et 0.5).
 *  @return Ream data source.
 *
 *  @par Example
 *  @snippet exemples/src/sdr/ex-sdr.cc ex_ohr
 *  @image html ohr.png width=600px
 *
 *  @sa source_ohc(), sigcos()
 */
extern sptr<Source<float, OHConfig>> source_ohr(float freq);



 /** @} */








}



