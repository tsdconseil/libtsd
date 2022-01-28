#pragma once

#include "tsd/tsd.hpp"
#include <vector>
#include <iosfwd>

#define ENABLE_ROOTS 1

#if ENABLE_ROOTS
#include <unsupported/Eigen/Polynomials>
#endif

namespace tsd {

/** @addtogroup filtrage
 *  @{
 */

/** @brief Polynôme (définit d'après ses racines ou ses coefficients) */
template<typename T>
struct Poly
{
  Vecteur<T> coefs;
  std::string vname;
  // Par défaut, définit d'après ses coefficients
  bool mode_racines = false;
  T mlt = 1.0f; // multiplieur si mode racine (= coefficient du monome le plus grand)



  Poly()
  {
    vname = "z";
  }

  /** @brief Retourne les coefficients du polynôme */
  Vecteur<T> get_coefs() const
  {
    Poly tmp = vers_coefs();
    return tmp.coefs;
  }

  /*Poly vers_racines() const
  {
    if(mode_racines)
      return *this;

    auto r = roots();

  }*/

  Poly vers_coefs() const
  {
    if(!mode_racines)
      return *this;

    Poly res;
    res.vname         = vname;
    res.mode_racines  = false;

    res.coefs.resize(1);
    res.coefs(0) = 1;

    for(auto i = 0; i < coefs.rows(); i++)
    {
      Poly monome;
      monome.coefs.resize(2);
      monome.coefs(0) = -coefs(i);
      monome.coefs(1) = 1;
      res = res * monome;
    }
    res = res * mlt;
    return res;
  }

  Poly<float> real() const
  {
    Poly<float> res;
    res.coefs = coefs.real();
    return res;
  }

  explicit Poly(const T &val)
  {
    mode_racines = false;
    vname = "z";
    coefs.resize(1);
    coefs(0) = val;
  }

  /** @brief Construction d'un polynôme d'après les coefficients */
  Poly(const Vecteur<T> &coefs)
  {
    mode_racines = false;
    this->coefs = coefs;
    vname = "z";
  }

  /** @brief Construction d'un polynôme d'après les racines */
  static Poly from_roots(const Vecteur<T> &racines)
  {
    Poly res;
    res.mode_racines = true;
    res.coefs = racines;
    res.mlt   = 1.0f;
    return res;
  }


  void clean(T precision)
  {
    for(auto k = 0u; k < coefs.length(); k++)
      if((coefs(k) < precision) && (coefs(k) > -precision))
        coefs(k) = 0;
  }

  /*void operator =(const Poly<T> &p)
  {
    vname = p.vname;
    coefs.resize(p.coefs.extent());
    coefs = p.coefs;
  }*/

  static const Poly<T> z;

  /** @brief Evaluation du polynôme */
  template<typename Trep>
    Trep horner(Trep val) const
  {
    auto n = coefs.rows();

    if(mode_racines)
    {
      if(n == 0)
        return Trep(mlt);

      Trep res(mlt);

      for(auto k = 0; k < n; k++)
        res = res * (val - Trep(coefs(k)));

      return res;

    }
    else
    {
      // Ex: P = a + bx + cx�
      // res = c
      // res = cx + b
      // res = cx² + bx + a

      if(n == 0)
        return Trep(0);//(Trep) 0;

      Trep res(coefs(n-1));

      for(auto k = 0; k < n - 1; k++)
        res = res * val + Trep(coefs(n-k-2));

      return res;
    }
  }

  void display_content(std::ostream &out) const
  {
    auto n = coefs.rows();

    if(mode_racines)
    {
      for(auto i = 0; i < n; i++)
      {
        if(coefs(i) == (T) 0)
          out << vname;
        else
          out << fmt::format("({} - {})", vname, coefs(i));
        if(i + 1 < n)
          out << " * ";
      }
      return;
    }

    // List des monomes non nuls;
    std::vector<std::pair<T,unsigned int>> lst;
    for(auto i = 0; i < n; i++)
    {
      if(coefs(i) != 0.0f)
      {
        auto p = std::pair<T,unsigned int>(coefs(i),i);
        lst.push_back(p);
      }
    }

    if(lst.empty())
      out << "0";

    for(auto k = 0u; k < lst.size(); k++)
    {
      T coef = lst[k].first;
      auto power = lst[k].second;

      //printf("[%f,%d]", coef, power);

      if((coef != 1.0f) || (power == 0))
        out << coef << " ";

      if(power == 1)
        out << vname;
      else if(power == 2)
        out << vname << "²";
      else if(power > 0)
        out << vname << "^" << power;

      if(k + 1 < lst.size())
      {
        out << " + ";
      }
    }
  }

  Poly<T> operator ^(unsigned int e) const
  {
    if(e == 0)
    {
      return Poly<T>::one();
    }
    if(mode_racines)
    {
      Poly<T> res = *this;
      for(auto k = 1u; k < e; k++)
        // Duplique les racines
        res.coefs = vconcat(res.coefs, coefs).eval();
      res.mlt = std::pow(mlt, e);
      return res;
    }
    else
    {
      Poly<T> res;
      res.vname = vname;
      res.coefs.setOnes(1);

      // TODO: algo indien
      for(auto k = 0u; k < e; k++)
        res = res * *this;
      return res;
    }
  }

  Poly<T> operator <<(int i) const
  {
    Poly<T> res;
    res.vname = vname;
    res.coefs.resize(coefs.rows() + i);
    res.coefs.setZero();
    for(auto k = 0; k < coefs.rows(); k++)
      res.coefs(i+k) = coefs(k);
    return res;
  }

  /*Poly<T> operator *(const float &s) const
  {
    auto res = *this;
    res.coefs *= s;
    return res;
  }*/

  Poly<T> operator *(const T &s) const
  {
    if(mode_racines)
    {
      Poly<T> res = *this;
      res.mlt *= s;
      return res;
    }
    else
    {
      auto res = *this;
      res.coefs *= s;
      return res;
    }
  }



  Poly<T> operator *(const Poly<T> &s) const
  {
    Poly<T> res;
    res.vname = vname;
    res.mode_racines = mode_racines;

    auto s1 = *this, s2 = s;

    if(s1.mode_racines && !s2.mode_racines)
    {
      s1 = s1.vers_coefs();
    }
    else if(!s1.mode_racines && s2.mode_racines)
    {
      s2 = s2.vers_coefs();
    }

    tsd_assert(s1.mode_racines == s2.mode_racines);

    if(s1.mode_racines)
    {
      res.coefs = vconcat(s1.coefs, s2.coefs);
      res.mlt = s1.mlt * s2.mlt;
    }
    else
    {
      for(auto k = 0; k < s2.coefs.rows(); k++)
      {
        Poly<T> prod = s1;
        // Shift prod and multiply by coef
        prod.coefs.setZero(k + s1.coefs.rows());
        for(auto i = 0; i < s1.coefs.rows(); i++)
          prod.coefs(i+k) = s1.coefs(i) * s2.coefs(k);
        res = res + prod;
      }
    }
    return res;
  }

  Poly<T> operator -(const T &s) const
  {
    return *this + (-s);
  }

  Poly<T> operator +(const T &s) const
  {
    auto res = vers_coefs();
    if(res.coefs.rows() == 0)
      res.coefs.resize(1);
    res.coefs(0) += s;
    return res;
  }

  Poly<T> operator +(const Poly<T> &s) const
  {
    auto s1 = vers_coefs(), s2 = s.vers_coefs();

    Poly<T> res;
    auto n1 = s1.coefs.rows(), n2 = s2.coefs.rows();
    res.vname = vname;

    if(n1 == n2)
      res.coefs = s1.coefs + s2.coefs;
    else if(n1 > n2)
    {
      res.coefs = s1.coefs;
      res.coefs.segment(0,n2) += s2.coefs;
    }
    else
    {
      res.coefs = s2.coefs;
      res.coefs.segment(0,n1) += s1.coefs;
    }
    return res;
  }

  Poly<T> operator -() const
  {
    auto res = *this;
    if(mode_racines)
    {
      res.mlt = -res.mlt;
    }
    else
    {
      res.coefs = -coefs;
    }
    return res;
  }

  Poly<T> operator -(const Poly<T> &s) const
  {
    return *this + (-s);
  }


  static Poly<T> basic()
  {
    Poly<T> res;
    res.coefs.resize(2);
    res.coefs(0) = 0;
    res.coefs(1) = 1;
    return res;
  }

  static Poly<T> one()
  {
    Poly<T> res;
    res.coefs.resize(1);
    res.coefs(0) = 1;
    return res;
  }

  T norm() const
  {
    return coefs.square().sum();
  }

  /*
  // Que faire si T est déjà complexe ?
  template<class S = T, std::enable_if_t<is_complex<S>{}>* = nullptr>
  Vecteur<T> roots() const
  {
    Eigen::PolynomialSolver<T, Eigen::Dynamic> solver;

    if(coefs.rows() <= 1)
      return Vecteur<T>(0);

    Vecteur<T> v = coefs;

    while(v(v.rows()-1) == 0.0f)
      v = v.head(v.rows()-1).eval();

    solver.compute(v);
    auto &r = solver.roots();
    return Vecteur<T>(r);
  }*/

  //template<class S = T, std::enable_if_t<!is_complex<S>{}>* = nullptr>
  /** @brief Retourne les racines du polynôme */
  /*Vecteur<std::complex<float>>*/ArrayXcf roots() const
  {
    if(mode_racines)
    {
      return coefs;
    }
    if constexpr(est_complexe<T>())
    {
      echec("Pas possible de calculer les racines d'un polynôme complexe avec cette version de Eigen.");
      return ArrayXcf();
    }
    else
    {
#     if ENABLE_ROOTS
      Eigen::PolynomialSolver<T, Eigen::Dynamic> solver;

      if(coefs.rows() <= 1)
        return Vecteur<std::complex<float>>();

      Vecteur<T> v = coefs;

      // Obligé, sinon la fonction Eigen se bloque (exemple : filtre Gaussien)
      while((v.rows() > 1) && (std::abs(v(v.rows()-1)) < 1e-8)) //== 0.0f))
        v = v.head(v.rows()-1).eval();

      // Polynome nul
      if((v.rows() == 1) && (v(0) == 0.0f))
        return Vecteur<std::complex<float>>();

      solver.compute(v);
      auto &r = solver.roots();
      return Vecteur<std::complex<float>>(r);
#     else
      erreur("TODO : FRat::roots()");
      return ArrayXcf(0);
#     endif
    }
  }


};


template<typename T>
  Poly<T> operator *(const T &s, const Poly<T> &p)
{
  return p * s;
}

template<typename T>
  Poly<T> operator +(const T &s, const Poly<T> &p)
{
  return p + s;
}

template<typename T>
  Poly<T> operator -(const T &s, const Poly<T> &p)
{
  return (-p) + s;
}

template<typename T>
  const Poly<T> Poly<T>::z = Poly<T>::basic();

template<typename T>
std::ostream &operator<< (std::ostream &out, const Poly<T> &m)
{
  Poly<T> cp = m;
  //cp.simplify();
  cp.display_content(out);
  return out;
}

template<typename T>
Poly<T> real(const Poly<std::complex<T>> &p)
{
  Poly<T> res;
  res.coefs.resize(p.coefs.rows());
  for(auto i = 0; i < p.coefs.rows(); i++)
    res.coefs(i) = p.coefs(i).real();
  return res;
}


/** @brief Fraction rationnelle polynomiale, avec coefficients réels ou complexes.
 *
 * La fraction peut être représentée sous forme coefficients :
 * @f[
 * H(z) = \frac{a_{K-1} + a_{K-2} z + ...  a_0 z^{K-1}}{b_{L-1} + b_{L-2} z + ...  b_0 z^{L-1}}
 * @f]
 * ou racines :
 * @f[
 * H(z) = G\cdot\frac{\prod z - z_i}{\prod z - p_i}
 * @f]
 */
template<typename T = float>
class FRat
{
public:

  /** @brief Constructeur, d'après un vecteur de coefficients (filtre RIF). */
  FRat(Eigen::Array<T,-1,1> x)
  {
    *this = rif(x);
  }

  /*FRat(const Eigen::Array<T,-1,1> &x)
  {
    *this = rif(x);
  }*/

  explicit FRat(const T &a)
  {
    numer = Poly<T>::one() * a;
    denom = Poly<T>::one();
  }


  /** @brief Constructeur par défaut */
  FRat()
  {
    numer = Poly<T>::one();
    denom = Poly<T>::one();
  }



  /** @brief Fraction rationnelle d'un filtre RIF.
   *
   * <h3>Fraction rationnelle d'un filtre RIF</h3>
   *
   *  Calculée à partir des coefficients @f$a_k@f$ de l'équation aux différences :
   *
   *  @f[
   *    y_n =  a_0 x_n + a_1 x_{n-1} + ... + a_{K-1} x_{n-K+1}
   *  @f]
   *
   *  La fonction de transfert résultante est :
   *
   *  @f[
   * H(z) = \frac{a_{K-1} + a_{K-2} z + ...  a_0 z^{K-1}}{z^{K-1}}
   *  @f]
   *
   *  @param a Liste des coefficients
   */
  static FRat<T> rif(const Vecteur<T> &a)
  {
    FRat<T> res;
    auto K = a.rows();
    if(K == 0)
    {
      msg_erreur("Filtre RIF : aucun coefficient.");
      return res;
    }

    res.numer = Poly<T>(a.reverse());
    Vecteur<T> c = Vecteur<T>::Zero(K);
    c(K - 1) = 1;
    res.denom = Poly<T>(c);
    return res;
  }

  /** @brief Calcul de @f$H(z^{-1})@f$
   *
   *  <h3>Calcul de @f$H(z^{-1})@f$</h3>
   *
   *  Cette fonction remplace @f$z@f$ par @f$z^{-1}@f$.
   */
  FRat<T> eval_inv_z() const
  {
    FRat<T> res;

    auto a = numer.vers_coefs().coefs;
    auto b = denom.vers_coefs().coefs;

    int nn = a.rows();
    int nd = b.rows();

    res.numer = Poly<T>(a.reverse());
    res.denom = Poly<T>(b.reverse());

    if(nn > nd)
      res.denom = res.denom * Poly<T>::z ^ (nn - nd);
    else if(nd > nn)
      res.numer = res.numer * Poly<T>::z ^ (nd - nn);

    return res;
  }


  /** @brief Fraction rationnelle d'un filtre RII.
   *
   * <h3>Fraction rationnelle d'un filtre RII</h3>
   *
   *  Calculée à partir des coefficients @f$a_k@f$ et @f$b_k@f$ de l'équation aux différences :
   *
   *  @f[
   *    b_0 y_n + b_1 y_{n-1} + \dots + b_{L-1} y_{n-L+1} =  a_0 x_n + a_1 x_{n-1} + ... + a_{K-1} x_{n-K+1}
   *  @f]
   *
   *  La fonction de transfert résultante est :
   *
   *  @f[
   * H(z) = \frac{a_{K-1} + a_{K-2} z + ...  a_0 z^{K-1}}{b_{L-1} + b_{L-2} z + ...  b_0 z^{L-1}} \cdot \frac{z^{L-1}}{z^{K-1}}
   *  @f]
   *
   *  @param a Liste des coefficients (numérateur)
   *  @param b Liste des coefficients (dénominateur)
   */
  static FRat<T> rii(const Vecteur<T> &a, const Vecteur<T> &b)
  {
    FRat<T> res;

    if((a.rows() == 0) || (b.rows() == 0))
    {
      msg_erreur("Filtre RII : aucun coefficient.");
      return res;
    }

    res.numer = Poly<T>(a);
    res.denom = Poly<T>(b);

    return res.eval_inv_z();
  }

  static FRat<T> un()
  {
    FRat<T> res;
    res.numer = Poly<T>::one();
    res.denom = Poly<T>::one();
    return res;
  }


  /** @brief Renvoie le monôme @f$z@f$
   *  <h3>Monôme @f$z@f$</h3>
   *
   *  Cette fonction peut être utile pour construire de manière algébrique une fonction de transfert.
   *
   *  @par Exemple
   *  @code
   *  auto z = FRat<float>::z();
   *
   *  auto H = (z-1) / (z.pow(2) - 0.5 * z + 2);
   *  @endcode
   *
   */
  static FRat<T> z()
  {
    FRat<T> res;
    res.numer = Poly<T>::z;
    res.denom = Poly<T>::one();
    return res;
  }

  static FRat<T> z_power(int n)
  {
    FRat<T> res;
    if(n == 0)
    {
      return un();
    }
    else if(n > 0)
    {
      res.numer.coefs = Vecteur<T>::Zero(n+1);
      res.numer.coefs(n) = 1;
      res.denom.coefs = Vecteur<T>::Ones(1);
    }
    else if(n < 0)
    {
      res.denom.coefs = Vecteur<T>::Zero(-n+1);
      res.denom.coefs(-n) = 1;
      res.numer.coefs = Vecteur<T>::Ones(1);
    }
    return res;
  }


  /** @brief Evaluation d'une fraction rationnelle.
   *
   *  <h3>Evaluation d'une fraction rationnelle</h3>
   *  Calcul de @f$H(x)@f$, @f$x@f$ pouvant être de type quelquonque.
   *
   *  @par Exemple
   *  @code
   *  // En supposant que H est de type FRat<float>
   *  float y = H.horner(5);
   *  @endcode
   *
   */
  template<typename Trep>
    Trep horner(Trep x) const
  {
    return numer.horner(x) / denom.horner(x);
  }

  FRat operator *(const T &s) const
  {
    FRat res;
    res.numer = numer * s;
    res.denom = denom;
    return res;
  }

  FRat<T> operator +(const T &s) const
  {
    FRat<T> res;
    res.numer = numer + denom * s;
    res.denom = denom;
    return res;
  }

  FRat<T> operator -(const T &s) const
  {
    FRat<T> res;
    res.numer = numer - denom * s;
    res.denom = denom;
    return res;
  }

  /** @brief Addition de deux fractions rationnelles. */
  FRat<T> operator +(const FRat<T> &s) const
  {
    FRat<T> res;
    res.numer = numer * s.denom + s.numer * denom;
    res.denom = denom * s.denom;
    return res;
  }

  /** @brief Soustraction de deux fractions rationnelles. */
  FRat<T> operator -(const FRat<T> &s) const
  {
    FRat<T> res;
    res.numer = numer * s.denom - s.numer * denom;
    res.denom = denom * s.denom;
    return res;
  }

  /** @brief Produit de deux fractions rationnelles. */
  FRat<T> operator *(const FRat<T> &s) const
  {
    FRat res;
    res.numer = numer * s.numer;
    res.denom = denom * s.denom;
    return res;
  }

  /** @brief Division de deux fractions rationnelles. */
  FRat<T> operator /(const FRat<T> &s) const
  {
    FRat<T> res;
    res.numer = numer * s.denom;
    res.denom = denom * s.numer;
    return res;
  }

  FRat<T> operator-() const
  {
    FRat<T> res;
    res.numer = -numer ;
    res.denom = denom;
    return res;
  }


  /** @brief Calcul de @f$H^N(z)@f$. */
  FRat<T> pow(unsigned int N) const
  {
    if(N == 0)
      return FRat<T>::un();
    else if(N == 1)
      return *this;
    return pow(N-1) * *this;
  }

  /** @brief Calcul de @f$1/H(z)@f$ */
  FRat<T> inv() const
  {
    FRat<T> res;
    res.numer = denom;
    res.denom = numer;
    return res;
  }

  bool est_fir() const
  {
    auto Nd = denom.coefs.rows();
    return denom.coefs.head(Nd-1).isZero() && (denom.coefs(Nd-1) == 1.0f);
  }

  void afficher(std::ostream &out) const
  {
    out << "(";
    numer.display_content(out);
    out << ") / (";
    denom.display_content(out);
    out << ")";
  }

  Poly<T> numer, denom;
};


// Quand le scalaire vient en premier
template<typename T1, typename T2>
  FRat<T2> operator /(const T1 &v, const FRat<T2> &s)
{
  return s.inv() * v;
}

template<typename T1, typename T2>
  FRat<T2> operator *(const T1 &v, const FRat<T2> &s)
{
  return s * v;
}

template<typename T1, typename T2>
  FRat<T2> operator +(const T1 &v, const FRat<T2> &s)
{
  return s + v;
}

template<typename T1, typename T2>
  FRat<T2> operator -(const T1 &v, const FRat<T2> &s)
{
  return s - v;
}

/** @} */


std::ostream& operator<<(std::ostream& os, const Poly<float> &p);
std::ostream& operator<<(std::ostream& os, const FRat<float> &f);


}


