﻿#pragma once

#include "tsd/tsd.hpp"
#include <iosfwd>



namespace tsd {

/** @addtogroup filtrage-poly
 *  @{
 */

/** @brief Polynôme (définit d'après ses racines ou ses coefficients) */
template<typename T>
struct Poly
{
  VecT<T> coefs;
  string vname = "z";
  // Par défaut, définit d'après ses coefficients
  bouléen mode_racines = non;
  T mlt = 1.0f; // multiplieur si mode racine (= coefficient du monome le plus grand)

  Poly()
  {
  }

  /** @brief Construction d'un polynôme d'après les coefficients */
  Poly(const Vecteur<T> &coefs)
  {
    mode_racines = non;
    this->coefs = coefs;
    vname = "z";
  }

  /** @brief Construction d'un polynôme d'après les racines */
  static Poly from_roots(const Vecteur<T> &racines)
  {
    Poly res;
    res.mode_racines  = oui;
    res.coefs         = racines.clone();
    res.mlt           = 1.0f;
    retourne res;
  }

  explicit Poly(const T &val)
  {
    mode_racines = non;
    vname = "z";
    coefs.resize(1);
    coefs(0) = val;
  }

  /** @brief Retourne les coefficients du polynôme */
  Vecteur<T> get_coefs() const
  {
    retourne vers_coefs().coefs;
  }

  entier degré() const
  {
    retourne coefs.rows() - 1;
  }

  // Nombre de coefficients successifs nuls, en commençant par 0
  entier nb_coefs_nuls() const
  {
    entier r;
    pour(r = 0; r < coefs.rows(); r++)
    {
      si(coefs(r) != (T) 0)
        break;
    }
    retourne r;
  }

  /*Poly &operator =(const Poly<T> &src)
  {
    coefs        = src.coefs.clone();
    mlt          = src.mlt;
    mode_racines = src.mode_racines;
    retourne *this;
  }*/

  // TODO template: cfloat ou double...
  Poly<cfloat> factorise() const
  {
    retourne Poly<cfloat>::from_roots(roots());
  }

  Poly vers_coefs() const
  {
    si(!mode_racines)
      retourne *this;

    Poly res;
    res.vname         = vname;
    res.mode_racines  = non;

    res.coefs.resize(1);
    res.coefs(0) = 1;

    pour(auto i = 0; i < coefs.rows(); i++)
    {
      soit monome = Poly(VecT<T>::valeurs({-coefs(i), 1}));
      res *= monome;
      /*Poly monome;
      monome.coefs.resize(2);
      monome.coefs(0) = -coefs(i);
      monome.coefs(1) = 1;
      res = res * monome;*/
    }
    res *= mlt;
    //res = res * mlt;
    retourne res;
  }

  Poly<float> real() const
  {
    retourne {tsd::real(coefs)};
    /*Poly<float> res;
    res.coefs = tsd::real(coefs);
    retourne res;*/
  }




  void clean(T precision)
  {
    pour(auto k = 0u; k < coefs.length(); k++)
      si((coefs(k) < precision) && (coefs(k) > -precision))
        coefs(k) = 0;
  }

  /*void operator =(const Poly<T> &p)
  {
    vname = p.vname;
    coefs.resize(p.coefs.extent());
    coefs = p.coefs;
  }*/

  static const Poly<T> z;


  /** @brief Evaluation sur un ensemble de points */
  Vecf éval(const Vecf &x)
  {
    retourne x.unaryExpr([&](float xf){return horner(xf);});
  }


  /** @brief Evaluation du polynôme */
  template<typename Trep>
    Trep horner(Trep val) const
  {
    soit n = coefs.rows();

    si(mode_racines)
    {
      si(n == 0)
        retourne Trep(mlt);

      Trep res(mlt);

      pour(auto k = 0; k < n; k++)
        res = res * (val - Trep(coefs(k)));

      retourne res;

    }
    sinon
    {
      // Ex: P = a + bx + cx�
      // res = c
      // res = cx + b
      // res = cx² + bx + a

      si(n == 0)
        retourne Trep(0);//(Trep) 0;

      Trep res(coefs(n-1));

      pour(auto k = 0; k < n - 1; k++)
        res = res * val + Trep(coefs(n-k-2));

      retourne res;
    }
  }

  void display_content(std::ostream &out) const
  {
    soit n = coefs.rows();

    si(mode_racines)
    {
      pour(auto i = 0; i < n; i++)
      {
        si(coefs(i) == (T) 0)
          out << vname;
        sinon
          out << sformat("({} - {})", vname, coefs(i));
        si(i + 1 < n)
          out << " * ";
      }
      retourne;
    }

    // List des monomes non nuls;
    vector<std::pair<T,unsigned int>> lst;
    pour(auto i = 0; i < n; i++)
      si(coefs(i) != 0.0f)
        lst.push_back({coefs(i),i});

    si(lst.empty())
      out << "0";

    pour(auto k = 0u; k < lst.size(); k++)
    {
      soit coef  = lst[k].first;
      soit power = lst[k].second;

      si((coef != 1.0f) || (power == 0))
        out << coef << " ";

      si(power == 1)
        out << vname;
      sinon si(power == 2)
        out << vname << "²";
      sinon si(power > 0)
        out << vname << "^" << power;

      si(k + 1 < lst.size())
        out << " + ";
    }
  }

  Poly<T> operator ^(unsigned int e) const
  {
    si(e == 0)
      retourne Poly<T>::one();

    si(mode_racines)
    {
      soit res = *this;
      pour(auto k = 1u; k < e; k++)
        // Duplique les racines
        res.coefs = vconcat(res.coefs, coefs);
      res.mlt = pow(mlt, e);
      retourne res;
    }
    sinon
    {
      Poly<T> res;
      res.vname = vname;
      res.coefs = VecT<T>::ones(1);

      // TODO: algo indien
      Pour(auto k = 0u; k < e; k++)
        res = res * *this;
      retourne res;
    }
  }

  Poly<T> operator <<(entier i) const
  {
    Poly<T> res;
    res.vname = vname;
    res.coefs = VecT<T>::zeros(coefs.rows() + i);
    pour(auto k = 0; k < coefs.rows(); k++)
      res.coefs(i+k) = coefs(k);
    retourne res;
  }

  Poly<T> operator *(const T &s) const
  {
    si(mode_racines)
    {
      Poly<T> res;
      res.mode_racines  = oui;
      res.coefs         = coefs.clone();
      res.mlt           = mlt * s;
      res.vname         = vname;
      retourne res;
    }
    sinon
    {
      Poly<T> res;
      res.coefs         = coefs * s;
      res.mode_racines  = non;
      res.vname         = vname;
      retourne res;
    }
  }

  Poly<T> &operator *=(const T &x)
  {
    *this = *this * x;
    retourne *this;
  }

  Poly<T> &operator +=(const T &x)
  {
    *this = *this + x;
    retourne *this;
  }

  Poly<T> &operator *=(const Poly<T> &x)
  {
    *this = *this * x;
    retourne *this;
  }

  Poly<T> &operator +=(const Poly<T> &x)
  {
    *this = *this + x;
    retourne *this;
  }

  Poly<T> operator *(const Poly<T> &s) const
  {
    Poly<T> res;
    res.vname = vname;
    res.mode_racines = mode_racines;

    Poly<T> s1 = *this, s2 = s;

    si(s1.mode_racines && !s2.mode_racines)
      s1 = s1.vers_coefs();
    sinon si(!s1.mode_racines && s2.mode_racines)
      s2 = s2.vers_coefs();

    assertion(s1.mode_racines == s2.mode_racines);

    si(s1.mode_racines)
    {
      res.coefs = vconcat(s1.coefs, s2.coefs);
      res.mlt = s1.mlt * s2.mlt;
    }
    sinon
    {
      soit n1 = s1.coefs.rows(),
           n2 = s2.coefs.rows();

      pour(auto k = 0; k < n2; k++)
      {
        Poly<T> prod;
        // Shift prod and multiply by coef
        prod.coefs.setZero(k + n1);

        assertion(prod.coefs.rows() == (k + n1));

        pour(auto i = 0; i < n1; i++)
          prod.coefs(i+k) = s1.coefs(i) * s2.coefs(k);

        res = res + prod;
      }
    }
    retourne res;
  }

  Poly<T> operator -(const T &s) const
  {
    return *this + (-s);
  }

  Poly<T> operator +(const T &s) const
  {
    soit res = vers_coefs();
    si(res.coefs.rows() == 0)
      res.coefs = VecT<T>::zeros(1);
    res.coefs(0) += s;
    retourne res;
  }

  Poly<T> operator +(const Poly<T> &s) const
  {
    soit s1 = vers_coefs(),
         s2 = s.vers_coefs();

    Poly<T> res;
    soit n1 = s1.coefs.rows(),
         n2 = s2.coefs.rows();
    res.vname = vname;

    si(n1 == n2)
      res.coefs = s1.coefs + s2.coefs;
    sinon si(n1 > n2)
    {
      res.coefs = s1.coefs.clone();
      res.coefs.head(n2) += s2.coefs;
    }
    sinon
    {
      res.coefs = s2.coefs.clone();
      res.coefs.head(n1) += s1.coefs;
    }
    retourne res;
  }

  Poly<T> operator -() const
  {
    soit res = *this;
    si(mode_racines)
      res.mlt = -res.mlt;
    sinon
      res.coefs = -coefs;
    retourne res;
  }

  Poly<T> operator -(const Poly<T> &s) const
  {
    retourne *this + (-s);
  }


  static Poly<T> basic()
  {
    /*Poly<T> res;
    res.coefs.resize(2);
    res.coefs(0) = 0;
    res.coefs(1) = 1;
    retourne res;*/
    retourne {VecT<T>::valeurs({(T) 0, (T) 1})};
  }

  static Poly<T> one()
  {
    Poly<T> res;
    res.coefs.resize(1);
    res.coefs(0) = 1;
    retourne res;
  }

  T norm() const
  {
    retourne coefs.square().somme();
  }

  /** @brief Retourne les racines du polynôme */

  Veccf roots() const;
};


/** @cond private  */
template<typename T>
  Poly<T> operator *(const T &s, const Poly<T> &p)
{
  retourne p * s;
}

template<typename T>
  Poly<T> operator +(const T &s, const Poly<T> &p)
{
  retourne p + s;
}

template<typename T>
  Poly<T> operator -(const T &s, const Poly<T> &p)
{
  retourne (-p) + s;
}

template<typename T>
  const Poly<T> Poly<T>::z = Poly<T>::basic();

template<typename T>
std::ostream &operator<< (std::ostream &out, const Poly<T> &m)
{
  Poly<T> cp = m;
  //cp.simplify();
  cp.display_content(out);
  retourne out;
}

template<typename T>
Poly<T> real(const Poly<std::complex<T>> &p)
{
  Poly<T> res;
  soit n = p.coefs.rows();
  res.coefs.resize(n);
  pour(auto i = 0; i < n; i++)
    res.coefs(i) = p.coefs(i).real();
  retourne res;
}
/** @endcond  */

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
struct FRat
{
  Poly<T> numer, denom;

  using VT = Vecteur<T>;
  using PT = Poly<T>;

  /** @brief Constructeur, d'après un vecteur de coefficients (filtre RIF). */
  FRat(const VT &x)
  {
    *this = rif(x);
  }

  explicit FRat(const T &a)
  {
    numer = PT::one() * a;
    denom = PT::one();
  }

  /** @brief Constructeur par défaut */
  FRat()
  {
    numer = PT::one();
    denom = PT::one();
  }

  FRat(const PT &numer, const PT &denom)
  {
    this->numer = numer;
    this->denom = denom;
  }

  static FRat<T> un()
  {
    retourne {PT::one(), PT::one()};
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
    retourne {PT::z, PT::one()};
  }

  static FRat<T> zinv()
  {
    retourne {PT::one(), PT::z};
  }


  tuple<Veccf,Veccf> roots() const
  {
    retourne {numer.roots(), denom.roots()};
  }

  tuple<entier,entier> degrés() const
  {
    retourne {numer.degré(), denom.degré()};
  }

  /** @brief Développement des numérateur et dénominateur sous la forme de polynômes définis par leurs coefficients. */
  FRat<T> developpe() const
  {
    retourne {numer.vers_coefs(), denom.vers_coefs()};
  }

  // TODO: cfloat ou cdouble...
  FRat<cfloat> factorise() const
  {
    FRat<cfloat> f{numer.factorise(), denom.factorise()};
    f.simplifier();
    retourne f;
  }


  //FRat<T> simplifier() const
  //{
    // TODO : supprimer facteurs communs...
    //retourne *this;
  //}

  // (a/b).real = a.r * (1/b).r - a.i * (1/b).i
  // @todo : renommer, et utiliser une autre fonction pour la partie réelle de la fraction
  /** @brief Extraction de la partie réelle des numérateurs et dénominateurs (attention ce n'est pas la partie réelle de la fraction !). */
  auto real() const
  {
    if constexpr(est_complexe<T>())
    {
      retourne FRat<typename T::value_type>{numer.real(), denom.real()};
    }
    sinon
      retourne *this;
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
  static FRat<T> rif(const VT &a)
  {
    soit K = a.rows();
    si(K == 0)
    {
      msg_erreur("Filtre RIF : aucun coefficient.");
      retourne {};
    }

    soit c = VT::zeros(K);
    c(K - 1) = 1;
    retourne {PT(a.reverse()), PT(c)};
  }

  VT coefs_rif() const
  {
    si(!est_rif())
      échec("FRat::coefs_rif(): ce filtre n'est pas un filtre RIF.");
    retourne numer.coefs.reverse();
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

    soit a = numer.vers_coefs().coefs,
         b = denom.vers_coefs().coefs;

    soit nn = a.rows(),
         nd = b.rows();

    res.numer = PT(a.reverse());
    res.denom = PT(b.reverse());

    Si(nn > nd)
      res.denom = res.denom * PT::z ^ (nn - nd);
    sinon si(nd > nn)
      res.numer = res.numer * PT::z ^ (nd - nn);

    retourne res;
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

    si((a.rows() == 0) || (b.rows() == 0))
    {
      msg_erreur("Filtre RII : aucun coefficient.");
      retourne {};
    }

    res.numer = PT(a);
    res.denom = PT(b);

    retourne res.eval_inv_z();
  }



  static FRat<T> z_power(entier n)
  {
    FRat<T> res;
    si(n == 0)
      retourne un();
    sinon si(n > 0)
    {
      res.numer.coefs = VT::zeros(n+1);
      res.numer.coefs(n) = 1;
      res.denom.coefs = VT::ones(1);
    }
    sinon
    {
      res.denom.coefs = VT::zeros(n+1);
      res.denom.coefs(n) = 1;
      res.numer.coefs = VT::ones(1);
    }
    retourne res;
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
    retourne numer.horner(x) / denom.horner(x);
  }

  FRat operator *(const T &s) const
  {
    retourne {numer * s, denom};
  }

  FRat<T> operator +(const T &s) const
  {
    retourne {numer + denom * s, denom};
  }

  FRat<T> operator -(const T &s) const
  {
    retourne (*this) + (-s);
  }

  /** @brief Addition de deux fractions rationnelles. */
  FRat<T> operator +(const FRat<T> &s) const
  {
    retourne {numer * s.denom + s.numer * denom, denom * s.denom};
  }

  /** @brief Soustraction de deux fractions rationnelles. */
  FRat<T> operator -(const FRat<T> &s) const
  {
    retourne (*this) + (-s);
  }

  /** @brief Produit de deux fractions rationnelles. */
  FRat<T> operator *(const FRat<T> &s) const
  {
    retourne {numer * s.numer, denom * s.denom};
  }

  // Simplification si facteur z^m commun
  void simplifier()
  {
    if(numer.mode_racines)
      return; // A faire
    soit nr = numer.coefs.rows(),
         nd = denom.coefs.rows();
    soit nr1 = numer.nb_coefs_nuls(),
         nd1 = denom.nb_coefs_nuls();

    soit m = min(nr1, nd1);
    si(m > 0)
    {
      numer.coefs = numer.coefs.tail(nr - m);
      denom.coefs = denom.coefs.tail(nd - m);
    }
  }

  /** @brief Division de deux fractions rationnelles. */
  FRat<T> operator /(const FRat<T> &s) const
  {
    FRat<T> res = {numer * s.denom, denom * s.numer};
    res.simplifier();
    retourne res;
  }

  FRat<T> operator-() const
  {
    retourne {-numer, denom};
  }

  FRat<T> &operator *=(const T &x)
  {
    *this = *this * x;
    retourne *this;
  }

  FRat<T> &operator +=(const T &x)
  {
    *this = *this + x;
    retourne *this;
  }

  FRat<T> &operator -=(const T &x)
  {
    *this = *this - x;
    retourne *this;
  }

  FRat<T> &operator /=(const T &x)
  {
    *this = *this / x;
    retourne *this;
  }


  /** @brief Calcul de @f$H^N(z)@f$. */
  FRat<T> pow(unsigned int n) const
  {
    si(n == 0)
      retourne FRat<T>::un();
    sinon si(n == 1)
      retourne *this;
    retourne pow(n-1) * *this;
  }

  /** @brief Calcul de @f$1/H(z)@f$ */
  FRat<T> inv() const
  {
    retourne {denom, numer};
  }

  bouléen est_rif() const
  {
    soit Nd = denom.coefs.rows();
    retourne denom.coefs.head(Nd-1).est_nul() && (denom.coefs(Nd-1) == 1.0f);
  }

  void afficher(std::ostream &out) const
  {
    out << "(";
    numer.display_content(out);
    out << ") / (";
    denom.display_content(out);
    out << ")";
  }

};

/** @} */

/** @cond private  */
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
  return -s + v;
}

std::ostream& operator<<(std::ostream& os, const Poly<float> &p);
std::ostream& operator<<(std::ostream& os, const FRat<float> &f);
std::ostream& operator<<(std::ostream& os, const Poly<cfloat> &p);
std::ostream& operator<<(std::ostream& os, const FRat<cfloat> &f);
/** @endcond */




}


ostream_formater(tsd::Poly<float>)
ostream_formater(tsd::FRat<float>)
ostream_formater(tsd::Poly<tsd::cfloat>)
ostream_formater(tsd::FRat<tsd::cfloat>)

