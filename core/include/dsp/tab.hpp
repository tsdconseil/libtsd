#ifndef DSP_TAB_H
#define DSP_TAB_H


#include "tsd/tableau.hpp"


namespace dsp {




using tsd::cfloat;
using tsd::cdouble;

template<typename T>
  constexpr bool is_complex() { retourne tsd::est_complex_t<T>(); }

using tsd::bool_t;

// Tableau de données évaluées, avec type spécifique, et nombre de dimensions spécifique
template<typename T, entier ndims>
struct TabT
{
  tsd::TabT<T, ndims> fr;


  TabT()
  {
  }


  TabT(tsd::Tab &&t): fr(t)
  {
  }

  TabT(const tsd::Tab &t): fr(t)
  {
  }

  /*TabT(tsd::TabT<T, ndims> fr)
  {
    this->fr = fr;
  }*/

  TabT(const tsd::TabT<T, ndims> &fr_): fr(fr_)
  {
    //msg("En::Tabt(const &fr) : just called TabT(const &).");
  }

  TabT(tsd::TabT<T, ndims> &&fr_): fr(std::move(fr_))
  {
    //msg("En::Tabt(&&fr) : just called TabT(&&).");
  }


  TabT(const TabT<T, ndims> &en): fr(en.fr)
  {
    //msg("En::Tabt(const &en) : just called TabT(const &).");
  }

  TabT(TabT<T, ndims> &&en): fr(std::move(en.fr))
  {
    //msg("En::Tabt(&&en) : just called TabT(move).");
  }

  TabT &operator =(TabT &&t)
  {
    //msg("TabEn:: operator =(&&) -> fr...");
    fr = std::move(t.fr);
    //msg("TabEn:: operator =(&&) -> fr ok.");
    return *this;
  }

  TabT &operator =(const TabT &t)
  {
    //msg("TabEn:: operator =() -> fr...");
    fr = t.fr;
    //msg("TabEn:: operator =() -> fr ok.");
    return *this;
  }



  template<typename T2, entier ndims2 = ndims>
    TabT<T2, ndims2> as() const
  {
    return fr.template as<T2,ndims2>();
  }


  template<typename T2, entier ndims2 = ndims>
    TabT<T2, ndims2> as_no_conv() const
  {
    return fr.template as_no_conv<T2,ndims2>();
  }

  entier rows()       const
  {
    return fr.rows();
  }

  entier cols()       const
  {
    return fr.cols();
  }

  T prod() const
  {
    return fr.produit();
  }

  TabT<T,ndims> cwiseSign() const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.cwiseSign();
  }


  TabT<T,ndims> cwiseMax(T x) const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.cwiseMax(x);
  }

  TabT<T,ndims> cwiseMax(const TabT<T,ndims> &x) const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.cwiseMax(x);
  }

  TabT<T,ndims> cwiseMin(const TabT<T,ndims> &x) const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.cwiseMin(x);
  }

  bool all_true() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    return fr.tous_vrai();
  }

  bool all_false() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    return fr.tous_faux();
  }

  bool is_any_true() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    return fr.au_moins_un_vrai();
  }

  bool is_any_false() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    return fr.au_moins_un_faux();
  }

  template<typename T2>
  TabT<T2,ndims> select(const TabT<T2,ndims> &x1, const TabT<T2,ndims> &x2) const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    return fr.select(x1, x2);
  }

  TabT<T,ndims> cwiseMin(T x) const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.cwiseMin(x);
  }

  entier index_max() const NOECLIPSE(requires(ndims == 1))
  {
    return fr.index_max();
  }

  entier index_min() const NOECLIPSE(requires(ndims == 1))
  {
    return fr.index_min();
  }

  T max_value() const
  {
    return fr.valeur_max();
  }

  T min_value() const
  {
    return fr.valeur_min();
  }

  tuple<T, T> minmax_values() const NOECLIPSE(requires(!is_complex<T>()))
  {
    return fr.valeurs_minmax();
  }

  tuple<entier, entier> max_index() const NOECLIPSE(requires((ndims == 2) && !is_complex<T>()))
  {
    return fr.index_max();
  }

  tuple<entier, entier> min_index() const NOECLIPSE(requires((ndims == 2) && !is_complex<T>()))
  {
    return fr.index_min();
  }


  tuple<T, entier> max() const NOECLIPSE(requires((ndims == 1) && !is_complex<T>()))
  {
    return fr.max();
  }

  tuple<T, entier, entier> max() const NOECLIPSE(requires((ndims == 2) && !is_complex<T>()))
  {
    return fr.max();
  }

  tuple<T, entier> min() const NOECLIPSE(requires((ndims == 1) && !is_complex<T>()))
  {
    return fr.min();
  }

  tuple<T, entier, entier> min() const NOECLIPSE(requires((ndims == 2) && !is_complex<T>()))
  {
    return fr.min();
  }



    // Conversion implicite si types différents,
    // mais seulement si le type cible n'est pas complexe
    // (sinon on pert implicitement la partie imaginaire !)
    template<typename T2>
    TabT(const TabT<T2,ndims> &t)
    NOECLIPSE(requires(!std::is_same<T,T2>()
             && ((is_complex<T>()) || (!is_complex<T2>()))))
    {
      fr = TabT<T,ndims>(t);
    }

    /*template<typename T2>
    TabT(TabT<T2,ndims> &&t)
    NOECLIPSE(requires(!std::is_same<T,T2>()
             && ((is_complex<T>()) || (!is_complex<T2>()))))
    {
      fr = TabT<T,ndims>(t);
    }*/



    /*TabT &operator =(const TabT &&t)
    {
      fr = t.fr;
      return *this;
    }*/

    explicit TabT(entier n, entier m = 1): fr(n, m)
    {
    }

    TabT<T,ndims> transpose() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.transpose();
    }

    void resize(entier n, entier m) NOECLIPSE(requires(ndims == 2))
    {
      fr.resize(n,m);
    }

    void resize(entier n) NOECLIPSE(requires(ndims == 1))
    {
      fr.resize(n);
    }

    void conservativeResize(entier n, entier m = 1)
    {
      fr.conservativeResize(n, m);
    }

    void set_ones(entier n = -1, entier m = 1)
    {
      fr.set_ones(n, m);
    }

    void setZero(entier n = -1, entier m = 1)
    {
      fr.setZero(n, m);
    }

    void setRandom()
    {
      fr.setRandom();
    }


    TabT<T, ndims> clone() const
    {
      return TabT<T,ndims>(fr.clone());
    }

    auto squared_norm() const
    {
      retourne fr.norme_carré();
    }

    auto norm() const
    {
      retourne fr.norme();
    }

    auto ms() const NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr.ms();
    }

    auto ms() const NOECLIPSE(requires(is_complex<T>()))
    {
      return fr.ms();
    }

    auto rms() const
    {
      return fr.rms();
    }


    T sum() const NOECLIPSE(requires(std::same_as<T,float>))
    {
      return fr.somme();
    }
  # if 1
    T sum() const NOECLIPSE(requires(std::same_as<T,std::complex<float>>))
    {
      return fr.somme();
    }

    T sum() const NOECLIPSE(requires(!tsd::base_float<T>))
    {
      return fr.somme();
    }


    TabT<T, 1> cwise_mean() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.moyenne_colonnes();
    }

    TabT<T, 1> rwise_mean() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.moyenne_lignes();
    }

    TabT<T, 1> cwise_sum() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.somme_colonnes();
    }

    TabT<T, 1> rwise_sum() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.somme_lignes();
    }

    TabT<T, ndims> add_row(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      return fr.ajoute_lignes(x);
    }

    TabT<T, ndims> mut_row(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      return fr.produit_lignes(x);
    }

    TabT<T, ndims> add_col(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      return fr.ajoute_colonnes(x);
    }

    TabT<T, ndims> mult_col(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      return fr.produit_colonnes(x);
    }


    T mean() const
    {
      return fr.moyenne();
    }

    TabT<T, 2> reverse_rows() const NOECLIPSE(requires(ndims == 2))
    {
      return fr.reverse_rows();
    }

    TabT<T, ndims> reverse() const NOECLIPSE(requires(ndims == 1))
    {
      return fr.reverse();
    }

    void reverseInPlace() NOECLIPSE(requires(ndims == 1))
    {
      fr.reverseInPlace();
    }


    TabT<T,2> reshape(entier m, entier n) const
    {
      return fr.reshape(m, n);
    }

    TabT<T,1> reshape(entier m) const
    {
      return fr.reshape(m);
    }

    entier dim() const NOECLIPSE(requires (ndims == 1))
    {
      return fr.dim();
    }

    // Extraction
    TabT<T,1> head(entier n) const NOECLIPSE(requires (ndims == 1))
    {
      return fr.head(n);
    }
    TabT<T,1> tail(entier n) const NOECLIPSE(requires (ndims == 1))
    {
      return fr.tail(n);
    }
    TabT<T,1> segment(entier i0, entier n) const NOECLIPSE(requires (ndims == 1))
    {
      return fr.segment(i0, n);
    }

    TabT<T,2> block(entier r0, entier nr, entier c0, entier nc) const NOECLIPSE(requires (ndims == 2))
    {
      return fr.block(r0, nr, c0, nc);
    }

    TabT<T,1> col(entier num) const NOECLIPSE(requires (ndims == 2))
    {
      return fr.col(num);
      /*msg("Tab/en::col -> fr::col...");
      soit fr_tab = fr.col(num);
      //msg("Tab/en::col -> fr infos:");
      fr_tab.dump_infos();
      msg("Tab/en::col -> fr -> en...");
      return fr_tab;*/
    }

    TabT<T,1> row(entier num) const NOECLIPSE(requires (ndims == 2))
    {
      return fr.row(num);
    }

    void set_row(entier num, const TabT<T, 1> &x) NOECLIPSE(requires (ndims == 2))
    {
      fr.set_row(num, x);
    }

    T dot(const TabT<T,ndims> &x) const NOECLIPSE(requires (ndims == 1))
    {
      return fr.dot(x);
    }

    void set_real(const auto &x) NOECLIPSE(requires(is_complex<T>()))
    {
      fr.set_real(x.fr);
    }

    void set_imag(const auto &x) const NOECLIPSE(requires(is_complex<T>()))
    {
      fr.set_imag(x.fr);
    }

    TabT<T,ndims> conjugate() const
    {
      return fr.conjugate();
    }

    TabT<std::complex<T>,ndims> as_complex() const NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr.as_complex();
    }

    TabT<T,ndims> as_complex() const NOECLIPSE(requires(is_complex<T>()))
    {
      return fr.as_complex();
    }

    static TabT empty()
    {
      return tsd::TabT<T,ndims>::vide();
    }

    bool is_empty() const
    {
      return fr.est_nul();
    }


    static TabT zeros(entier i, entier j = 1)
    {
      return tsd::TabT<T, ndims>::zeros(i, j);
    }

    void setConstant(T c)
    {
      fr.setConstant(c);
    }

    TabT &operator =(T c)
    {
      fr = c;
      return *this;
    }

    static TabT constant(entier n, entier m, T valeur) NOECLIPSE(requires(ndims == 2))
    {
      return tsd::TabT<T, ndims>::constant(n, m, valeur);
    }

    static TabT constant(entier n, T valeur) NOECLIPSE(requires(ndims == 1))
    {
      return tsd::TabT<T, ndims>::constant(n, valeur);
    }

    static TabT values(std::initializer_list<T> vals) NOECLIPSE(requires(ndims == 1))
    {
      return tsd::TabT<T, ndims>::valeurs(vals);
    }

    static TabT values(int nrows, int ncols, std::initializer_list<T> vals) NOECLIPSE(requires(ndims == 2))
    {
      return tsd::TabT<T, ndims>::valeurs(nrows, ncols, vals);
    }

    static TabT map(vector<T> &x)
    {
      return tsd::TabT<T, ndims>::map(x);
    }

    static TabT map(T *ptr, int n, int m = 1)
    {
      return tsd::TabT<T, ndims>::map(ptr, n, m);
    }

    static const TabT map(const T *ptr, int n, int m = 1)
    {
      return tsd::TabT<T, ndims>::map(ptr, n, m);
    }

    static TabT ones(entier i, entier j = 1)
    {
      return tsd::TabT<T, ndims>::ones(i, j);
    }

    static TabT eye(entier i)
    {
      return tsd::TabT<T, ndims>::eye(i);
    }

    static TabT random(entier n) NOECLIPSE(requires(ndims == 1))
    {
      return tsd::TabT<T, ndims>::random(n);
    }

    static TabT random(entier n, entier m) NOECLIPSE(requires(ndims == 2))
    {
      return tsd::TabT<T, ndims>::random(n, m);
    }

    TabT unaryExpr(auto f) const
    {
      return fr.unaryExpr(f);
    }

    static TabT int_expr(entier n, auto f) NOECLIPSE(requires(ndims == 1))
    {
      return tsd::TabT<T,ndims>::int_expr(n, f);
    }

    static TabT int_expr(entier n, entier m, auto f) NOECLIPSE(requires(ndims == 2))
    {
      return tsd::TabT<T,ndims>::int_expr(n, m, f);
    }

//  #   define IMAP(BB) [&](int i) {return BB;}

    inline T &operator()(entier i, entier j)
    {
      return fr(i, j);
    }

    inline T &operator()(entier i)
    {
      return fr(i);
    }

    inline const T &operator()(entier i, entier j) const
    {
      return fr(i, j);
    }

    inline const T &operator()(entier i) const
    {
      return fr(i);
    }

    inline T *data()
    {
      return fr.data();
    }

    inline const T *data() const
    {
      return fr.data();
    }


    TabT<T,ndims> operator-( ) const
    {
      return -fr;
    }

    TabT<T,ndims> inv() const
    {
      return fr.inv();
    }

    TabT<T, 1> matprod(const TabT<T, 1> &x) const
    {
      return fr.matprod(x.fr);
    }

    TabT<T, 2> matprod(const TabT<T, 2> &x) const
    {
      return fr.matprod(x.fr);
    }

    TabT<T,ndims> &operator *=(const T &x)
    {
      fr *= x;
      return *this;
    }
    TabT<T,ndims> &operator *=(const TabT<T,ndims> &x)
    {
      fr *= x.fr;
      return *this;
    }
    TabT<T,ndims> &operator /=(const T &x)
    {
      fr /= x;
      return *this;
    }
    TabT<T,ndims> &operator /=(const TabT<T,ndims> &x)
    {
      fr /= x.fr;
      return *this;
    }
    TabT<T,ndims> &operator +=(const T &x)
    {
      fr += x;
      return *this;
    }
    TabT<T,ndims> &operator +=(const TabT<T,ndims> &x)
    {
      fr += x.fr;
      return *this;
    }
    TabT<T,ndims> &operator -=(const T &x)
    {
      fr -= x;
      return *this;
    }
    TabT<T,ndims> &operator -=(const TabT<T,ndims> &x)
    {
      fr -= x.fr;
      return *this;
    }


    TabT<T,ndims> operator &&(const TabT<T,ndims> &x) const NOECLIPSE(requires(std::same_as<T,char>))
    {
      return fr && x.fr;
    }

    TabT<T,ndims> operator ||(const TabT<T,ndims> &x) const NOECLIPSE(requires(std::same_as<T,char>))
    {
      return fr || x.fr;
    }


    TabT<T,ndims> operator *(const T &x) const
    {
      return fr * x;
    }

    TabT<T,ndims> operator +(const T &x) const
    {
      return fr + x;
    }

    TabT<T,ndims> operator /(const T &x) const
    {
      return fr / x;
    }

    TabT<T,ndims> operator -(const T &x) const
    {
      return fr - x;
    }

    TabT<std::complex<T>,ndims> operator +(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr + x.fr;
    }

    TabT<T,ndims> operator +(const TabT<T,ndims> &x) const
    {
      return fr + x.fr;
    }

    TabT<std::complex<T>,ndims> operator *(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr * x.fr;
    }

    TabT<T,ndims> operator *(const TabT<T,ndims> &x) const
    {
      return fr * x.fr;
    }

    TabT<std::complex<T>,ndims> operator -(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr - x.fr;
    }

    TabT<T,ndims> operator -(const TabT<T,ndims> &x) const
    {
      return fr - x.fr;
    }

    TabT<std::complex<T>,ndims> operator /(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!is_complex<T>()))
    {
      return fr / x.fr;
    }

    TabT<T,ndims> operator /(const TabT<T,ndims> &x) const
    {
      return fr / x.fr;
    }

    TabT<bool_t,ndims> operator >(const TabT<T,ndims> &t) const
    {
      return fr > t.fr;
    }
    TabT<bool_t,ndims> operator <(const TabT<T,ndims> &t) const
    {
      return fr < t.fr;
    }
    TabT<bool_t,ndims> operator ==(const TabT<T,ndims> &t) const
    {
      return fr == t.fr;
    }
    TabT<bool_t,ndims> operator >=(const TabT<T,ndims> &t) const
    {
      return fr >= t.fr;
    }
    TabT<bool_t,ndims> operator <=(const TabT<T,ndims> &t) const
    {
      return fr <= t.fr;
    }

    TabT<bool_t,ndims> operator >(double x) const
    {
      return fr > x;
    }
    TabT<bool_t,ndims> operator <(double x) const
    {
      return fr < x;
    }
    TabT<bool_t,ndims> operator ==(double x) const
    {
      return fr == x;
    }
    TabT<bool_t,ndims> operator >=(double x) const
    {
      return fr >= x;
    }
    TabT<bool_t,ndims> operator <=(double x) const
    {
      return fr <= x;
    }



  #   endif
  };


template<typename T, entier ndims, typename T2>
  TabT<T,ndims> operator *(const T2 &x, const TabT<T,ndims> &x1) NOECLIPSE(requires(std::is_scalar_v<T2>))
{
  retourne x1 * ((T) x);
}

template<typename T, entier ndims, typename T2>
  TabT<T,ndims> operator +(const T2 &x, const TabT<T,ndims> &x1) NOECLIPSE(requires(std::is_scalar_v<T2>))
{
  retourne x1 + ((T) x);
}


template<typename T, entier ndims, typename T2>
  TabT<T,ndims> operator /(const T2 &x, const TabT<T,ndims> &x1) NOECLIPSE(requires(std::is_scalar_v<T2>))
{
  retourne ((T) x) * x1.inv();
}


template<typename T, entier ndims, typename T2>
  TabT<T,ndims> operator -(const T2 &x, const TabT<T,ndims> &x1) NOECLIPSE(requires(std::is_scalar_v<T2>))
{
  retourne (-x1) + ((T) x);
}



template<typename T>
using VecT = TabT<T, 1>;
using Vecf = VecT<float>;
using Vecd = VecT<double>;
using Veccf = VecT<cfloat>;
using Veccd = VecT<cdouble>;
using Vecb  = VecT<char>;
using Veci  = VecT<int32_t>;


using Tabi  = TabT<int32_t, 2>;
using Tabf  = TabT<float, 2>;
using Tabd  = TabT<double, 2>;
using Tabcf  = TabT<cfloat, 2>;
using Tabcd  = TabT<cdouble, 2>;

template<typename T, int ndims>
inline std::ostream& operator<<(std::ostream& strm, const TabT<T,ndims> &t){strm << t.fr; return strm;}




#define EN_DOP1(AA,BB)\
  template<typename T, int ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x)\
  {return tsd::BB(x.fr);}


#define EN_DOP1C2R(AA,BB)\
  template<typename T, int ndims>\
  TabT<T, ndims> AA(const TabT<std::complex<T>, ndims> &x)\
  {return tsd::BB(x.fr);}\
  template<typename T, int ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x)\
  NOECLIPSE(requires (!is_complex<T>()))\
  {return tsd::BB(x.fr);}


#define EN_DOP2(AA,BB)\
  template<typename T, int ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x, const TabT<T, ndims> &y)\
  {return tsd::BB(x.fr,y.fr);}

EN_DOP1(sqrt,sqrt)
EN_DOP1(square,square)
EN_DOP1(sin,sin)
EN_DOP1(cos,cos)
EN_DOP1(log,log)
EN_DOP1(log10,log10)
EN_DOP1(exp,exp)
EN_DOP1(sign,signe)
// complex -> float
EN_DOP1C2R(abs,abs)
EN_DOP1C2R(abs2,abs2)
EN_DOP1C2R(arg,arg)

EN_DOP2(pow,pow)


template <typename T, int ndims>
TabT<T,ndims> real(const TabT<std::complex<T>, ndims> &x)
{
  return real(x.fr);
}


template <typename T, int ndims>
TabT<T,ndims> real(const TabT<T, ndims> &x)NOECLIPSE( requires(!is_complex<T>()))
{
  return x;
}

template <typename T, entier ndims>
TabT<T,ndims> imag(const TabT<std::complex<T>, ndims> &x)
{
  return imag(x.fr);
}


}


//ostream_formater(dsp::TabT<>)

template <typename T, int ndims> struct fmt::formatter<dsp::TabT<T,ndims>> { \
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();} \
  template <typename FormatContext> \
  auto format(const dsp::TabT<T,ndims>& t, FormatContext& ctx) const -> decltype(ctx.out())  \
  { \
    std::ostringstream ss; \
    ss << t; \
    return fmt::format_to(ctx.out(), "{}", ss.str()); \
  } \
};




#endif
