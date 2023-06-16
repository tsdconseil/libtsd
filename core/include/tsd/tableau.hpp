#ifndef TAB_H
#define TAB_H

#include "tsd/commun.hpp"

#include <complex>
#include <cmath>
#include <ostream>
#include <sstream>

#define TSD_SAFE 1



namespace dsp
{
template<typename T, entier ndims>
  struct TabT;
}

namespace tsd
{

using bool_t = char;

/** @brief Raccourci pour nombre complexe flottant 32 bits */
typedef std::complex<float> cfloat;

/** @brief Raccourci pour nombre complexe flottant 64 bits */
typedef std::complex<double> cdouble;


template<typename T>
struct est_complex_t : public std::false_type {};

template<typename T>
struct est_complex_t<std::complex<T>> : public std::true_type {};

template<typename T>
constexpr bouléen est_complexe() { retourne est_complex_t<T>::value; }


#ifdef __CDT_PARSER__
#define NOECLIPSE(AA)
#else
#define NOECLIPSE(AA) AA
#endif


template<typename T, entier ndims>
  struct TabT;

typedef enum scalaire_enum
{
  ℝ, ℂ, ℕ, ℤ, B//, 𝔹
} Scalaire;

struct NTenseurDim
{
  vector<entier> dims;

  static NTenseurDim dim1(entier n);
  static NTenseurDim dim2(entier n, entier m);

  inline const entier &operator[](entier i) const{retourne dims[i];}
  inline entier &operator[](entier i){retourne dims[i];}
  inline entier operator()(entier i) const{si(i >= ndims()) retourne 1; retourne dims[i];}
  inline entier ndims() const{retourne dims.size();}

  entier total() const;

  bouléen operator ==(const NTenseurDim &d2) const;

  NTenseurDim operator *(const NTenseurDim &d2) const;
};

  struct Tab
  {

    template<typename T, entier ndims>
    Tab(const TabT<T, ndims> &t);
    Tab();
    Tab(Scalaire s, entier reso, const NTenseurDim &dim);
    Tab(Scalaire s, entier reso, entier n, entier m = 1);

    const NTenseurDim &get_dims() const;

    static Tab ones (Scalaire s, entier reso, const NTenseurDim &dim);
    static Tab zeros(Scalaire s, entier reso, const NTenseurDim &dim);
    static Tab diag (const Tab &src, entier ndims);

    static Tab ones(Scalaire s, entier reso, entier n, entier m = 1);
    static Tab zeros(Scalaire s, entier reso, entier n, entier m = 1);
    static Tab eye(Scalaire s, entier reso, entier n);
    static Tab constant(Scalaire s, entier reso, entier n, entier m, double valeur);
    static Tab random(Scalaire s, entier reso, entier n, entier m);
    static Tab map(Scalaire s, entier reso, entier n, entier m, void *data);

    Tab même_format() const;

    // Extraction
    Tab head(entier n) const;
    Tab tail(entier n) const;
    Tab segment(entier i0, entier n) const;


    bouléen est_reference() const;

    Tab block(entier r0, entier nr, entier c0, entier nc) const;

    // TODO: restrict type not complex ?
    double maxCoeff(entier *i = nullptr, entier *j = nullptr) const;
    double minCoeff(entier *i = nullptr, entier *j = nullptr) const;

    Scalaire tscalaire() const;
    entier dim_scalaire() const;
    entier nbits() const;

    Tab clone() const;

    Tab cwiseMax(double x) const;

    bouléen est_nul() const;


    // TODO: 1e-12 pour double
    bouléen est_approx(const Tab &x, double tol = 1e-5) const;
    bouléen est_approx(double x, double tol = 1e-5) const;

    void resize(entier n, entier m = 1);
    void conservativeResize(entier n, entier m = 1);

    void setZero(entier n = -1, entier m = 1);
    void set_ones(entier n = -1, entier m = 1);

    Tab transpose_int() const;
    Tab conjugate() const;

    inline entier get_ndims()      const{retourne get_dims().ndims();}
    entier rows()       const;
    entier cols()       const;
    entier nelems()     const;

    // ASUPPRIMER
    void evaluer() const;
    // ASUPPRIMER
    Tab eval()   const;

    bouléen hasNaN() const;

    Tab col(entier num) const;
    Tab row(entier num) const;

    // Réductions
    double somme()   const;
    double produit() const;
    double moyenne() const;

    bool est_de_même_dimensions(const Tab &t) const;

    void setConstant(double c);
    //void setOnes(entier n, entier m = 1);

    Tab reshape(entier m, entier n = 1) const;

    template<typename T, entier ndims = 1>
      TabT<T, ndims> as() const;


    template<typename T, entier ndims = 1>
      TabT<T, ndims> as_no_conv() const;

    Tab reverse() const;
    void reverseInPlace();

    Tab reverse_rows() const;

    Tab &operator =(const Tab &t);

    void copie(const Tab &source);


    Tab operator +(const Tab &t) const;
    Tab operator -(const Tab &t) const{retourne operator_minus(t);}
    Tab operator_minus(const Tab &t) const;
    Tab operator *(const Tab &t) const;
    Tab operator /(const Tab &t) const;

    Tab matprod(const Tab &t) const;

    Tab operator-( ) const;
    Tab inv() const;
    Tab lsq(const Tab &x) const;

    Tab inverse_matricielle() const;

#   define DEC_COMP(AA) \
    Tab operator AA(const Tab &t) const;\
    Tab operator AA(double x) const;

    DEC_COMP(>);
    DEC_COMP(<);
    DEC_COMP(==);
    DEC_COMP(>=);
    DEC_COMP(<=);

#   define DEC_OP(TYPE_SCAL)\
    Tab operator *(const TYPE_SCAL &x)   const;\
    Tab &operator *=(const TYPE_SCAL &x);\
    Tab operator /(const TYPE_SCAL &x)   const;\
    Tab &operator /=(const TYPE_SCAL &x);\
    Tab operator +(const TYPE_SCAL &x)   const;\
    Tab &operator +=(const TYPE_SCAL &x);\
    Tab operator -(const TYPE_SCAL &x)   const;\
    Tab &operator -=(const TYPE_SCAL &x);

    DEC_OP(entier);
    DEC_OP(float);
    DEC_OP(double);
    DEC_OP(cfloat);
    DEC_OP(cdouble);


    Tab &operator *=(const Tab &x);
    Tab &operator /=(const Tab &x);
    Tab &operator +=(const Tab &x);
    Tab &operator -=(const Tab &x);

    void *rawptr();
    const void *rawptr() const;


    void dump_infos() const;

    Tab select(const Tab &x1, const Tab &x2) const;

  //private:
    struct Impl;
    sptr<Impl> impl;
  };




  /** Retourne le type et la dimension sous forme dynamique à partir du type */
  template<typename T>
    tuple<Scalaire,entier> T2S()
  {
    if constexpr(std::is_same<T, float>())
      retourne {Scalaire::ℝ, 32};
    else if constexpr(std::is_same<T, double>())
      retourne {Scalaire::ℝ, 64};
    else if constexpr(std::is_same<T, cfloat>())
      retourne {Scalaire::ℂ, 32*2};
    else if constexpr(std::is_same<T, cdouble>())
      retourne {Scalaire::ℂ, 64*2};
    else if constexpr(std::is_same<T, int32_t>())
      retourne {Scalaire::ℤ, 32};
    else if constexpr(std::is_same<T, char>())
      retourne {Scalaire::B, 8};
    else
      retourne {Scalaire::ℝ, 32};
  }

  template<typename T>
    concept base_float = (std::same_as<T,float> || std::same_as<T,std::complex<float>>);



  using NTenseur = Tab;

// Tableau de données évaluées, avec type spécifique, et nombre de dimensions spécifique
template<typename T, entier ndims>
struct TabT: Tab
{
  TabT(const dsp::TabT<T,ndims> &x): TabT(x.fr)
  {

  }

  T produit() const
  {
    entier n = nelems();
    T x = 1;
    Pour(auto i = 0; i < n; i++)
      x *= raw_data[i];
    retourne x;
  }

  TabT<T,ndims> cwiseSign() const NOECLIPSE(requires(!est_complexe<T>()))
  {
    TabT<T,ndims> y(rows(), cols());
    soit n = nelems();

    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = raw_data[i] >= 0 ? 1 : -1;

    retourne y;
  }


  TabT<T,ndims> cwiseMax(T x) const NOECLIPSE(requires(!est_complexe<T>()))
  {
    TabT<T,ndims> y(rows(), cols());
    soit n = nelems();

    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = std::max(x, raw_data[i]);

    retourne y;
  }

  TabT<T,ndims> cwiseMax(const TabT<T,ndims> &x) const NOECLIPSE(requires(!est_complexe<T>()))
  {
    TabT<T,ndims> y(rows(), cols());
    soit n = nelems();

    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = std::max(x.raw_data[i], raw_data[i]);

    retourne y;
  }

  TabT<T,ndims> cwiseMin(const TabT<T,ndims> &x) const NOECLIPSE(requires(!est_complexe<T>()))
  {
    TabT<T,ndims> y(rows(), cols());
    soit n = nelems();

    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = std::min(x.raw_data[i], raw_data[i]);

    retourne y;
  }

  bool tous_vrai() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    retourne valeur_min() > 0;
  }

  bool tous_faux() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    retourne valeur_max() == 0;
  }

  bool au_moins_un_vrai() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    retourne valeur_max() > 0;
  }

  bool au_moins_un_faux() const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    retourne valeur_min() == 0;
  }

  template<typename T2>
  TabT<T2,ndims> select(const TabT<T2,ndims> &x1, const TabT<T2,ndims> &x2) const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    TabT<T2,ndims> y(rows(), cols());
    soit n = nelems();

    // TODO: itérateur qui marche sur les blocs
    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = raw_data[i] ? x1.raw_data[i] : x2.raw_data[i];

    retourne y;
  }

  template<typename T2>
  TabT<T2,ndims> select(const TabT<T2,ndims> &x1, T2 x2) const NOECLIPSE(requires(std::is_same<T, char>::value))
  {
    TabT<T2,ndims> y(rows(), cols());
    soit n = nelems();

    // TODO: itérateur qui marche sur les blocs
    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = raw_data[i] ? x1.raw_data[i] : x2;

    retourne y;
  }

  TabT<T,ndims> cwiseMin(T x) const NOECLIPSE(requires(!est_complexe<T>()))
  {
    TabT<T,ndims> y(rows(), cols());
    soit n = nelems();

    pour(auto i = 0; i < n; i++)
      y.raw_data[i] = std::min(x, raw_data[i]);

    retourne y;
  }

  entier index_max() const NOECLIPSE(requires(ndims == 1))
  {
    entier i;
    maxCoeff(&i);
    retourne i;
  }

  entier index_min() const NOECLIPSE(requires(ndims == 1))
  {
    entier i;
    minCoeff(&i);
    retourne i;
  }

  T valeur_max() const
  {
    retourne maxCoeff();
  }

  T valeur_min() const
  {
    retourne minCoeff();
  }

  tuple<T, T> valeurs_minmax() const NOECLIPSE(requires(!est_complexe<T>()))
  {
    retourne {minCoeff(), maxCoeff()};
  }

  tuple<entier, entier> index_max() const NOECLIPSE(requires((ndims == 2) && !est_complexe<T>()))
  {
    entier i, j;
    maxCoeff(&i, &j);
    retourne {i, j};
  }

  tuple<entier, entier> index_min() const NOECLIPSE(requires((ndims == 2) && !est_complexe<T>()))
  {
    entier i, j;
    minCoeff(&i, &j);
    retourne {i, j};
  }


  tuple<T, entier> max() const NOECLIPSE(requires((ndims == 1) && !est_complexe<T>()))
  {
    entier i;
    soit y = maxCoeff(&i);
    retourne {y, i};
  }

  tuple<T, entier, entier> max() const NOECLIPSE(requires((ndims == 2) && !est_complexe<T>()))
  {
    entier i, j;
    soit y = maxCoeff(&i, &j);
    retourne {y, i, j};
  }

  tuple<T, entier> min() const NOECLIPSE(requires((ndims == 1) && !est_complexe<T>()))
  {
    entier i;
    soit y = minCoeff(&i);
    retourne {y, i};
  }

  tuple<T, entier, entier> min() const NOECLIPSE(requires((ndims == 2) && !est_complexe<T>()))
  {
    entier i, j;
    soit y = minCoeff(&i, &j);
    retourne {y, i, j};
  }


    TabT()
    {
      soit [K, nb]= T2S<T>();
      Tab tmp(K, nb, 0, 0);
      impl      = tmp.impl;
      raw_data  = (T *) rawptr();
    }


    TabT(const Tab &t)
    {
      // Vérification dynamique de la compatibilité
      soit [K, nb]= T2S<T>();

      if((t.nbits() != nb) || (t.tscalaire() != K))
      {
        échec("Conversion Tab[{},{}] -> TabT[{},{}]: type non compatible", t.tscalaire(), t.nbits(), K, nb);
      }
      if((ndims == 1) && (t.cols() > 1))
      {
        échec("Conversion Tab[ndims=2] -> TabT[ndims=1]: nombres de dimensions non compatibles");
      }
      impl = t.impl;
      raw_data  = (T *) rawptr();
    }

    // Conversion implicite si types différents,
    // mais seulement si le type cible n'est pas complexe
    // (sinon on pert implicitement la partie imaginaire !)
    template<typename T2>
    TabT(const TabT<T2,ndims> &t)
    NOECLIPSE(requires(!std::is_same<T,T2>()
             && ((est_complexe<T>()) || (!est_complexe<T2>()))))
    {
      *this = Tab(t).as<T>();
    }

    /*TabT &operator =(const TabT &t)
    {
      Tab::operator =(t);
      this->raw_data = (T *) rawptr();
      retourne *this;
    }*/

    TabT(TabT &&t)
    {
      //msg("Tabt::TabT(&&) : move");
      this->impl      = t.impl;
      this->raw_data  = t.raw_data;

      t.impl      = sptr<Tab::Impl>();
      t.raw_data  = nullptr;
    }

    bool est_reference() const
    {
      return Tab::est_reference();
    }

    // From RHV
    TabT &operator =(TabT &&t)
    {
      //Tab::operator =(t);
      //this->raw_data = (T *) rawptr();

      // TODO here:
      // if self is a subvector and dimensions are compatible, do a move.
      si(est_reference() && est_de_même_dimensions(t))
      {
        //msg("TabT::operator =(&&) : copie.");
        copie(t);
      }
      sinon
      {
        //msg("TabT::operator =(&&) : move.");
        this->impl      = t.impl;
        this->raw_data  = t.raw_data;

        t.impl      = sptr<Impl>();
        t.raw_data  = nullptr;
      }

      retourne *this;
    }

    TabT(const TabT &t)
    {
      //msg("Tabt::TabT(const &) : clone");
      *this = t.clone();
    }

    // From other than RHV -> deep copy ?
    TabT &operator =(const TabT &t)
    {
      //msg("TabT::operator =(&) : est_reference = {}.", est_reference());
      si(est_reference() && est_de_même_dimensions(t))
      {
        //msg("TabT::operator =(&) : copie.");
        copie(t);
      }
      sinon
      {
        //msg("TabT::operator =(&) : clone.");
        *this = t.clone();
      }
      retourne *this;
    }

    explicit TabT(entier n, entier m = 1)
    {
      soit [K, nb]= T2S<T>();
      Tab res(K, nb, n, m);
      this->impl      = res.impl;
      this->raw_data  = (T *) res.rawptr();
    }

    TabT<T,ndims> transpose() const NOECLIPSE(requires(ndims == 2))
    {
      retourne Tab(*this).transpose_int();
    }

    TabT<T,ndims> inverse_matricielle() const NOECLIPSE(requires(ndims == 2))
    {
      retourne Tab(*this).inverse_matricielle();
    }

    void resize(entier n, entier m) NOECLIPSE(requires(ndims == 2))
    {
      Tab::resize(n, m);
      raw_data  = (T *) rawptr();
    }

    void resize(entier n) NOECLIPSE(requires(ndims == 1))
    {
      Tab::resize(n, 1);
      raw_data  = (T *) rawptr();
    }

    void conservativeResize(entier n, entier m = 1)
    {
      Tab::conservativeResize(n, m);
      raw_data  = (T *) rawptr();
    }

    void set_ones(entier n = -1, entier m = 1)
    {
      Tab::set_ones(n, m);
      raw_data  = (T *) rawptr();
    }

    void setZero(entier n = -1, entier m = 1)
    {
      Tab::setZero(n, m);
      raw_data  = (T *) rawptr();
    }

    void setRandom()
    {
      if constexpr(ndims == 2)
        *this = TabT<T, ndims>::random(rows(), cols());
      else
        *this = TabT<T, ndims>::random(rows());
    }


    TabT<T, ndims> clone() const
    {
      retourne Tab(*this).clone();
    }

    template<typename Ti>
    T somme_impl() const
    {
      soit n = nelems();
      Ti x = 0;
      Pour(auto i = 0; i < n; i++)
        x += raw_data[i];
      retourne (T) x;
    }


    auto ms() const NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne (square(*this)).moyenne();
    }

    auto ms() const NOECLIPSE(requires(est_complexe<T>()))
    {
      retourne (abs2(*this)).moyenne();
    }

    auto rms() const
    {
      retourne std::sqrt(ms());
    }

    auto norme() const
    {
      retourne std::sqrt(norme_carré());
    }

    auto norme_carré() const NOECLIPSE(requires(est_complexe<T>()))
    {
      retourne (abs2(*this)).somme();
    }

    auto norme_carré() const NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne (square(*this)).somme();
    }

    T somme() const NOECLIPSE(requires(std::same_as<T,float>))
    {
      retourne somme_impl<double>();
    }

    T somme() const NOECLIPSE(requires(std::same_as<T,std::complex<float>>))
    {
      retourne somme_impl<std::complex<double>>();
    }

    T somme() const NOECLIPSE(requires(!base_float<T>))
    {
      retourne somme_impl<T>();
    }


    TabT<T, 1> moyenne_colonnes() const NOECLIPSE(requires(ndims == 2))
    {
      soit n = cols();
      TabT<T,1> y(n);

      pour(auto i = 0; i < n; i++)
        y(i) = col(i).moyenne();

      retourne y;
    }

    TabT<T, 1> moyenne_lignes() const NOECLIPSE(requires(ndims == 2))
    {
      soit n = rows();
      TabT<T,1> y(n);

      pour(auto i = 0; i < n; i++)
        y(i) = row(i).moyenne();

      retourne y;
    }

    TabT<T, 1> somme_colonnes() const NOECLIPSE(requires(ndims == 2))
    {
      soit n = cols();
      TabT<T,1> y(n);

      pour(auto i = 0; i < n; i++)
        y(i) = col(i).somme();

      retourne y;
    }

    TabT<T, 1> somme_lignes() const NOECLIPSE(requires(ndims == 2))
    {
      soit n = rows();
      TabT<T,1> y(n);

      pour(auto i = 0; i < n; i++)
        y(i) = row(i).somme();

      retourne y;
    }

    TabT<T, ndims> ajoute_lignes(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      soit n = cols();
      assertion_msg(n == x.rows(),
          "ajoute_lignes(): le nombre de colonnes ({}) doit être égale à la dimension du vecteur ({}).",
          cols(), n);

      TabT<T, ndims> y(rows(), cols());

      pour(auto i = 0; i < n; i++)
        y.col(i) = col(i) + x(i);

      retourne y;
    }

    TabT<T, ndims> produit_lignes(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      soit n = cols();
      assertion_msg(n == x.rows(),
          "produit_lignes(): le nombre de colonnes ({}) doit être égale à la dimension du vecteur ({}).",
          n, x.rows());

      TabT<T, ndims> y(rows(), cols());

      pour(auto i = 0; i < n; i++)
        y.col(i) = col(i) * x(i);

      retourne y;
    }

    TabT<T, ndims> ajoute_colonnes(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      soit n = rows(), m = cols();
      assertion_msg(n == x.rows(),
          "ajoute_colonnes(): le nombre de lignes ({}) doit être égale à la dimension du vecteur ({}).",
          rows(), n);

      TabT<T, ndims> y(rows(), cols());

      pour(auto j = 0; j < m; j++)
        pour(auto i = 0; i < n; i++)
          y(i, j) = (*this)(i,j) + x(i);
        //y.col(i) = col(i) + x;

      retourne y;
    }

    TabT<T, ndims> produit_colonnes(const TabT<T, 1> &x) const NOECLIPSE(requires(ndims == 2))
    {
      soit n = rows(), m = cols();
      assertion_msg(n == x.rows(),
          "produit_colonnes(): le nombre de lignes ({}) doit être égale à la dimension du vecteur ({}).",
          rows(), n);

      TabT<T, ndims> y(rows(), cols());

      //pour(auto i = 0; i < n; i++)
        //y.col(i) = col(i) * x;

      pour(auto j = 0; j < m; j++)
        pour(auto i = 0; i < n; i++)
          y(i, j) = (*this)(i,j) * x(i);

      retourne y;
    }



    T moyenne() const
    {
      retourne somme() / ((T) nelems());
    }

    TabT<T, 2> reverse_rows() const NOECLIPSE(requires(ndims == 2))
    {
      retourne Tab(*this).reverse_rows();
    }

    TabT<T, ndims> reverse() const NOECLIPSE(requires(ndims == 1))
    {
      retourne Tab(*this).reverse();
    }

    void reverseInPlace() NOECLIPSE(requires(ndims == 1))
    {
      Tab(*this).reverseInPlace();
    }


    TabT<T,2> reshape(entier m, entier n) const
    {
      if(m * n != nelems())
      {
        échec("Tab::reshape : la dim totale doit être identique.");
        retourne {};
      }
      TabT<T,2> x(m, n);
      memcpy(x.rawptr(), this->rawptr(), m * n * dim_scalaire());
      retourne x;
    }

    TabT<T,1> reshape(entier m) const
    {
      if(m != nelems())
      {
        échec("Tab::reshape : la dim totale doit être identique.");
        retourne {};
      }
      TabT<T,1> x(m);
      memcpy(x.rawptr(), this->rawptr(), m * dim_scalaire());
      retourne x;
    }


    //template <entier nd2 = ndims>
    //typename std::enable_if<(nd2 == ndims) && (nd2 == 1)>::type dim()
    entier dim() const NOECLIPSE(requires (ndims == 1))
    {
      retourne rows();
    }


    // Extraction
    TabT<T,1> head(entier n) const NOECLIPSE(requires (ndims == 1))
    {
      retourne segment(0, n);
    }
    TabT<T,1> tail(entier n) const NOECLIPSE(requires (ndims == 1))
    {
      retourne segment(dim() - n, n);
    }
    TabT<T,1> segment(entier i0, entier n) const NOECLIPSE(requires (ndims == 1))
    {
      retourne Tab(*this).segment(i0, n);
    }

    TabT<T,2> block(entier r0, entier nr, entier c0, entier nc) const NOECLIPSE(requires (ndims == 2))
    {
      retourne Tab(*this).block(r0, nr, c0, nc);
    }

    TabT<T,1> col(entier num) const NOECLIPSE(requires (ndims == 2))
    {
      retourne Tab(*this).col(num);
    }

    TabT<T,1> row(entier num) const NOECLIPSE(requires (ndims == 2))
    {
      soit n = cols();
      TabT<T,1> y(n);

      pour(auto i = 0; i < n; i++)
        y(i) = (*this)(num, i);

      retourne y;
      //échec("TODO: row(...)");
      //retourne TabT<T,1>::zeros(cols());
    }

    void set_row(entier num, const TabT<T, 1> &x) NOECLIPSE(requires (ndims == 2))
    {
      soit n = cols();
      pour(auto i = 0; i < n; i++)
        (*this)(num, i) = x(i);
    }

    T dot(const TabT<T,ndims> &x) const NOECLIPSE(requires (ndims == 1))
    {
      T res = (T) 0;
      entier n = dim();
      assertion_msg(x.dim() == n, "Vec::dot(): les dimensions devraient être identiques ({} vs {})", n, x.dim());
      Pour(auto i = 0; i < n; i++)
        res += (*this)(i) * x(i);
      retourne res;
    }


    void set_real(const auto &x) NOECLIPSE(const requires(est_complexe<T>()))
    {
      entier n = nelems();
      Pour(auto i = 0; i < n; i++)
        raw_data[i].real(x.raw_data[i]);
    }


    void set_imag(const auto &x) const NOECLIPSE(requires(est_complexe<T>()))
    {
      entier n = nelems();
      Pour(auto i = 0; i < n; i++)
        raw_data[i].imag(x.raw_data[i]);
    }


    TabT<T,ndims> conjugate() const
    {
      retourne Tab(*this).conjugate();
    }

    TabT<std::complex<T>,ndims> as_complex() const NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne Tab(*this).as<std::complex<T>>();
    }

    TabT<T,ndims> as_complex() const NOECLIPSE(requires(est_complexe<T>()))
    {
      retourne *this;
    }

    bouléen est_vide() const
    {
      retourne nelems() == 0;
    }

    inline static TabT fromTab(const Tab &t)
    {
      retourne t.as_no_conv<T,ndims>();
    }

    static TabT zeros(entier i, entier j = 1)
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::zeros(K, nb, i, j));
    }

    void setConstant(T c)
    {
      if constexpr (ndims == 1)
      {
        pour(auto k = 0; k < dim(); k++)
          (*this)(k) = c;
      }
      else
      {
        pour(auto i = 0; i < cols(); i++)
          pour(auto j = 0; j < rows(); j++)
            (*this)(j,i) = c;
      }
    }

    TabT &operator =(T c)
    {
      setConstant(c);
      retourne *this;
    }

    static TabT diagonal(const TabT<T, 1> &v) NOECLIPSE(requires(ndims == 2))
    {
      soit n = v.rows();
      soit x = TabT<T,ndims>::zeros(n, n);
      pour(auto i = 0; i < n; i++)
        x(i,i) = v(i);
      retourne x;
    }

    static TabT constant(entier n, entier m, T valeur) NOECLIPSE(requires(ndims == 2))
    {
      TabT<T, ndims> x(n, m);
      x.setConstant(valeur);
      retourne x;
    }

    static TabT constant(entier n, T valeur) NOECLIPSE(requires(ndims == 1))
    {
      TabT<T, ndims> x(n);
      x.setConstant(valeur);
      retourne x;
    }

    static TabT valeurs(std::initializer_list<T> vals) NOECLIPSE(requires(ndims == 1))
    {
      soit [K, nb]= T2S<T>();
      soit res = TabT::fromTab(Tab(K, nb, vals.size(), 1));
      soit i = 0;
      pour(auto itr = vals.begin(); itr != vals.end(); itr++, i++)
        res(i) = *itr;
      retourne res;
    }

    static TabT valeurs(int nrows, int ncols, std::initializer_list<T> vals) NOECLIPSE(requires(ndims == 2))
    {
      soit [K, nb]= T2S<T>();

      assertion_msg(nrows * ncols == (int) vals.size(),
          "TabT::valeurs(rows={},cols={},...): le nombre d'élément doit être égal à rows*cols.", nrows, ncols);

      soit res = TabT::fromTab(Tab(K, nb, nrows, ncols));
      pour(auto i = 0; i < ncols; i++)
        pour(auto j = 0; j < nrows; j++)
          res(j,i) = *(vals.begin() + i+j*ncols);
      retourne res;
    }

    static TabT map(vector<T> &x)
    {
      soit [K, nb]= T2S<T>();
      soit res = TabT::fromTab(Tab(K, nb, x.size(), 1));
      memcpy(res.raw_data, x.data(), x.size() * sizeof(T));
      retourne res;
    }

    static TabT map(T *ptr, int n, int m = 1)
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::map(K, nb, n, m, ptr));
    }

    static const TabT map(const T *ptr, int n, int m = 1)
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::map(K, nb, n, m, (T *) ptr));
    }

    static TabT ones(entier i, entier j = 1)
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::ones(K, nb, i, j));
    }

    static TabT eye(entier i)
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::eye(K, nb, i));
    }

    static TabT random(entier n) NOECLIPSE(requires(ndims == 1))
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::random(K, nb, n, 1));
    }

    static TabT random(entier n, entier m) NOECLIPSE(requires(ndims == 2))
    {
      soit [K, nb]= T2S<T>();
      retourne TabT::fromTab(Tab::random(K, nb, n, m));
    }

    TabT unaryExpr(auto f) const
    {
      TabT y = même_format();
      soit ne = nelems();
      pour(auto i = 0; i < ne; i++)
        y(i) = f(raw_data[i]);
      retourne y;
    }

    static TabT int_expr(entier n, auto f) NOECLIPSE(requires(ndims == 1))
    {
      soit [K, nb]= T2S<T>();
      soit x = TabT::fromTab(Tab(K, nb, n, 1));
      pour(auto i = 0; i < n; i++)
        x.raw_data[i] = f(i);
      retourne x;
    }

    static TabT int_expr(entier n, entier m, auto f) NOECLIPSE(requires(ndims == 2))
    {
      soit [K, nb]= T2S<T>();
      soit x = TabT::fromTab(Tab(K, nb, n, 1));
      pour(auto j = 0; j < m; j++)
        pour(auto i = 0; i < n; i++)
          x.raw_data[i+j*n] = f(i,j);
      retourne x;
    }

#   define IMAP(BB) [&](entier i) {retourne BB;}


    inline T &operator()(entier i, entier j)
    {
      retourne raw_data[i + j * rows()];
    }

    inline T &operator()(entier i)
    {
#     if TSD_SAFE
      Si((i < 0) || (i >= rows()))
      {
        échec("Tab::operator()(i = {}) : hors borne (dim = {})", i, rows());
        //retourne raw_data[0];
      }
#     endif
      retourne raw_data[i];
    }

    inline const T &operator()(entier i, entier j) const
    {
#     if TSD_SAFE
      Si((i < 0) || (i >= rows()) || (j < 0) || (j >= cols()))
      {
        échec("Tab::operator()(i={}, j={}) : hors borne (dim = {}x{})", i, j, rows(), cols());
        //retourne raw_data[0];
      }
#     endif
      retourne raw_data[i + j * rows()];
    }

    inline const T &operator()(entier i) const
    {
#     if TSD_SAFE
      Si((i < 0) || (i >= rows()))
      {
        échec("Tab::operator()(i = {}) : hors borne (dim = {})", i, rows());
        //retourne raw_data[0];
      }
#     endif
      retourne raw_data[i];
    }

    inline T *data()
    {
      retourne raw_data;
    }

    inline const T *data() const
    {
      retourne raw_data;
    }


    TabT<T,ndims> operator-( ) const
    {
      retourne -Tab(*this);
    }

    TabT<T,ndims> inv() const
    {
      retourne Tab(*this).inv();
    }

    TabT<T,1> lsq(const TabT<T, 1> &x) const
    {
      retourne Tab(*this).lsq(x);
    }

    TabT<T, 1> matprod(const TabT<T, 1> &x) const
    {
      retourne Tab(*this).matprod(x);
    }

    TabT<T, 2> matprod(const TabT<T, 2> &x) const
    {
      retourne Tab(*this).matprod(x);
    }

    TabT<T,ndims> &operator *=(const T &x)
    {
      Tab(*this) *= x;
      retourne *this;
    }
    TabT<T,ndims> &operator *=(const TabT<T,ndims> &x)
    {
      Tab(*this) *= x;
      retourne *this;
    }
    TabT<T,ndims> &operator /=(const T &x)
    {
      Tab(*this) /= x;
      retourne *this;
    }
    TabT<T,ndims> &operator /=(const TabT<T,ndims> &x)
    {
      Tab(*this) /= x;
      retourne *this;
    }
    TabT<T,ndims> &operator +=(const T &x)
    {
      Tab(*this) += x;
      retourne *this;
    }
    TabT<T,ndims> &operator +=(const TabT<T,ndims> &x)
    {
      Tab(*this) += x;
      retourne *this;
    }
    TabT<T,ndims> &operator -=(const T &x)
    {
      Tab(*this) -= x;
      retourne *this;
    }
    TabT<T,ndims> &operator -=(const TabT<T,ndims> &x)
    {
      Tab(*this) -= x;
      retourne *this;
    }


    TabT<T,ndims> operator &&(const TabT<T,ndims> &x) const NOECLIPSE(requires(std::same_as<T,char>))
    {
      soit [K, nb]= T2S<T>();
      soit y = TabT::fromTab(Tab(K, nb, x.rows(), x.cols()));
      soit n = nelems();

      pour(auto i = 0; i < n; i++)
        y(i) = x(i) && (*this)(i);

      retourne y;
    }

    TabT<T,ndims> operator ||(const TabT<T,ndims> &x) const NOECLIPSE(requires(std::same_as<T,char>))
    {
      soit [K, nb]= T2S<T>();
      soit y = TabT::fromTab(Tab(K, nb, x.rows(), x.cols()));
      soit n = nelems();

      pour(auto i = 0; i < n; i++)
        y(i) = x(i) || (*this)(i);

      retourne y;
    }


    TabT<T,ndims> operator *(const T &x) const
    {
      retourne Tab(*this) * x;
    }

    TabT<T,ndims> operator +(const T &x) const
    {
      retourne Tab(*this) + x;
    }

    TabT<T,ndims> operator /(const T &x) const
    {
      retourne Tab(*this) / x;
    }

    TabT<T,ndims> operator -(const T &x) const
    {
      retourne Tab(*this) - x;
    }

    TabT<std::complex<T>,ndims> operator +(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne this->as<std::complex<T>>() + x;
    }

    TabT<T,ndims> operator +(const TabT<T,ndims> &x) const
    {
      retourne Tab(*this) + x;
    }

    TabT<std::complex<T>,ndims> operator *(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne this->as<std::complex<T>>() * x;
    }

    TabT<T,ndims> operator *(const TabT<T,ndims> &x) const
    {
      retourne Tab(*this) * x;
    }

    TabT<std::complex<T>,ndims> operator -(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne this->as<std::complex<T>>() - x;
    }

    TabT<T,ndims> operator -(const TabT<T,ndims> &x) const
    {
      retourne Tab(*this).operator_minus(x);
    }

    TabT<std::complex<T>,ndims> operator /(const TabT<std::complex<T>,ndims> &x) const
    NOECLIPSE(requires(!est_complexe<T>()))
    {
      retourne this->as<std::complex<T>>() / x;
    }

    TabT<T,ndims> operator /(const TabT<T,ndims> &x) const
    {
      retourne Tab(*this) / x;
    }

    TabT<bool_t,ndims> operator >(const TabT<T,ndims> &t) const
    {
      retourne Tab(*this) > Tab(t);
    }
    TabT<bool_t,ndims> operator <(const TabT<T,ndims> &t) const
    {
      retourne Tab(*this) < Tab(t);
    }
    TabT<bool_t,ndims> operator ==(const TabT<T,ndims> &t) const
    {
      retourne Tab(*this) == Tab(t);
    }
    TabT<bool_t,ndims> operator >=(const TabT<T,ndims> &t) const
    {
      retourne Tab(*this) >= Tab(t);
    }
    TabT<bool_t,ndims> operator <=(const TabT<T,ndims> &t) const
    {
      retourne Tab(*this) <= Tab(t);
    }

    TabT<bool_t,ndims> operator >(double x) const
    {
      retourne Tab(*this) > x;
    }
    TabT<bool_t,ndims> operator <(double x) const
    {
      retourne Tab(*this) < x;
    }
    TabT<bool_t,ndims> operator ==(double x) const
    {
      retourne Tab(*this) == x;
    }
    TabT<bool_t,ndims> operator >=(double x) const
    {
      retourne Tab(*this) >= x;
    }
    TabT<bool_t,ndims> operator <=(double x) const
    {
      retourne Tab(*this) <= x;
    }


    T *raw_data = nullptr;
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









template <typename T, entier ndims>
TabT<T,ndims> real(const TabT<std::complex<T>, ndims> &x)
{
  TabT<T,ndims> y(x.rows(), x.cols());
  entier n = x.nelems();
  pour(auto i = 0; i < n; i++)
    y.raw_data[i] = real(x.raw_data[i]);
  retourne y;
}


template <typename T, entier ndims>
TabT<T,ndims> real(const TabT<T, ndims> &x)NOECLIPSE( requires(!est_complexe<T>()))
{
  retourne x;
}

template <typename T, entier ndims>
TabT<T,ndims> imag(const TabT<std::complex<T>, ndims> &x)
{
  TabT<T,ndims> y(x.rows(), x.cols());
  entier n = x.nelems();
  pour(auto i = 0; i < n; i++)
    y.raw_data[i] = imag(x.raw_data[i]);
  retourne y;
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


#ifndef ostream_formater
#define ostream_formater(T) \
template <> struct fmt::formatter<T> { \
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();} \
  template <typename FormatContext> \
  auto format(const T& t, FormatContext& ctx) const -> decltype(ctx.out())  \
  { \
    std::ostringstream ss; \
    ss << t; \
    return fmt::format_to(ctx.out(), "{}", ss.str()); \
  } \
};
#endif


extern std::ostream& operator<<(std::ostream& strm, const NTenseurDim &t);
extern std::ostream& operator<<(std::ostream& strm, const Tab &t);
extern std::ostream& operator<<(std::ostream& strm, const Scalaire &s);



/////////////////////////////////////////////////////////////////////////
// Opérations mathématiques standards
/////////////////////////////////////////////////////////////////////////




#define DOP1_gen(AA) extern Tab AA ## _i(const Tab &x);\
  extern Tab AA ## i(const Tab &x);\
  inline Tab AA(const Tab &x){retourne AA ## _i(x);}

#define DOP1(AA) DOP1_gen(AA)\
  template<typename T, entier ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x)\
  {retourne AA ## _i(x);}\
  template<typename T>\
  VecT<T> AA(const VecT<T> &x)\
  {retourne AA ## _i(x);}


#define DOP1C2R(AA) DOP1_gen(AA)\
  template<typename T, entier ndims>\
  TabT<T, ndims> AA(const TabT<std::complex<T>, ndims> &x)\
  {retourne AA ## _i(x);}\
  template<typename T, entier ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x)\
  NOECLIPSE(requires (!est_complexe<T>()))\
  {retourne AA ## _i(x);}


#define DOP2(AA) extern Tab AA ## _i(const Tab &x, const Tab &y);\
  extern Tab AA ## i(const Tab &x, const Tab &y);\
  inline Tab AA(const Tab &x, const Tab &y){retourne AA ## _i(x,y);}\
  template<typename T, entier ndims>\
  TabT<T, ndims> AA(const TabT<T, ndims> &x, const TabT<T, ndims> &y)\
  {retourne AA ## _i(x,y);}\
  template<typename T>\
  VecT<T> AA(const VecT<T> &x, const VecT<T> &y)\
  {retourne AA ## _i(x,y);}

DOP1(sqrt)
DOP1(square)
DOP1(sin)
DOP1(cos)
DOP1(log)
DOP1(log10)
DOP1(exp)
DOP1(signe)
// complex -> float
DOP1C2R(abs)
DOP1C2R(abs2)
DOP1C2R(arg)

DOP2(pow)



template<typename T, entier ndims>
  TabT<T, ndims> Tab::as() const
{
  TabT<T, ndims> res(rows(), cols());
  ((Tab &) res) = *this;
  res.raw_data  = (T *) res.rawptr();
  retourne res;
}

template<typename T, entier ndims>
  TabT<T, ndims> Tab::as_no_conv() const
{
  TabT<T, ndims> res;
  // TODO : check type ok
  res.impl      = this->impl;
  res.raw_data  = (T *) rawptr();
  retourne res;
}


template<typename T, entier ndims>
Tab::Tab(const TabT<T, ndims> &t)
{
  impl = t.impl;
}




inline bouléen TRI(const Tab &v, auto f)
{
  soit K     = v.tscalaire();
  soit nbits = v.nbits();

  si(K == ℝ)
  {
    si(nbits == 32)
    {
      using T = float;
      soit t = v.as_no_conv<T>();
      f(t);
      retourne oui;
    }
    sinon si(nbits == 64)
    {
      using T = double;
      soit t = v.as_no_conv<T>();
      f(t);
      retourne oui;
    }
  }
  retourne non;
}



inline bouléen TCI(const Tab &v, auto f)
{
  soit K     = v.tscalaire();
  soit nbits = v.nbits();

  si(K == ℂ)
  {
    si(nbits == 64)
    {
      using T = cfloat;
      soit t = v.as_no_conv<T>();
      f(t);
      retourne oui;
    }
    sinon si(nbits == 128)
    {
      using T = cdouble;
      soit t = v.as_no_conv<T>();
      f(t);
      retourne oui;
    }
  }
  retourne non;
}

inline bouléen TZI(const Tab &v, auto f)
{
  soit K     = v.tscalaire();
  soit nbits = v.nbits();

  si(K == ℤ)
  {
    si(nbits == 32)
    {
      using T = int32_t;
      soit t = v.as_no_conv<T>();
      f(t);
      retourne oui;
    }
  }
  retourne non;
}

inline bouléen TBI(const Tab &v, auto f)
{
  soit K     = v.tscalaire();

  si(K == B)
  {
    using T = char;
    soit t = v.as_no_conv<T>();
    f(t);
    retourne oui;
  }
  retourne non;
}


inline void TR(const Tab &v, auto f)
{
  si(!TRI(v, f))
    échec("Type réel attendu.");
}

inline void TC(const Tab &v, auto f)
{
  si(!TCI(v, f))
    échec("Type complexe attendu.");
}

inline void TB(const Tab &v, auto f)
{
  si(!TBI(v, f))
    échec("Type bouléen attendu.");
}


inline void TRZ(const Tab &v, auto f)
{
  si(!TRI(v, f) && !TZI(v, f))
    échec("Type réel ou entier attendu.");
}

inline void TRC(const Tab &v, auto f)
{
  si(!TRI(v, f) && !TCI(v, f))
    échec("Type réel ou complexe attendu.");
}

inline void TG(const Tab &v, auto f)
{
  si(!TRI(v, f) && !TCI(v, f) && !TZI(v, f) && !TBI(v, f))
    échec("Type invalide.");
}





}


ostream_formater(tsd::NTenseurDim)
ostream_formater(tsd::Tab)
ostream_formater(tsd::Scalaire)
ostream_formater(tsd::cfloat)
ostream_formater(tsd::cdouble)

template <typename T, entier ndims> struct fmt::formatter<tsd::TabT<T,ndims>> { \
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();} \
  template <typename FormatContext> \
  auto format(const tsd::TabT<T,ndims>& t, FormatContext& ctx) const -> decltype(ctx.out())  \
  { \
    std::ostringstream ss; \
    ss << t; \
    return fmt::format_to(ctx.out(), "{}", ss.str()); \
  } \
};

#endif


