#include "tsd/tsd.hpp"
#include "tsd/tableau.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <fmt/core.h>
//#include <fmt/ostream.h>
//#include <fmt/color.h>

#include <vector>
#include <cstdlib>
#include <memory>


using namespace std;
using Eigen::ArrayXi;

namespace tsd {




template<typename t>
  using EVecteur = Eigen::Array<t, Eigen::Dynamic, 1, Eigen::ColMajor>;




// Comment faire une évaluation efficace
// Par exemple de :
//  x.carre().abs() * 4 - y.abs2() + Vecteur::ones(n)

template<typename D>
  void mult_i(D &x, entier i)
{
  x *= i;
}

template<typename D>
  void add_i(D &x, entier i)
{
  x += i;
}




static const string types[] =
{
    "VALEUR", "SEGMENT", "BLOCK", "EXPR", "MAP"
};

struct Tab::Impl: enable_shared_from_this<Tab::Impl>
{
  // Pas forcément contenu peut être une expression !
  enum Type
  {
    VALEUR,
    SEGMENT,
    BLOCK,
    EXPR,
    MAP
  };

  static string type2str(Type t)
  {
    retourne ((entier) t < 5) ? types[(entier) t] : string("?");
  }


  Type type       = VALEUR;
  entier seg_i0 = 0;
  entier blk_r0 = 0, blk_c0 = 0;
  Scalaire Tscalaire = ℝ;
  entier nbits       = 32;
  void *vals         = nullptr;
  ArrayXi dims;
  sptr<Impl> enfant;


  Impl()
  {

  }

  ~Impl()
  {
    si((vals != nullptr) && (type == VALEUR))
    {
      free(vals);
    }
    vals = nullptr;
  }

  bouléen est_de_même_dimensions(const Tab &t) const
  {
    //si(dims.size() != t.impl->dims.size())
    //  retourne non;
    pour(auto i = 0; i < max(dims.size(), t.impl->dims.size()); i++)
    {
      si(i >= t.impl->dims.size())
      {
        si(dims(i) != 1)
          retourne non;
      }
      sinon si (i >= dims.size())
      {
        si(t.impl->dims(i) != 1)
          retourne non;
      }
      sinon si(dims(i) != t.impl->dims(i))
        retourne non;
    }

    retourne oui;
  }

  bouléen est_de_même_format(const Tab &t) const
  {
    si(!est_de_même_dimensions(t))
      retourne non;
    retourne (Tscalaire == t.impl->Tscalaire) && (nbits == t.impl->nbits);
  }

  inline void *get_ptr_at(entier r) const
  {
    retourne ((char *) vals) + dim_scalaire() * r;
  }

  inline void *get_ptr_at(entier r, entier c) const
  {
    si(type != BLOCK)
      retourne ((char *) vals) + dim_scalaire() * (r + c * dims(0));
    retourne enfant->get_ptr_at(blk_r0 + r, blk_c0 + c);
  }


  sptr<Impl> clone() const
  {
    sptr<Impl> res = make_shared<Impl>();
    res->type        = Impl::Type::VALEUR;
    res->Tscalaire   = Tscalaire;
    res->nbits       = nbits;
    res->dims        = dims;

    si(vals != nullptr)
    {
      soit ne = nb_éléments() * dim_scalaire();
      si(type != BLOCK)
      {
        res->vals = malloc(ne);
        memcpy(res->vals, vals, ne);
      }
      sinon
      {
        soit nre = dims(0) * dim_scalaire();
        soit nc = dims(1);
        res->vals = malloc(ne);
        for(auto col = 0; col < nc; col++)
          memcpy(res->get_ptr_at(0, col), get_ptr_at(0, col), nre);
      }
    }
    retourne res;
  }

  entier nb_éléments() const
  {
    retourne dims.prod();
  }



  entier dim_scalaire() const
  {
    retourne nbits / 8;
  }





  sptr<Impl> eval() const
  {
    si((type == VALEUR) || (type == MAP))
      retourne ((Impl *) this)->shared_from_this();

    soit res = make_shared<Impl>();
    res->type       = VALEUR;
    res->dims       = dims;
    res->Tscalaire  = Tscalaire;
    res->nbits      = nbits;

    soit ds = dim_scalaire();
    soit ne = nb_éléments();

    switch(type)
    {

    case SEGMENT:
    {
      soit e = enfant;
      si(e->type != VALEUR)
        e = e->eval();
      res->vals = malloc(ne * ds);
      memcpy(res->vals,
             ((char *) e->vals) + seg_i0 * ds,
             dims(0) * ds);
      break;
    }

    case BLOCK:
    {
      soit e = enfant;
      si(e->type != VALEUR)
        e = e->eval();

      res->vals = malloc(ne * ds);

      soit nre  = e->dims(0);
      soit iptr = ((char *) e->vals) + blk_r0 * ds + blk_c0 * ds * nre;
      soit optr = (char *) res->vals;

      pour(auto c = 0; c < dims(1); c++)
      {
        memcpy(optr + c * ds * dims(0),
              iptr + c * nre * ds, dims(0) * ds);
      }
      break;
    }

#   if 0
    case EXPO:
      //res->vals = enfants2[0]->vals.exp();
      break;
    case CONST_1:
    {
      res->vals = malloc(nb_éléments() * dim_scalaire());
      APPLIQUE_METHODE((*res), setOnes());
      break;
    }
    case CONST_0:
    {
      res->vals = malloc(nb_éléments() * dim_scalaire());
      APPLIQUE_METHODE((*res), setZero());
      break;
    }
    case MULT_INT:
    {
      tsd_assert(enfants2.size() == 1);

      *res = *(enfants2[0]);
      //res = enfants2[0]; // PB ICI, pas de copie, et modification

      // Doit convertir en entier
      //si(enfants2[0]->is_boolean)
        //res = enfants2[0]->cast_int();

      APPLIQUE_OP((*res), mult_i, this->s_int);
      break;
    }
    case ADD_INT:
    {
      tsd_assert(enfants2.size() == 1);

      res = enfants2[0];

      // Doit convertir en entier
      //si(enfants2[0]->is_boolean)
        //res = enfants2[0]->cast_int();

      APPLIQUE_OP((*res), add_i, this->s_int);
      break;
    }
#   endif
    default:
    {
      msg_erreur("Evaluation expression type invalide : {} ({})", (entier) type, type2str(type));
    }
    }
    retourne res;

    //->type = Impl::VALEUR;
    // Détermine une manière efficace d'évaleur l'arbre d'expression
    //  -> Opérations scalaires, VS réductions partielles VS opérations matricielles
    // Puis à la fin, stockage du tableau évalué pour ne pas avoir à le réévaluer si nécessaire

  }

  const Tab tab() const
  {
    Tab res;
    res.impl = ((Impl *) this)->shared_from_this();
    retourne res;
  }

  string chaine() const
  {
    soit t = eval();
    string s;

    string dim;
    entier ndims = t->dims.rows();
    si(ndims == 0)
      dim = "[]";
    //sinon si(ndims == 1)
    //  dim = sformat("{}", t->dims(0));
    //sinon si(ndilms)
    sinon
    {
      pour(auto i = 0; i < ndims; i++)
      {
        dim += sformat("{}", t->dims(i));
        si(i + 1 < ndims)
          dim += "x";
      }
    }

    s = sformat("{} {}[{}] = ", dim, Tscalaire, nbits);

    si(ndims >= 2)
      s += "\n";

    soit tmp = tab();
    TG(tmp, [&]<typename T, entier ndims2>(TabT<T,ndims2> &)
    {
      soit m1 = emap<T>(tmp);
      si(ndims == 1)
        s += sformat("{}", m1.transpose());
      sinon
        s += sformat("{}", m1);
    });

    retourne s;
  }
};




void *Tab::rawptr()
{
  /*si(impl->type == Tab::Impl::Type::SEGMENT)
  {
    impl->enfants[0]->evalue();
    retourne ((char *) impl->enfants[0]->vals) + impl->seg_i0 * impl->dim_scalaire();
  }*/

  /*si(impl->type != Tab::Impl::Type::VALEUR)
  {
    echec("rawptr(): le tableau est une expression.");
  }*/
  retourne impl->vals;
}

const void *Tab::rawptr() const
{
  /*si(impl->type == Tab::Impl::Type::SEGMENT)
  {
    impl->enfants[0]->evalue();
    retourne ((char *) impl->enfants[0]->vals) + impl->seg_i0 * impl->dim_scalaire();
  }

  si(impl->type != Tab::Impl::Type::VALEUR)
  {
    echec("rawptr(): le tableau est une expression.");
  }*/
  retourne impl->vals;
}

Tab format_as(const Tab &x)
{
  Tab y;

  y.impl              = make_shared<Tab::Impl>();
  y.impl->type        = Tab::Impl::Type::VALEUR;
  y.impl->Tscalaire   = x.impl->Tscalaire;
  y.impl->nbits       = x.impl->nbits;
  y.impl->dims        = x.impl->dims;
  y.impl->vals        = malloc(y.impl->nb_éléments() * y.impl->dim_scalaire());

  retourne y;
}



entier Tab::nelems()     const
{
  retourne impl->dims.prod();
}


bool Tab::est_reference() const
{
  if(!impl)
    return false;
  retourne impl->type != Impl::Type::VALEUR;
}


Tab Tab::head(entier n) const
{
  retourne segment(0, n);
}

Tab Tab::tail(entier n) const
{
  retourne segment(nelems() - n, n);
}

Tab Tab::col(entier num) const
{
  Tab res;
  res.impl = make_shared<Impl>();

  //msg("Tab::col() -> reference.");

  tsd_assert_msg(num <= cols(), "Tab::col({}): dépassement ({} colonnes).", num, cols());

  res.impl->type        = Impl::Type::SEGMENT;
  res.impl->Tscalaire   = impl->Tscalaire;
  res.impl->nbits       = impl->nbits;
  res.impl->seg_i0      = num * rows();
  res.impl->dims.resize(1);
  res.impl->dims(0)     = rows();
  res.impl->vals        = impl->get_ptr_at(0, num);
  res.impl->enfant     = impl;
  retourne res;
}

Tab Tab::segment(entier i0, entier n) const
{
  Tab res;
  res.impl = make_shared<Impl>();

  tsd_assert_msg(i0 + n <= rows(), "Tab::segment({},{}): dépassement ({} éléments).", i0, n, rows());

  res.impl->type        = Impl::Type::SEGMENT;
  res.impl->Tscalaire   = impl->Tscalaire;
  res.impl->nbits       = impl->nbits;
  res.impl->seg_i0      = i0;
  res.impl->dims.resize(1);
  res.impl->dims(0)     = n;
  res.impl->vals        = impl->get_ptr_at(i0);
  res.impl->enfant      = impl;
  retourne res;
}



Tab Tab::reshape(entier m, entier n) const
{
  tsd_assert_msg(m * n == nelems(), "Tab::reshape : la dim totale doit être identique.");
  Tab x = Tab(impl->Tscalaire, impl->nbits, m, n);
  memcpy(x.rawptr(), this->rawptr(), m * n * impl->dim_scalaire());
  retourne x;
}


Tab Tab::block(entier r0, entier nr, entier c0, entier nc) const
{
  Tab res;
  res.impl = make_shared<Impl>();
  res.impl->type = Impl::Type::BLOCK;
  res.impl->blk_r0 = r0;
  res.impl->blk_c0 = c0;
  res.impl->dims.resize(2);
  res.impl->dims(0) = nr;
  res.impl->dims(1) = nc;
  res.impl->enfant = impl;
  res.impl->vals = impl->get_ptr_at(r0, c0);
  retourne res;
}

bool Tab::est_de_même_dimensions(const Tab &t) const
{
  retourne impl->est_de_même_dimensions(t);
}

void Tab::setConstant(double c)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1.setConstant((T) c);
  });
}

Tab Tab::reverse_rows() const
{
  Tab res = format_as(*this);
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(res);
    m2 = m1.rowwise().reverse();
  });
  retourne res;
}

Tab Tab::reverse() const
{
  Tab res = format_as(*this);
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(res);
    m2 = m1.reverse();
  });
  retourne res;
}

void Tab::reverseInPlace()
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1.reverseInPlace();
  });
}

double Tab::minCoeff(entier *i, entier *j) const
{
  double res = 0;
  entier i1 = 0, j1 = 0;

  if(nelems() == 0)
  {
    msg_avert("Tab::minCoeff(): tableau vide.");
  }
  else
  {
    TRZ(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
    {
      soit m1 = emap<T>(*this);
      res = m1.minCoeff(&i1,&j1);
    });
  }

  si(i != nullptr)
    *i = i1;
  si(j != nullptr)
    *j = j1;

  retourne res;
}

double Tab::maxCoeff(entier *i, entier *j) const
{
  double res = 0;
  entier i1 = 0, j1 = 0;

  if(nelems() == 0)
  {
    msg_avert("Tab::maxCoeff(): tableau vide.");
  }
  else
  {
    TRZ(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
    {
      soit m1 = emap<T>(*this);
      res = m1.maxCoeff(&i1,&j1);
    });
  }

  si(i != nullptr)
    *i = i1;
  si(j != nullptr)
    *j = j1;

  retourne res;
}

Scalaire Tab::tscalaire() const
{
  retourne impl->Tscalaire;
}

entier Tab::dim_scalaire() const
{
  retourne impl->dim_scalaire();
}

entier Tab::nbits() const
{
  retourne impl->nbits;
}

Tab Tab::clone() const
{
  Tab res;
  res.impl = impl->clone();
  retourne res;
}

Tab Tab::même_format() const
{
  retourne format_as(*this);
}

Tab::Tab()
{
  impl = make_shared<Impl>();
}


void setup(Tab::Impl &impl, entier n, entier m)
{
  si(m > 1)
  {
    impl.dims.resize(2);
    impl.dims(0) = n;
    impl.dims(1) = m;
  }
  sinon
  {
    impl.dims.resize(1);
    impl.dims(0) = n;
  }
  si(impl.nb_éléments() > 0)
    impl.vals = malloc(impl.nb_éléments() * impl.dim_scalaire());
  sinon
    impl.vals = nullptr;
}


void Tab::resize(entier n, entier m)
{
  si((rows() == n) && (cols() == m))
    retourne;
  //*this = Tab(impl->Tscalaire, impl->nbits, n, m);

  si((impl->vals != nullptr) && (impl->type == Impl::VALEUR))
    free(impl->vals);

  setup(*impl, n, m);
}

Tab::Tab(Scalaire s, entier reso, entier n, entier m)
{
  impl = make_shared<Impl>();
  impl->type      = Impl::VALEUR;
  impl->Tscalaire = s;
  impl->nbits     = reso;

  setup(*impl, n, m);
}

Tab Tab::map(Scalaire s, entier reso, entier n, entier m, void *data)
{
  Tab res;
  res.impl = make_shared<Impl>();
  res.impl->type      = Impl::MAP;
  res.impl->Tscalaire = s;
  res.impl->nbits     = reso;

  si(m > 1)
  {
    res.impl->dims.resize(2);
    res.impl->dims(0) = n;
    res.impl->dims(1) = m;
  }
  sinon
  {
    res.impl->dims.resize(1);
    res.impl->dims(0) = n;
  }

  res.impl->vals = data;

  return res;
}

entier Tab::rows() const
{
  si(impl->dims.size() == 0)
    retourne 0;
  retourne impl->dims(0);
}

entier Tab::cols() const
{
  si(impl->dims.size() < 2)
    retourne (rows() == 0) ? 0 : 1;
  retourne impl->dims(1);
}




Tab Tab::ones(Scalaire s, entier reso, entier n, entier m)
{
  Tab t{s, reso, n, m};

  TG(t, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    emap<T>(t).setOnes();
  });

  retourne t;
}

void Tab::setZero(entier n, entier m)
{
  si((n >= 0) || (impl->dims.rows() == 0))
    resize(n, m);
  memset(rawptr(), 0, nelems() * impl->dim_scalaire());
}

void Tab::set_ones(entier n, entier m)
{
  si((n >= 0) || (impl->dims.rows() == 0))
    resize(n, m);
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    emap<T>(*this).setOnes();
  });
}

void Tab::conservativeResize(entier n, entier m)
{
  tsd_assert(impl->type == Impl::VALEUR);

  soit l1 = n * m * impl->dim_scalaire();
  soit l2 = nelems() * impl->dim_scalaire();
  soit optr = impl->vals;

  setup(*impl, n, m);

  memcpy(impl->vals, optr, min(l1, l2));

  free(optr);
}

Tab Tab::transpose_int() const
{
  Tab y;
  y.impl->Tscalaire   = impl->Tscalaire;
  y.impl->nbits       = impl->nbits;


  si(impl->dims.size() == 1)
    setup(*(y.impl), 1, impl->dims(0));
  sinon
    setup(*(y.impl), impl->dims(1), impl->dims(0));


  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>>((T*)this->rawptr(), rows(), cols());
    soit m2 = Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>>((T*)y.rawptr(), cols(), rows());
    m2 = m1.transpose();
  });
  retourne y;
}

Tab Tab::conjugate() const
{
  Tab res = format_as(*this);
  TC(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>>((T*)this->rawptr(), rows(), cols());
    soit m2 = Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>>((T*)res.rawptr(), rows(), cols());
    m2 = m1.conjugate();
  });
  retourne res;
}

template<typename T>
auto emap(auto &x)
{
  retourne Eigen::Map<Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>((T*)x.rawptr(), x.rows(), x.cols());
}

template<typename T>
auto emap_mat(auto &x)
{
  retourne Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>((T*)x.rawptr(), x.rows(), x.cols());
}

template<typename T>
auto cemap(auto &x)
{
  retourne Eigen::Map<const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>((T*)x.rawptr(), x.rows(), x.cols());
}


Tab Tab::eye(Scalaire s, entier reso, entier n)
{
  soit x = Tab::zeros(s, reso, n, n);

  TG(x, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    pour(auto i = 0; i < n; i++)
      xt(i,i) = (T) 1;
  });

  //pour(auto i = 0; i < n; i++)
    //x(i,i).set(1);
  retourne x;
}

Tab Tab::zeros(Scalaire s, entier reso, entier n, entier m)
{
  Tab t(s, reso, n, m);
  t.setZero(n, m);
  retourne t;
}


void Tab::evaluer() const
{
}

Tab Tab::eval() const
{
  /*Tab res;
  res.impl = impl->eval();
  retourne res;*/
  return clone();
}





Tab Tab::constant(Scalaire s, entier reso, entier n, entier m, double valeur)
{
  Tab x(s, reso, n, m);

  TG(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(x);
    m1.setConstant(valeur);
  });
  retourne x;
}

Tab Tab::random(Scalaire s, entier reso, entier n, entier m)
{
  Tab x(s, reso, n, m);

  TG(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(x);
    m1.setRandom();
  });
  retourne x;
}

Tab Tab::matprod(const Tab &t) const
{
  Tab y(impl->Tscalaire, impl->nbits, rows(), t.cols());

  tsd_assert_msg(impl->Tscalaire == t.impl->Tscalaire, "Tab::Produit matriciel : types incompatibles ({} et {})", impl->Tscalaire, t.impl->Tscalaire);

  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(t);
    soit m3 = emap<T>(y);

    m3 = (m1.matrix() * m2.matrix()).array();
  });
  retourne y;
}


Tab arg_i(const Tab &x)
{
  si(x.impl->Tscalaire != ℂ)
  {
    echec("arg(): type complexe attendu.");
    retourne {};
  }

  Tab y(ℝ, x.impl->nbits / 2, x.rows(), x.cols());

  TC(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(x);
    soit m2 = emap<typename T::value_type>(y);

    m2 = m1.arg();
  });
  retourne y;
}


template<typename T1, typename T2>
  T1 conversion_scalaire(const T2 &t2)
{
  si constexpr (est_complexe<T1>())
  {
    si constexpr (est_complexe<T2>())
    {
      // Source et destination complexes
      retourne (T1) t2;
    }
    sinon
    {
      // Source réelle, destination complexe
      retourne (T1) t2;
    }
  }
  sinon si constexpr(est_complexe<T2>())
  {
    si constexpr (est_complexe<T1>())
    {
      // Source complexe, destination complexe
      retourne (T1) t2;
    }
    sinon
    {
      // Source complexe, destination réelle
      retourne (T1) real(t2);
    }
  }
  sinon
  {
    // Aucun des deux n'est complexe, simple cast
    retourne (T1) t2;
  }
}

void Tab::dump_infos() const
{
  msg("Infos pour tab: type = {}, dims={}, tscal={}, nbits={}, a enfant={}",
      Impl::type2str(impl->type),
      impl->dims.transpose(),
      impl->Tscalaire,
      impl->nbits,
      impl->enfant ? "oui" : "non");
}

void Tab::copie(const Tab &source)
{
  soit &i1 = *impl,
       &i2 = *(source.impl);

  si(i1.est_de_même_format(source))
  {
    //msg("Tab::copie -> même format.");
    soit ds = i1.dim_scalaire();
    // Attention : la source comme la destination peuvent être des segments !
    soit ne = i1.nb_éléments();
    tsd_assert(ne == i2.nb_éléments());

    si((i1.type != Impl::BLOCK) && ((i2.type != Impl::BLOCK)))
    {
      //msg("  -> copie non bloc.");
      memcpy(i1.vals, i2.vals, ne * ds);
    }
    // BLOC <- CONT
    // CONT <- BLOC
    // BLOC <- BLOC
    sinon
    {
      //msg("  -> copie bloc.");
      soit nr = i1.dims(0), nc = i1.dims(1);

      for(auto col = 0; col < nc; col++)
        memcpy(i1.get_ptr_at(0, col), i2.get_ptr_at(0, col), nr * ds);
      //echec("TODO: operator = pour block.");
    }
  }
  sinon si(i1.est_de_même_dimensions(source))
  {
    //msg("Tab::copie -> même dims.");

    // Idem, il faut juste convertir les types...
    // Attention : les cas de figure se croisent.

    TG(*this, [&]<typename T1, entier ndims1>(TabT<T1,ndims1> &)
    {
      TG(source, [&]<typename T2, entier ndims2>(TabT<T2,ndims2> &)
      {
        soit m1 = emap<T1>(*this);
        soit m2 = emap<T2>(source);

        pour(auto c = 0; c < cols(); c++)
        {
          pour(auto r = 0; r < rows(); r++)
          {
            m1(r,c) = conversion_scalaire<T1,T2>(m2(r,c));
            //m1(r,c) = (T1) m2(r,c);
          }
        }

      });
    });
  }
  sinon
  {
    msg("Dim self = {}, source = {}", impl->dims, source.impl->dims);
    echec("Tab::copie : source et destination non compatibles.");
  }
}

Tab &Tab::operator =(const Tab &source)
{
  //soit &i1 = *impl,
  //     &i2 = *(source.impl);
  // Deux cas possible:
  // (1) dimensions compatibles -> copie valeur
  // (2) dimensions non compatibles -> copie référence

  si(impl->est_de_même_dimensions(source))
  {
    copie(source);
  }
  sinon
  {
    //impl = source.impl;
    impl = source.impl->clone();
    //*this = source.clone();
  }
  retourne *this;
}

void essai()
{
  Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> a1(5,4);
  Eigen::Array<entier,Eigen::Dynamic,Eigen::Dynamic> a2(5,4);
  a1 = a2.cast<float>();
  a2 = a1.cast<entier>();
}

bouléen Tab::est_approx(const Tab &x, double tol) const
{
  bouléen res = non;

  si(x.tscalaire() != this->tscalaire())
  {
    echec("est_approx(): types scalaires différents ({} vs {}).", x.tscalaire(), tscalaire());
    retourne non;
  }

  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = cemap<T>(*this);
    soit m2 = cemap<T>(x);
    res = m1.isApprox(m2, tol);
  });
  retourne res;
}

bouléen Tab::est_nul() const
{
  retourne est_approx(0.0, 0.0);
}

bouléen Tab::est_approx(double x, double tol) const
{
  bouléen res = non;
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = cemap<T>(*this);
    res = m1.isApproxToConstant((T) x, tol);
  });
  retourne res;
}



/*void myessai(void)
{
  Eigen::ArrayXf m1;
  Eigen::Array<char, Eigen::Dynamic, 1> m2;

  //soit tmp = (m1 == 1.0f);
  soit tmp = (m1 == Eigen::ArrayXf::Ones(4));

  m2 = tmp;
}*/


#define OPERATEUR_COMP(XX)\
    Tab Tab::operator XX(double x) const\
    {\
      Tab y(B, 8, rows(), cols());\
      TR(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
      {\
        soit m1 = emap<T>(*this);\
        soit m2 = emap<bouléen>(y);\
        m2 = m1 XX /*x*/ (x * Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>::Ones(m1.rows(), m1.cols()));\
      });\
      retourne y;\
    }

OPERATEUR_COMP(>)
OPERATEUR_COMP(<)
OPERATEUR_COMP(==)
OPERATEUR_COMP(>=)
OPERATEUR_COMP(<=)


Tab Tab::cwiseMax(double x) const
{
  soit y = format_as(*this);

  TR(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(y);

    m2 = m1.cwiseMax(x);
  });
  retourne y;
}

Tab &Tab::operator *=(const entier &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 *= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator +=(const entier &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 += (T) x;
  });

  retourne *this;
}



Tab &Tab::operator *=(const float &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 *= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator *=(const double &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 *= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator *=(const cfloat &x)
{
  TC(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 *= (T) x;
  });

  retourne *this;
}


Tab &Tab::operator /=(const float &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 /= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator /=(const double &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 /= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator /=(const cfloat &x)
{
  TC(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 /= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator +=(const double &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 += (T) x;
  });

  retourne *this;
}

Tab &Tab::operator +=(const float &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 += (T) x;
  });

  retourne *this;
}

Tab &Tab::operator -=(const double &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 -= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator -=(const float &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 -= (T) x;
  });

  retourne *this;
}

Tab &Tab::operator -=(const cfloat &x)
{
  TC(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 -= (T) x;
  });

  retourne *this;
}




#define DEF_OP_SCAL(MD,TYPE,OP) Tab Tab::operator OP(const TYPE &x) const\
{\
  Tab y = format_as(*this);\
  MD(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(*this);\
    soit m2 = emap<T>(y);\
    m2 = m1 OP (T) x;\
  });\
  retourne y;\
}

DEF_OP_SCAL(TG, float,   -)
DEF_OP_SCAL(TG, double,  -)
DEF_OP_SCAL(TG, entier,  -)
DEF_OP_SCAL(TC, cfloat,  -)
DEF_OP_SCAL(TC, cdouble, -)

DEF_OP_SCAL(TG, float,   +)
DEF_OP_SCAL(TG, double,  +)
DEF_OP_SCAL(TG, entier,  +)
DEF_OP_SCAL(TC, cfloat,  +)
DEF_OP_SCAL(TC, cdouble, +)

DEF_OP_SCAL(TG, float,   *)
DEF_OP_SCAL(TG, double,  *)
DEF_OP_SCAL(TG, entier,  *)
DEF_OP_SCAL(TC, cfloat,  *)
DEF_OP_SCAL(TC, cdouble, *)

DEF_OP_SCAL(TG, float,   /)
DEF_OP_SCAL(TG, double,  /)
DEF_OP_SCAL(TG, entier,  /)
DEF_OP_SCAL(TC, cfloat,  /)
DEF_OP_SCAL(TC, cdouble, /)


#define DEF_OP_TAB(MD,OP) Tab Tab::operator OP(const Tab &x) const\
{\
  Tab y = format_as(*this);\
  MD(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(*this);\
    soit m2 = emap<T>(x);\
    soit m3 = emap<T>(y);\
    si((x.impl->Tscalaire != impl->Tscalaire))\
    {\
      echec("Tab::Opérateur {} : types incompatibles ({} et {})", #OP, x.impl->Tscalaire, impl->Tscalaire);\
    }\
    si((!x.impl->dims.isApprox(impl->dims)))\
    {\
      echec("Tab::Opérateur {} : dimensions incompatibles ({} et {})", #OP, x.impl->dims.transpose(), impl->dims.transpose());\
    }\
    m3 = m1 OP m2;\
  });\
  retourne y;\
}

#define DEF_OP_TAB2(MD,OP,OP2) Tab Tab::OP2(const Tab &x) const\
{\
  Tab y = format_as(*this);\
  MD(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(*this);\
    soit m2 = emap<T>(x);\
    soit m3 = emap<T>(y);\
    si((x.impl->Tscalaire != impl->Tscalaire))\
    {\
      echec("Tab::Opérateur {} : types incompatibles ({} et {})", #OP, x.impl->Tscalaire, impl->Tscalaire);\
    }\
    m3 = m1 OP m2;\
  });\
  retourne y;\
}

DEF_OP_TAB(TG,+)
DEF_OP_TAB2(TG,-,operator_minus)
DEF_OP_TAB(TG,*)
DEF_OP_TAB(TG,/)

Tab &Tab::operator -=(const entier &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    m1 -= (T) x;
  });

  //Eigen::VectorXf x1,x2;
  //x1.operator *()

  retourne *this;
}


Tab &Tab::operator +=(const Tab &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = emap<T>(*this);
    soit m2 = cemap<T>(x);
    m1 += m2;
  });

  retourne *this;
}

Tab &Tab::operator -=(const Tab &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = emap<T>(*this);
    soit m2 = cemap<T>(x);
    m1 -= m2;
  });

  retourne *this;
}

Tab &Tab::operator *=(const Tab &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = emap<T>(*this);
    soit m2 = cemap<T>(x);
    m1 *= m2;
  });

  retourne *this;
}

Tab &Tab::operator /=(const Tab &x)
{
  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &xt)
  {
    soit m1 = emap<T>(*this);
    soit m2 = cemap<T>(x);
    m1 /= m2;
  });

  retourne *this;
}

Tab Tab::inverse_matricielle() const
{
  soit y = format_as(*this);

  TRC(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap_mat<T>(*this);
    soit m2 = emap_mat<T>(y);
    m2 = m1.inverse();
  });

  retourne y;
}

Tab Tab::inv() const
{
  soit y = format_as(*this);

  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(y);
    m2 = ((T) 1) / m1;
  });

  retourne y;
}

Tab Tab::operator-( ) const
{
  soit y = format_as(*this);

  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    soit m2 = emap<T>(y);
    m2 = -m1;
  });

  retourne y;
}

bouléen Tab::hasNaN() const
{
  soit res = non;

  TG(*this, [&]<typename T, entier ndims>(TabT<T,ndims> &)
  {
    soit m1 = emap<T>(*this);
    res = m1.hasNaN();
  });

  retourne res;
}






ostream& operator<<(ostream& s, const Tab &t)
{
  s  << t.impl->chaine();
  retourne s;
}

std::ostream& operator<<(std::ostream& os, const Scalaire &s)
{
  switch(s)
  {
  case ℝ:
    os << "ℝ"; break;
  case ℂ:
    os << "ℂ"; break;
  case ℕ:
    os << "ℤ"; break;
  case B:
    os << "𝔹"; break;
  default:
    os << "?";
  }
  retourne os;
}




#define EIGEN_UNARY_OP(OP,OPE)\
Tab OP(const Tab &x)\
{\
  Tab y = format_as(x);\
  \
  TG(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(x);\
    soit m2 = emap<T>(y);\
    m2 = m1.OPE();\
  });\
  retourne y;\
}

#define EIGEN_UNARY_OP_C2R(OP,OPE)\
Tab OP(const Tab &x)\
{\
  Tab y;\
  \
  TCI(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(x);\
    y = Tab(ℝ, x.impl->nbits/2, x.rows(), x.cols());\
    soit m2 = emap<typename T::value_type>(y);\
    m2 = m1.OPE();\
  });\
  TRI(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
  {\
    soit m1 = emap<T>(x);\
    y = format_as(x);\
    soit m2 = emap<T>(y);\
    m2 = m1.OPE();\
  });\
  retourne y;\
}




#define EIGEN_BINARY_OP(OP,OPE)\
Tab OP(const Tab &x, const Tab &y)\
{\
  Tab z = format_as(x);\
\
TG(x, [&]<typename T, entier ndims>(TabT<T,ndims> &)\
{\
  soit m1 = emap<T>(x);\
  soit m2 = emap<T>(y);\
  soit m3 = emap<T>(z);\
  m3 = m1.OPE(m2);\
});\
\
  retourne z;\
}


EIGEN_UNARY_OP(sqrt_i,sqrt)
EIGEN_UNARY_OP(sin_i,sin)
EIGEN_UNARY_OP(cos_i,cos)
EIGEN_UNARY_OP(square_i,square)
EIGEN_UNARY_OP(log_i,log)
EIGEN_UNARY_OP(log10_i,log10)
EIGEN_UNARY_OP(exp_i,exp)
EIGEN_UNARY_OP(signe_i,sign)

EIGEN_UNARY_OP_C2R(abs_i,abs)
EIGEN_UNARY_OP_C2R(abs2_i,abs2)


EIGEN_BINARY_OP(pow_i,pow)




}
