#include "tsd/filtrage/spline.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/divers.hpp"

using namespace std;

namespace tsd::filtrage {


template<typename T>
struct InterpolateurSinc: InterpolateurRIF<T>
{
  InterpolateurSincConfig config;
  Tabf lut;

  Vecf coefs(float τ)
  {
    assertion_msg((τ >= 0) && (τ <= 1), "InterpolateurSinc::coefs(τ={}) : délais invalide.", τ);
    entier lut_index = τ * config.nphases;
    assertion((lut_index >= 0) && (lut_index <= config.nphases));
    retourne lut.col(lut_index);
  }

  Vecf coefs_calcule(float τ)
  {
    soit nc = config.ncoefs;
    soit h = Vecf::int_expr(nc, [&](entier i){retourne sinc(2 * config.fcut, i-nc/2-τ);});

    // Fenêtre de Hann, décalée d'un retard fractionnaire -τ
    si(config.fenetre == "hn")
    {
      soit a  = 0.5f, b = 0.25f;
      soit t  = (linspace(-nc/2,(nc-1)/2,nc) - τ) * (2*π/nc),
           r1 = cos(t),
           r2 = a + 2 * b * r1;
      h *= r2;
    }
    retourne h;
  }

  InterpolateurSinc(const InterpolateurSincConfig &config)
  {
    this->nom           = format("sinc - ncoefs={}, nphases={}, fcut={}, fen={}",
        config.ncoefs, config.nphases, config.fcut, config.fenetre);
    this->config        = config;
    this->K             = config.ncoefs;
    this->delais        = 0.5 * this->K;

    verifie_frequence_normalisee(config.fcut, "interpolateur sinc");

    lut.resize(this->K, config.nphases + 1);
    pour(auto i = 0; i <= config.nphases; i++)
      lut.col(i) = coefs_calcule((1.0 * i) / config.nphases);
  }
};


template<typename T>
struct InterpolateurCSpline: InterpolateurRIF<T>
{
  Tabf lut;

  Vecf coefs(float τ)
  {
    soit lut_index = (entier) (τ * 256);
    assertion((lut_index >= 0) && (lut_index <= 256));
    retourne lut.col(lut_index);
  }

  /** @param c Tension parameter (0 : Catmull-Rom spline)
   *  @param n Number of delayed version of the filter (time resolution).*/
  InterpolateurCSpline(entier n = 256, float c = 0)
  {
    this->nom     = "cspline";
    this->K       = 4;
    this->delais  = 1.5;
    this->lut     = cspline_calc_lut(n, c);
  }
};

template<typename T>
struct InterpolateurLineaire: InterpolateurRIF<T>
{
  Vecf coefs(float τ)
  {
    retourne Vecf::valeurs({1-τ, τ});
  }

  InterpolateurLineaire()
  {
    this->nom     = "linéaire";
    this->K       = 2; // requiert 2 échantillon
    this->delais  = 0.5;
  }
};

template<typename T>
struct InterpolateurLagrange: InterpolateurRIF<T>
{
  entier d;

  InterpolateurLagrange(entier d)
  {
    // pour un polynome de degré d,
    // on peut passer par d+1 point
    this->nom     = sformat("Lagrange degré {}", d);
    this->K       = d + 1; // nécessite d+1 échantillons
    this->d = d;
    this->delais  = 0.5 * d;
  }

  Vecf coefs(float τ)
  {
    Vecf h(d+1);

    // Points: 0 1 2 ... d
    // interpole à (d-1)/2 + mu
    // Exemple: d = 1
    // interpole à mu
    soit t = ((d - 1.0f) / 2) + τ;

    pour(auto j = 0; j <= d; j++)
    {
      soit p = 1.0f;
      pour(auto k = 0; k <= d; k++)
        si(k != j)
          p *= (t - k) / (j - k);
      h(j) = p;
    }
    retourne h;
  }
};


template<typename T>
sptr<InterpolateurRIF<T>> itrp_cspline()
{
  retourne make_shared<InterpolateurCSpline<T>>();
}

template<typename T>
sptr<InterpolateurRIF<T>> itrp_lineaire()
{
  retourne make_shared<InterpolateurLineaire<T>>();
}

template<typename T>
sptr<InterpolateurRIF<T>> itrp_lagrange(unsigned int degre)
{
  retourne make_shared<InterpolateurLagrange<T>>(degre);
}

template<typename T>
sptr<InterpolateurRIF<T>> itrp_sinc(const InterpolateurSincConfig &config)
{
  retourne make_shared<InterpolateurSinc<T>>(config);
}



// Interpolation linéaire, points non équidistants
template<typename T>
Vecteur<T> interp_lin(const Vecf &x, const Vecteur<T> &y, const Vecf &x2)
{
  soit n = x.rows(), n2 = x2.rows(), j = 0;
  Vecteur<T> y2(n2);

  pour(auto i = 0; i < n2; i++)
  {
    si(x2(i) <= x(0))
      y2(i) = y(0);
    sinon si(x2(i) >= x(n-1))
      y2(i) = y(n-1);
    sinon
    {
      pour(; x(j) < x2(i); j++)
        ;

      soit a  = (x2(i) - x(j-1)) / (x(j) - x(j-1));
      y2(i) = (1-a) * y(j-1) + a * y(j);
    }
  }
  retourne y2;
}



/** Interpolation spline
 *
 * Families of cubic spline interpolation:
 *    - Natural cubic splines: fixing f''(ti) = f''(ti+1) = 0
 *    - Clamped cubic spline: fixing f'(ti) = ki
 *      - Cardinal (Hermit) splines: Setting f'(ti) = (1-c) (pi+1 - pi-1) / (ti+1 - ti-1)
 *        - Catmull rom splines: with c = 0.
 *        - Zero-tangents: c = 1
 *      - Kochanek�Bartels spline: more general settings pour f'(ti)
 *
 *  Voir http://en.wikipedia.org/wiki/Cubic_Hermite_spline */
template<typename T>
Vecteur<T> interp_cspline(const Vecf &x, const Vecteur<T> &y, const Vecf &x2)
{
  entier i, n = x.rows() - 1, n2 = x2.rows();

  // (1) Calcul des coefficients des polynômes :
  //    Ak + Bk * t + Ck * t² + Dk * t^3;
  // Note : Ak = yk

  // Utiliser plutôt un vecteur de Poly<> ?

  Vecteur<T> B(n), C(n+1), D(n), Y(n), L(n+1), U(n+1), Z(n+1), R(n2);
  soit &A = y;

  soit H = diff(x);

  pour(i = 1; i < n; i++)
    Y(i) = ((T) (3.0 / H(i))) * (y(i+1) - y(i)) - ((T) (3.0 / H(i-1))) * (y(i) - y(i-1));

  L(0) = U(0) = Z(0) = 0;
  pour(i = 1; i < n; i++)
  {
    L(i) = ((T) 2) * (x(i+1) - x(i-1)) - H(i-1) * U(i-1);
    U(i) = H(i) / L(i);
    Z(i) = (Y(i) - H(i-1) * Z(i-1)) / L(i);
  }
  L(n) = Z(n) = C(n) = 0;

  pour(i = n-1; i >= 0; i--)
  {
    C(i) = Z(i) - U(i) * C(i+1);
    B(i) = (y(i+1) - y(i)) / H(i) - H(i) * (C(i+1) + ((T) 2) * C(i)) / ((T) 3);
    D(i) = (C(i+1) - C(i)) / (3 * H(i));
  }

  // Rééchantillonnage
  soit id = 0;
  pour(i = 0; i < n2; i++)
  {
    si(x2(i) <= x(0))
      R(i) = y(0);
    sinon si(x2(i) >= x(n))
      R(i) = y(n);
    sinon
    {
      pour(; x(id) < x2(i); id++)
        ;
      // x2(i) compris dans [x(id-1) ; x(id)]
      soit k = id - 1;
      soit t = x2(i) - x(k);
      R(i) = A(k) + B(k) * t + C(k) * t * t + D(k) * t * t * t;
    }
  }
  retourne R;
}




template<typename T>
Vecteur<T> interp(const Vecf &x, const Vecteur<T> &y, const Vecf &x2, InterpOption mode)
{
  assertion_msg(x.rows() == y.rows(),
      "interp: x et y doivent avoir le même nombre d'éléments (ici {} et {}).", x.rows(), y.rows());

  si(mode == InterpOption::LINEAIRE)
    retourne interp_lin(x, y, x2);

  retourne interp_cspline(x, y, x2);
}

/** @brief Compute the coefficients of an interpolation polynomial in p0,m0,p1,m1.
 *         Compute the third degree polynomial p(t)
 *         with the following constraints on the values of zeroth
 *         and first derivative at 0, 1:
 *           - p(0) = p0,
 *           - p'(0) = m0,
 *           - p(1) = p1,
 *           - p'(1) = m1.
 *
 *  @param[out] h : h[0] : coefficient pour p(0), h[1] pour p'(0), h[2] pour p(1) and h[3] pour p'(1).
 *  @param t interpolation point between 0 and 1.0
 *
 *  p(t) = (2t^3-3t²+1) p0
 *       + (t^3-2t²+t)  m0
 *       + (-2t^3+3t²)  p1
 *       + (t^3-t²)     m1
 *
 * Après factorisation:
 *  p(t) = (1+2t)(1-t)² p0
 *       + t(t-1)²      m0
 *       + t²(3-2t)     p1
 *       + t²(t-1)      m1 */
// Coefficient sur (p0,m0,p1,m1)
Vecf cspline_calc(float t)
{
  retourne Vecf::valeurs({
    (1 + 2 * t) * (t - 1) *(t - 1),
    t * (t - 1) * (t - 1),
    t * t * (3 - 2 * t),
    t * t * (t - 1)});
}

// Coefficients sur (p-1,p0,p1,p2)
Vecf cspline_filtre(float t, float c)
{
  soit h = cspline_calc(t);

  return Vecf::valeurs({
    - (1 - c) * h(1) / 2,
    h(0) - (1-c) * h(3) / 2,
    h(2) + (1-c) * h(1) / 2,
    (1 - c) * h(3) / 2});
}

Tabf cspline_calc_lut(entier n, float c)
{
  Tabf lut(4, n+1);
  pour(auto i = 0; i <= n; i++)
    lut.col(i) = cspline_filtre(((float) i) / n, c);
  retourne lut;
}



namespace hidden {

soit interp_float = interp<float>;
soit interp_double = interp<double>;
soit interp_cfloat = interp<cfloat>;

soit itrp_cspline_float = itrp_cspline<float>;
soit itrp_cspline_cfloat = itrp_cspline<cfloat>;
soit itrp_lagrange_float = itrp_lagrange<float>;
soit itrp_lagrange_cfloat = itrp_lagrange<cfloat>;
soit itrp_lineaire_float = itrp_lineaire<float>;
soit itrp_lineaire_cfloat = itrp_lineaire<cfloat>;
soit itrp_sinc_float = itrp_sinc<float>;
soit itrp_sinc_cfloat = itrp_sinc<cfloat>;

}

} // namespace



