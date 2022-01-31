#include "tsd/filtrage/spline.hpp"
#include "tsd/filtrage.hpp"

namespace tsd::filtrage {


  template<typename T>
  struct InterpolateurSinc: InterpolateurRIF<T>
  {
    InterpolateurSincConfig config;
    ArrayXXf lut;

    ArrayXf coefs(float phase)
    {
      tsd_assert((phase >= 0) && (phase <= 1));
      int lut_index = phase * config.nphases;
      tsd_assert(lut_index >= 0);
      tsd_assert(lut_index <= config.nphases);
      return lut.col(lut_index);
    }

    ArrayXf coefs_calcule(float τ)
    {
      auto nc = config.ncoefs;
      //assert((tau >= 0) && (tau <= 1+1e-7));

      ArrayXf h = ArrayXf::Zero(nc);
      for(auto i = 0; i < nc; i++)
        h(i) = tsd::filtrage::sinc2(i-nc/2-τ, config.fcut);

      // PB : il faut décaler la fenêtre de -tau !!
      //ArrayXf fen = fenetre(type_fenetre, npts);

      if(config.fenetre == "hn")
      {
        float a = 0.5, b = 0.25;

        ArrayXf t = linspace(-nc/2,(nc-1)/2,nc);
        t -= τ;
        t *= (2*π/nc);
        ArrayXf r1 = t.cos();
        ArrayXf r2 = a + 2 * b * r1;
        h *= r2;
      }
      return h;
    }

    InterpolateurSinc(const InterpolateurSincConfig &config)
    {
      this->nom           = format("sinc - ncoefs={}, nphases={}, fcut={}, fen={}",
          config.ncoefs, config.nphases, config.fcut, config.fenetre);
      this->config        = config;
      this->delais        = 0.5 * config.ncoefs;
      this->K             = config.ncoefs;

      verifie_frequence_normalisee(config.fcut, "interpolateur sinc");

      lut.resize(config.ncoefs, config.nphases + 1);
      for(auto i = 0; i <= config.nphases; i++)
        lut.col(i) = coefs_calcule((1.0 * i) / config.nphases);
    }
  };


  template<typename T>
  struct InterpolateurCSpline: InterpolateurRIF<T>
  {
    ArrayXf coefs(float phase)
    {
      int lut_index = phase * 256;
      tsd_assert((lut_index >= 0) && (lut_index <= 256));
      return lut.col(lut_index);
    }

    /** @param c Tension parameter (0 : Catmull-Rom spline)
     *  @param n Number of delayed version of the filter (time resolution).*/
    InterpolateurCSpline(unsigned int n = 256, float c = 0)
    {
      this->nom     = "cspline";
      this->K       = 4;
      this->delais  = 1.5;
      tsd::filtrage::cspline_calc_lut(lut, n, c);
    }
    ArrayXXf lut;
  };

  template<typename T>
  struct InterpolateurLineaire: InterpolateurRIF<T>
  {
    ArrayXf coefs(float τ)
    {
      ArrayXf h(2);
      h(0) = 1-τ;
      h(1) = τ;
      return h;
    }

    InterpolateurLineaire()
    {
      this->nom     = "linéaire";
      this->K       = 2; // needs 2 samples
      this->delais  = 0.5;
    }
  };

  template<typename T>
  struct InterpolateurLagrange: InterpolateurRIF<T>
  {
    int d;

    ArrayXf coefs(float τ)
    {
      ArrayXf h(d+1);

      // Points: 0 1 2 ... d
      // interpole à (d-1)/2 + mu
      // Exemple: d = 1
      // interpole à mu
      float t = ((d-1.0f) / 2) + τ;

      for(auto j = 0; j <= d; j++)
      {
        float p = 1.0f;
        for(auto k = 0; k <= d; k++)
        {
          if(k != j)
            p *= (t - k) / (j - k);
        }
        h(j) = p;
      }
      return h;
    }

    InterpolateurLagrange(unsigned int d)
    {
      // Pour un polynome de degré d,
      // on peut passer par d+1 point
      this->nom     = fmt::format("Lagrange degré {}", d);
      this->K       = d+1; // needs d+1 samples
      this->d       = d;
      this->delais  = 0.5f * d;//1 + d - ((d-1) / 2);
    }

  };


  template<typename T>
  sptr<InterpolateurRIF<T>> itrp_cspline()
  {
    return std::make_shared<InterpolateurCSpline<T>>();
  }

  template<typename T>
  sptr<InterpolateurRIF<T>> itrp_lineaire()
  {
    return std::make_shared<InterpolateurLineaire<T>>();
  }

  template<typename T>
  sptr<InterpolateurRIF<T>> itrp_lagrange(unsigned int degre)
  {
    return std::make_shared<InterpolateurLagrange<T>>(degre);
  }

  template<typename T>
  sptr<InterpolateurRIF<T>> itrp_sinc(const InterpolateurSincConfig &config)
  {
    return std::make_shared<InterpolateurSinc<T>>(config);
  }



template<typename T>
Vecteur<T> interp_lin(const ArrayXf &x, const Vecteur<T> &y, const ArrayXf &x2)
{
  int n = x.rows(), n2 = x2.rows();
  Vecteur<T> y2(n2);

  int j = 0;
  for(auto i = 0; i < n2; i++)
  {
    if(x2(i) <= x(0))
      y2(i) = y(0);
    else if(x2(i) >= x(n-1))
      y2(i) = y(n-1);
    else
    {
      for(; x(j) < x2(i); j++)
        ;

      auto a  = (x2(i) - x(j-1)) / (x(j) - x(j-1));
      y2(i) = (1-a) * y(j-1) + a * y(j);

      // Premier élément j tel que x2(i) < x(j)
      //auto j = find(x2(i) < x);

      //if(j.empty()) // x2 > tous les x(j)
      //  y2(i) = y(y.rows() - 1);
      //else if(j[0] == 0) // x2 < x(0)
      //  y2(i) = y(0);
      //else
      //{
      //auto a  = (x2(i) - x(j[0]-1)) / (x(j[0]) - x(j[0]-1));
      //y2(i) = (1-a) * y(j[0]-1) + a * y(j[0]);
      //}
    }
  }
  return y2;
}



/** @brief Spline interpolation
 *
 * Families of cubic spline interpolation:
 *    - Natural cubic splines: fixing f''(ti) = f''(ti+1) = 0
 *    - Clamped cubic spline: fixing f'(ti) = ki
 *      - Cardinal (Hermit) splines: Setting f'(ti) = (1-c) (pi+1 - pi-1) / (ti+1 - ti-1)
 *        - Catmull rom splines: with c = 0.
 *        - Zero-tangents: c = 1
 *      - Kochanek�Bartels spline: more general settings for f'(ti)
 *
 *  See http://en.wikipedia.org/wiki/Cubic_Hermite_spline */
template<typename T>
Vecteur<T> interp_cspline(const ArrayXf &x, const Vecteur<T> &y, const ArrayXf &x2)
{
  int i, n = x.size() - 1, n2 = x2.rows();
  Vecteur<T> B(n), C(n+1), D(n), Y(n), L(n+1), U(n+1), Z(n+1), R(n2);
  ArrayXf H(n);
  auto &A = y;

  H = diff(x);

  /*for(i = 0; i < n; i++)
  {
    H(i) = x(i+1) - x(i);
    tsd_assert_msg(H(i) > 0, "interp_cspline: vecteur x croissant attendu");
  }*/

  //for(i = 0; i <= n; i++)
    //A(i) = vy(i);
  //A = y;
  //auto dy = diff(y);
  //Y = (3 / H) * dy - (3 / H) * dy;

  for(i = 1; i < n; i++)
    Y(i) = ((T) (3.0 / H(i))) * (y(i+1) - y(i)) - ((T) (3.0 / H(i-1))) * (y(i) - y(i-1));

  L(0) = U(0) = Z(0) = 0;

  for(i = 1; i < n; i++)
  {
    L(i) = ((T) 2) * (x(i+1) - x(i-1)) - H(i-1) * U(i-1);
    U(i) = H(i) / L(i);
    Z(i) = (Y(i) - H(i-1) * Z(i-1)) / L(i);
  }
  L(n) = Z(n) = C(n) = 0;

  for(i = n-1; i >= 0; i--)
  {
    C(i) = Z(i) - U(i) * C(i+1);
    B(i) = (y(i+1) - y(i)) / H(i) - H(i) * (C(i+1) + ((T) 2) * C(i)) / ((T) 3.0);
    D(i) = (C(i+1) - C(i)) / (3 * H(i));
  }
  // Resample
  int id = 0;
  for(i = 0; i < x2.rows(); i++)
  {
    if(x2(i) <= x(0))
      R(i) = y(0);
    else if(x2(i) >= x(n))
      R(i) = y(n);
    else
    {
      for(; x(id) < x2(i); id++)
        ;
      // x2(i) compris dans [x(id-1) ; x(id)]
      float t = x2(i) - x(id-1);
      R(i) = A(id-1) + B(id-1) * t + C(id-1) * t * t + D(id-1) * t * t * t;
    }
  }
  return R;
}




template<typename T>
Vecteur<T> interp(const ArrayXf &x, const Vecteur<T> &y, const ArrayXf &x2, InterpOption mode)
{
  tsd_assert_msg(x.rows() == y.rows(), "interp: x et y doivent avoir le même nombre d'éléments (ici {} et {}).", x.rows(), y.rows());

  if(mode == InterpOption::LINEAIRE)
    return interp_lin(x, y, x2);
  return interp_cspline(x, y, x2);
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
 *  @param[out] h : h[0] : coefficient for p(0), h[1] for p'(0), h[2] for p(1) and h[3] for p'(1).
 *  @param t interpolation point between 0 and 1.0
 *
 *  p(t) = (2t^3-3t²+1) p0
 *       + (t^3-2t²+t)  m0
 *       + (-2t^3+3t²)  p1
 *       + (t^3-t²)     m1
 *
 * After factorization:
 *  p(t) = (1+2t)(1-t)² p0
 *       + t(t-1)²      m0
 *       + t²(3-2t)     p1
 *       + t²(t-1)      m1 */
// Coefficient sur (p0,m0,p1,m1)
void calc_cubic_spline(float h[4], float t)
{
  h[0] = (1 + 2 * t) * (t - 1)*(t - 1);
  h[1] = t * (t - 1) * (t - 1);
  h[2] = t * t * (3 - 2 * t);
  h[3] = t * t * (t - 1);
  //dsp_trace("cspline[%.2f]: %f %f %f %f", t, h[0], h[1], h[2], h[3]);
}

// Coefficients sur (p-1,p0,p1,p2)
void calc_filter(float g[4], float t, float c)
{
  float h[4];

  calc_cubic_spline(h, t);

  g[0] = - (1 - c) * h[1] / 2;
  g[1] = h[0] - (1-c) * h[3] / 2;
  g[2] = h[2] + (1-c) * h[1] / 2;
  g[3] = (1 - c) * h[3] / 2;

  //dsp_trace("cfilter[%.2f,%.2f]: %f %f %f %f", t, c, g[0], g[1], g[2], g[3]);
}

void cspline_calc_lut(ArrayXXf &lut, int n, float c)
{
  lut.resize(4, n+1);

  for(auto i = 0; i <= n; i++)
  {
    float mu = ((float) i) / n;
    float g[4];

    calc_filter(g, mu, c);

    for(auto j = 0; j < 4; j++)
      lut(j,i) = g[j];
  }
}



namespace hidden {

auto interp_float = interp<float>;
auto interp_double = interp<double>;
auto interp_cfloat = interp<cfloat>;

auto itrp_cspline_float = itrp_cspline<float>;
auto itrp_cspline_cfloat = itrp_cspline<cfloat>;
auto itrp_lagrange_float = itrp_lagrange<float>;
auto itrp_lagrange_cfloat = itrp_lagrange<cfloat>;
auto itrp_lineaire_float = itrp_lineaire<float>;
auto itrp_lineaire_cfloat = itrp_lineaire<cfloat>;
auto itrp_sinc_float = itrp_sinc<float>;
auto itrp_sinc_cfloat = itrp_sinc<cfloat>;

}

} // namespace



