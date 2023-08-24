#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <set>

using namespace std;
using namespace tsd::fourier;

namespace tsd::filtrage
{


template<typename T>
struct LigneARetard: FiltreGen<T>
{
  Vecteur<T> fenetre;
  entier ri = 0, retard = 0;

  LigneARetard(entier retard)
  {
    this->retard = retard;
    fenetre.setZero(retard);
    msg("ligne à retard : retard = {} échantillons.", retard);
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    si(retard == 0)
    {
      y = x;
      retourne;
    }

    soit n = x.rows();
    y.resize(n);

    pour(auto i = 0; i < n; i++)
    {
      y(i) = fenetre(ri);
      fenetre(ri) = x(i);
      ri = (ri + 1) % retard;
    }
  }

};

template<typename T>
  sptr<FiltreGen<T>> ligne_a_retard(entier retard)
{
  retourne make_shared<LigneARetard<T>>(retard);
}

template<typename T, typename Tc>
struct FiltreRIF: FiltreGen<T>
{
  Vecteur<T> fen;
  entier index = 0, K = 0;
  Vecteur<Tc> coefs;

  FiltreRIF(const Vecteur<Tc> &c)
  {
    this->coefs = c;
    K = c.rows();
    fen.setZero(K);
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    assertion(K > 0);
    soit n    = x.rows();
    soit iptr = x.data();
    soit optr = y.data();

    //msg("filtre rif: iptr={}, optr={}", (int64_t) iptr, (int64_t) optr);

    si(iptr != optr)
    {
      y.resize(n);
      optr = y.data();
    }

    pour(auto j = 0; j < n; j++)
    {
      soit cptr = coefs.data() + (K - 1);

      T somme = 0;

      // Remplace le plus ancien élément
      fen(index) = *iptr++;
      index = (index + 1) % K;
      soit fptr = fen.data() + index;

      // Number of samples at end of delay line
      soit K1 = K - index;
      // Number at begin of delay line
      soit K2 = K - K1;

      pour(auto i = 0; i < K1; i++)
        somme += *fptr++ * *cptr--;

      // Restart at beginning of delay line
      fptr = fen.data();
      pour(auto i = 0; i < K2; i++)
        somme += *fptr++ * *cptr--;

      *optr++ = somme;
    }
  }
};

template<typename T>
struct FiltreId: FiltreGen<T>
{
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    y = x;
  }
};

template<typename T>
  sptr<FiltreGen<T>> filtre_id()
{
  retourne make_shared<FiltreId<T>>();
}


template<typename T>
struct Decimateur: FiltreGen<T, T>
{
  entier R = 0, cnt = 0;
  Decimateur(entier R)
  {
    this->R = R;
  }
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();

    y.resize((n + R - 1 - cnt) / R);

    entier i, j = 0;
    pour(i = cnt; i < n; i += R)
    {
      si(j >= y.rows())
      {
        échec("j={},y.rows()={},i={},n={},cnt={},R={}", j, y.rows(), i, n, cnt, R);
      }

      y(j++) = x(i);
    }

    i -= R;
    assertion(j == y.rows());

    // si i = n-R : cnt = 0
    // si i = n-1 : cnt = R-1
    cnt = i - (n-R);

    // R - nombre d'échantillons non traités
    //cnt = (R - (n - 1 - i)) % R;
    //cnt = R - ((cnt + n) % R);
  }
};

template<typename T>
sptr<FiltreGen<T>> decimateur(entier R)
{
  retourne std::make_shared<Decimateur<T>>(R);
}

template<typename Tc, typename T>
sptr<FiltreGen<T>> filtre_rif(const Vecteur<Tc> &c)
{
  retourne make_shared<FiltreRIF<T,Tc>>(c);
}

template<typename T, typename Tc>
struct FiltreRII: FiltreGen<T>
{
  Vecteur<Tc> numer, denom;
  Vecteur<T> wndx, wndy;
  entier index = 0, index_y = 0;
  entier Kx = 0, Ky = 0;

  FiltreRII(const FRat<Tc> &h)
  {
    // Attention, ici, il faut passer en puissance de z-1
    soit h2 = h.eval_inv_z();
    soit num = h2.numer.vers_coefs().coefs;
    soit den = h2.denom.vers_coefs().coefs;

    init(num, den);
  }
  // numer = a0 + a1 z-1 + ...
  // denom = b0 + b1 z-1 + ...
  FiltreRII(const Vecteur<Tc> &numer, const Vecteur<Tc> &denom)
  {
    init(numer, denom);
  }
  void init(const Vecteur<Tc> &numer, const Vecteur<Tc> &denom)
  {
    index_y = index = 0;
    Kx = numer.rows();
    Ky = denom.rows() - 1;
    wndx.setZero(Kx);
    wndy.setZero(Ky);
    this->numer = numer;
    this->denom = denom;
  }
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();
    soit iptr = x.data();
    soit optr = y.data();
    si(optr != iptr)
    {
      y.resize(n);
      optr = y.data();
    }

    // Forme directe I : numérateur puis dénominateur

    // (1) Filtre non récursif
    pour(auto j = 0; j < n; j++)
    {
      soit nptr = numer.data() + (Kx - 1);
      T somme = 0;

      // Remplace le plus ancien élément
      wndx(index) = *iptr++;
      index = (index + 1) % Kx;
      soit wptr = wndx.data() + index;

      // Number of samples at end of delay line
      soit KK1 = Kx - index;
      // Number at begin of delay line
      soit KK2 = Kx - KK1;

      pour(auto i = 0; i < KK1; i++)
        somme += *wptr++ * *nptr--;

      // Redémarre au début de la ligne à retard
      wptr = wndx.data();
      pour(auto i = 0; i < KK2; i++)
        somme += *wptr++ * *nptr--;

      *optr++ = somme;
    }

    // (2) Filtre récursif
    optr = y.data();
    pour(auto j = 0; j < n; j++)
    {
      soit dptr = denom.data() + 1; // Premier élément ignoré
      T somme = *optr;
      // Remplace le plus ancien élément
      // wndy(index) = *ioptr;
      // index = (index + 1) % Ky;
      soit wptr = wndy.data() + index_y;

      // Number of samples at end of delay line
      soit KK1 = Ky - index_y;
      // Number at begin of delay line
      soit KK2 = Ky - KK1;

      pour(auto i = 0; i < KK1; i++)
        somme -= *wptr++ * *dptr++;

      // Redémarre au début de la ligne à retard
      wptr = wndy.data();
      pour(auto i = 0; i < KK2; i++)
        somme -= *wptr++ * *dptr++;

      *optr = somme / denom(0);

      // Remplace le plus ancien élément
      index_y = (index_y + Ky - 1) % Ky;
      wndy(index_y) = *optr++;
    }
  }


};

template<typename Tc, typename T>
  sptr<FiltreGen<T>> filtre_rii(const FRat<Tc> &h)
  {
    retourne make_shared<FiltreRII<T,Tc>>(h);
  }





struct SOISConfig
{
  Vecf coefs;
  RIIStructure structure = FormeDirecte2;
};



/** @brief Second Order IIR Section */
template<typename T = float, typename Tcoef = float, typename Ty = float>
struct SOIS: Filtre<T, T, SOISConfig>
{
  Ty x1 = (T) 0, x2 = (T) 0;
  Ty y0 = (T) 0, y1 = (T) 0, y2 = (T) 0;
  Tcoef b0 = 1, b1 = 1, b2 = 1, a0 = 1, a1 = 1, a2 = 1;
  bouléen premier_appel = oui;
  RIIStructure structure = RIIStructure::FormeDirecte2;

  /**       b0 + b1 z-1 + b2 z-2
   *  y/x = ---------------------
   *        a0 + a1 z-1 + a2 z-2
   */
  void configure(float b0, float b1, float b2, float a0, float a1, float a2)
  {
    assertion_msg(a0 != 0, "Section SOIS : a0 doit être non nul.");

    this->b0 = b0 / a0;
    this->b1 = b1 / a0;
    this->b2 = b2 / a0;
    this->a0 = 1.0;
    this->a1 = a1 / a0;
    this->a2 = a2 / a0;
    premier_appel = oui;
  }
  void configure_impl(const SOISConfig &config)
  {
    soit coefs = config.coefs;
    structure = config.structure;
    assertion_msg(coefs.rows() == 6,
              "Invalid number of coefficients pour SOIC initialization (6 required).");
    configure(coefs(0), coefs(1), coefs(2), coefs(3), coefs(4), coefs(5));
  }
  void configure(const FRat<float> &h)
  {
    assertion_msg(h.numer.coefs.rows() == 3,
              "Invalid number of coefficients pour SOIC initialization (3 required).");
    assertion_msg(h.denom.coefs.rows() == 3,
              "Invalid number of coefficients pour SOIC initialization (3 required).");
    configure(h.numer.coefs(0), h.numer.coefs(1), h.numer.coefs(2),
          h.denom.coefs(0), h.denom.coefs(1), h.denom.coefs(2));
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();

    si(n == 0)
      retourne;

    si(x.data() != y.data())
      y.resize(n);

    Ty x0;
    soit iptr = x.data();
    soit optr = y.data();

    si(premier_appel)
    {
      y0 = y1 = y2 = x2 = x1 = x(0);
      premier_appel = non;
    }

    si(structure == RIIStructure::FormeDirecte2)
    {
      pour(auto i = 0; i < n; i++)
      {
        x0 = *iptr++;

        soit d2 = x0 - a1 * y1 - a2 * y0;
        y2 = b0 * d2 + b1 * y1 + b2 * y0;

        *optr++ = (T) y2;

        y0 = y1;
        y1 = d2;
      }
    }
    sinon si(structure == RIIStructure::FormeDirecte1)
    {
      pour(auto i = 0; i < n; i++)
      {
        x0 = *iptr++;
        y0  = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = y0;
        x2 = x1;
        x1 = x0;
        *optr++ = (T) y0;
      }
    }
    sinon
      échec("SOIS : structure non implémentée.");
  }


};


/**       b0 + b1 z-1         z + b1/b0
 *  y/x = ----------- = b0 * ----------  =
 *        1 + a1 z-1            z + a1
 */
template<typename T, typename Tcoef>
struct RIIFoS: FiltreGen<T>
{
  Tcoef b0 = 0, b1 = 0, a1 = 0;
  T x1 = 0, y1 = 0;

  void configure(Tcoef b0, Tcoef b1, Tcoef a1)
  {
    this->b0 = b0;
    this->b1 = b1;
    this->a1 = a1;
    this->x1 = 0;
    this->y1 = 0;
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();
    si(x.data() != y.data())
      y.resize(n);
    pour(auto i = 0; i < n; i++)
    {
      soit x0 = x(i);
      y1 = -a1 * y1 + b0 * x0 + b1 * x1;
      y(i) = y1;
      x1 = x0;

      //msg("k={}, x={}, y={}", i, x(i), y(i));
    }
  }
};


template<typename T, typename Tcoef, typename Ty>
struct ChaineSOIS: FiltreGen<T>
{
  float gain = 1.0f;
  vector<SOIS<T, Tcoef, Ty>> sections;
  RIIFoS<T, Tcoef> rii1;
  bouléen avec_rii1 = non;

  ChaineSOIS(const FRat<cfloat> &h, RIIStructure structure)
  {
    soit [z,p] = h.roots();
    soit nz = z.rows(),
         np = p.rows();

    //msg("ChaineSOIS: nz={}, np={}", nz, np);

    // pour l'instant
    assertion_msg(nz == np,
        "SOIS : numérateur et dénominateur doivent avoir la même dimension (nz={}, np={})",
        nz, np);

    // Pool qui contient tous les index a priori, on vient puiser dedans
    set<entier> pool;
    pour(auto i = 0; i < nz; i++)
      pool.insert(i);

    // Groupe les racines et poles conjugés
    entier i;
    pour(i = 0; i + 1 < nz; i += 2)
    {
      soit k = *(pool.begin());
      pool.erase(pool.begin());

      // (z-(z(i)) * (z-z(i+1))
      // z² -(z(i)+z(i+1))*z + z(i)*z(i+1)

      // 1 -(z(i)+z(i+1))*z^-1 + z(i)*z(i+1)*z^-2

      SOIS<T, Tcoef, Ty> s;

      /**       b0 + b1 z-1 + b2 z-2
       *  y/x = ---------------------
       *        a0 + a1 z-1 + a2 z-2
       */
      // entier setup(float b0, float b1, float b2, float a0, float a1, float a2)

      float berr = 1e9;

      cfloat sz, pz, sp, pp;
      soit bj = pool.begin();

      pour(auto j = 0; j < nz; j++)
      {
        si(pool.count(j) == 0)
          continue;

        soit sz0 = -(z(j) + z(k)), pz0 = z(j) * z(k),
             sp0 = -(p(j) + p(k)), pp0 = p(j) * p(k);

        soit err = abs(sz0.imag()) + abs(pz0.imag()) + abs(sp0.imag()) + abs(pp0.imag());
        si(err < berr)
        {
          bj = pool.find(j);
          berr = err;
          sz = sz0;
          pz = pz0;
          sp = sp0;
          pp = pp0;
        }
      }

      pool.erase(bj);
      //msg("Section %d : err = %e", i, berr);

      si(berr > 1e-5)
      {
        msg_erreur("Factorisation SOIS : erreur = {}", berr);
      }

      soit coefs = Vecf::valeurs(
          {1.0f, sz.real(), pz.real(),
           1.0f, sp.real(), pp.real()});

      s.Configurable<SOISConfig>::configure({coefs, structure});

      sections.push_back(s);

    }

    pour(; i < nz; i++)
    {
      //msg("Ordre impair : insertion d'une section du premier ordre");
      assertion(pool.size() == 1);

      soit id = *(pool.begin());

      soit zer = z(id),
           pol = p(id);

      si(abs(pol.imag()) + abs(zer.imag()) > 1e-5f)
      {
        msg("pol.imag = {:e}, zer.imag = {:e}", pol.imag(), zer.imag());
        msg_erreur("Chaine RII: section du premier ordre restante : nombres imaginaires.");
      }

      /**       b0 + b1 z-1         z + b1/b0
       *  y/x = ----------- = b0 * ----------
       *        1 + a1 z-1            z + a1
       */

      soit a1 = - pol.real();
      soit b0 = h.numer.mlt.real() / h.denom.mlt.real();
      soit b1 = - zer.real() * b0;

      //msg("Cfg rii1: {}, {}, {}", b0, b1, a1);

      rii1.configure(b0, b1, a1);
      avec_rii1 = oui;
    }

    si(!avec_rii1)
      gain = h.numer.mlt.real() / h.denom.mlt.real();
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    y = x.clone();
    pour(auto &s: sections)
      y = s.FiltreGen<T>::step(y);
    si(avec_rii1)
      y = rii1.FiltreGen<T>::step(y);
    sinon
      y *= gain;
  }
};

template<typename T>
sptr<FiltreGen<T>> filtre_sois(const FRat<cfloat> &h, RIIStructure structure)
{
  retourne make_shared<ChaineSOIS<T,T,T>>(h, structure);
}

template<typename T>
sptr<FiltreGen<T>> filtre_sois(const FRat<float> &h, RIIStructure structure)
{
  soit &numer = h.numer, &denom = h.denom;

  FRat<cfloat> h2;
  h2.numer.mlt              = numer.mlt;
  h2.numer.vname            = numer.vname;
  h2.numer.coefs            = numer.coefs.clone();
  h2.denom.mlt              = denom.mlt;
  h2.denom.vname            = denom.vname;
  h2.denom.coefs            = denom.coefs.clone();

  si(!numer.mode_racines)
    h2.numer = Poly<cfloat>::from_roots(numer.roots())
                * numer.coefs(numer.coefs.rows() - 1);
  si(!denom.mode_racines)
    h2.denom = Poly<cfloat>::from_roots(denom.roots())
                * denom.coefs(denom.coefs.rows() - 1);


  retourne make_shared<ChaineSOIS<T,T,T>>(h2, structure);
}

template<typename T = float>
struct FiltreDC: FiltreGen<T>
{
  float α;
  T xp, yp;

  FiltreDC(float fc)
  {
    α = 1 - lexp_coef(Fréquence(fc));
    xp = yp = (T) 0;
  }
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();

    si(y.data() != x.data())
      y.resize(n);

    pour(auto i = 0; i < n; i++)
    {
      y(i) = (T) (α * ((x(i) - xp) + yp));
      xp = x(i);
      yp = y(i);
    }
  }
};



template<typename T = float, typename Tacc = float, typename Tinv = float>
struct MoyenneGlissante: FiltreGen<T>
{
  entier index = 0, K = 0;
  Tacc accu = 0;
  Vecteur<T> fenetre;
  // Inverse de la longueur
  Tinv K_inv;

  MoyenneGlissante(entier K)
  {
    this->K = K;
    fenetre.setZero(K);
    K_inv = (Tinv) (1.0 / ((double) K));
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();

    assertion(K > 0);

    si(y.data() != x.data())
      y.resize(n);

    pour(auto i = 0; i < n; i++)
    {
      accu += x(i);
      accu -= fenetre(index);
      fenetre(index) = x(i);
      index = (index + 1) % K;
      y(i) = ((T) accu) * K_inv;
    }
  }
};


// TODEL
#if 0
/** @brief First order exponential high-pass filter */
template<typename T = float, typename Tacc = float, typename Tgamma = float>
struct ExpHpFilter
{
  Tacc   acc = 0;
  // Coefficient du filtre
  Tgamma γ = 0;
  bouléen   init_fait = non;

  void setup(float γ)
  {
    this->γ = γ;
    reset();
  }

  void reset()
  {
    init_fait = non;
    acc  = 0;
  }

  /** @param[in] in   Input matrix.
   *  @param[out] out Output matrix (can be the same object as input si in place computation is needed). */
  void step(Vecteur<T> &out, const Vecteur<T> &in)
  {
    entier i, nsples = in.rows();

    si(nsples == 0)
      retourne;

    out.resize(nsples);

    si(!init_fait)
      acc = in(0); // Will automatically shift input data of n bits, si acc
                   // is fixed point.

    pour(i = 0; i < nsples; i++)
    {
      Tacc x = in(i);
      acc += γ * (x - acc);
      out(i) = (T) (x - acc);
    }

    init_fait = oui;
  }

};
#endif

/** @brief First order exponential low-pass filter
 *  Canonical form:
 *
 *  y_n+1 = (1-g) y_n + g x_n e.g. y_n+1 = y_n + g(x_n - y_n) */
template<typename T = float, typename Tacc = float, typename Tγ = float>
struct FiltreLExp : FiltreGen<T>
{
  Tacc acc     = 0;
  Tγ   γ       = 0;
  bouléen init = non;

  FiltreLExp(float γ)
  {
    this->γ = γ;
    init = non;
    acc  = 0;
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit n = x.rows();

    si(n == 0)
      retourne;

    y.resize(n);

    si(!init)
      acc = x(0); // Will automatically shift input data of n bits, si acc
                   // is fixed point.

    // acc_i = xpi + (1-g)*acc_{i-1}

    pour(auto i = 0; i < n; i++)
    {
      Tacc xi = x(i);
      acc += γ * (xi - acc);
      y(i) = (T) acc;
    }

    init = oui;
  }
};


template<typename T> sptr<FiltreGen<T>> filtre_lexp(float γ)
{
  retourne make_shared<FiltreLExp<T,T,T>>(γ);
}

template<typename T>
  sptr<FiltreGen<T>> filtre_dc(float fc)
{
  retourne make_shared<FiltreDC<T>>(fc);
}

template<typename T, typename Tacc>
  sptr<FiltreGen<T>> filtre_mg(entier n)
{
  retourne make_shared<MoyenneGlissante<T,Tacc,T>>(n);
}





// Instanciations explicites.
template sptr<FiltreGen<float>> filtre_lexp(float γ);
template sptr<FiltreGen<cfloat>> filtre_lexp(float γ);

template sptr<FiltreGen<float>> filtre_dc(float fc);
template sptr<FiltreGen<cfloat>> filtre_dc(float fc);

template sptr<FiltreGen<float>> filtre_rii(const FRat<float> &h);
template sptr<FiltreGen<cfloat>> filtre_rii(const FRat<cfloat> &h);
template sptr<FiltreGen<float>> ligne_a_retard<float>(entier retard);
template sptr<FiltreGen<cfloat>> ligne_a_retard<cfloat>(entier retard);
template sptr<FiltreGen<float>> filtre_mg<float,double>(entier n);
template sptr<FiltreGen<cfloat>> filtre_mg<cfloat,cdouble>(entier n);


template sptr<FiltreGen<float>> filtre_sois<float>(const FRat<cfloat> &coefs, RIIStructure structure);
template sptr<FiltreGen<float>> filtre_sois<float>(const FRat<float> &coefs, RIIStructure structure);
template sptr<FiltreGen<cfloat>> filtre_sois<cfloat>(const FRat<cfloat> &coefs, RIIStructure structure);
template sptr<FiltreGen<cfloat>> filtre_sois<cfloat>(const FRat<float> &coefs, RIIStructure structure);



namespace hidden {
soit filtre_id1  = filtre_id<float>;
soit filtre_id2  = filtre_id<cfloat>;

soit decimateur1  = decimateur<float>;
soit decimateur2  = decimateur<cfloat>;

soit filtre_rif1 = filtre_rif<float, float>;
soit filtre_rif2 = filtre_rif<float, cfloat>;
soit filtre_rif4 = filtre_rif<cfloat, cfloat>;




}


}

