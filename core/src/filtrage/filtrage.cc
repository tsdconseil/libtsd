#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <iostream>
#include <set>

using namespace std;


namespace tsd
{
  ostream& operator<<(ostream& os, const Poly<float> &p)
  {
    p.display_content(os);
    return os;
  }
  ostream& operator<<(ostream& os, const FRat<float> &f)
  {
    f.afficher(os);
    return os;
  }

}

namespace tsd::filtrage
{

  ArrayXf design_rif_prod(const ArrayXf &h1, const ArrayXf &h2)
  {
    int n1 = h1.rows();//, n2 = h2.rows();
    ArrayXf h2p = vconcat(h2, ArrayXf::Zero(n1-1));
    ArrayXf h3 = filtrer(h1, h2p);
    return h3;
  }

  void verifie_frequence_normalisee(float f, const string &msg)
  {
    string s = "F";
    if(!msg.empty())
      s = msg + " : f";
    if((f < 0) || (f > 0.5))
      echec("{}réquence normalisée attendue, f = {} (devrait être comprise entre 0 et 0.5).", s, f);
  }

  template<typename T>
  Tableau<T> forme_polyphase(const Eigen::Ref<const Vecteur<T>> x, unsigned int M)
  {
    auto n = x.rows();
    auto r = n % M;
    Vecteur<T> x2;
    if(r != 0)
    {
      x2 = vconcat(x, Vecteur<T>::Zero(M-r));
    }
    else
      x2 = x;
    n = x2.rows();

    // n / M = nb échantillons par canal
    auto X = Eigen::Map<Tableau<T>>(x2.data(), M, n / M);
    // Inversion des lignes
    return X.rowwise().reverse();
  }


  template<typename T>
  Vecteur<T> iforme_polyphase(const Eigen::Ref<const Tableau<T>> x)
  {
    return Eigen::Map<const Vecteur<T>>(x.data(), x.cols() * x.rows());
  }




  ArrayXcf hilbert_tfd(IArrayXf x)
  {
    auto N = x.rows();
    ArrayXcf X = tsd::fourier::fft(x);
    X.tail(N/2) = 0;
    return tsd::fourier::ifft(2 * X);
  }


  template<typename T>
  tuple<ArrayXf, ArrayXf> frphase(const FRat<T> &h, int npts)
  {
    auto fr = linspace(0, 0.5 - (0.5 / npts), npts);
    auto dfr = polar(2.0f * π * fr);

    ArrayXf xm(npts);
    for(auto i = 0; i < npts; i++)
      xm(i) = arg(h.horner(dfr(i)));

    xm = unwrap(xm, π);
    return make_tuple(fr, xm);
  }



  template<typename T> tuple<ArrayXf, ArrayXf> frgroup(const Vecteur<T> &h, int npts)
    {
      return frgroup(FRat<T>::rif(h), npts);
    }

  template<typename T> tuple<ArrayXf, ArrayXf> frmag(const Vecteur<T> &h, int npts)
  {
    Vecteur<T> h2 = Vecteur<T>::Zero(2*npts);
    h2.head(h.rows()) = h;
    ArrayXf fr = linspace(0, 0.5 - (0.5 / npts), npts);
    return {fr, tsd::fourier::fft(h2).abs().head(npts) * sqrt(2*npts)};
  }

  template<typename T> tuple<ArrayXf, ArrayXf> frphase(const Vecteur<T> &h, int npts)
  {
    return frphase(FRat<T>::rif(h), npts);
  }

  template<typename T>
  tuple<ArrayXf, ArrayXf> frgroup(const FRat<T> &h, int npts)
  {
    auto [fr, ϕ] = frphase(h, npts);

    // Pas de phase pour la dérivée
    float delta_ϕ = ϕ(1) - ϕ(0);//2 * π / npts;

    ArrayXf ϕp(npts);
    ϕp.head(npts-1) = diff(ϕ) / delta_ϕ;
    ϕp(npts-1)      = ϕp(npts-2);


    /*ϕp = ϕp.unaryExpr([](float x)
    {
      return wrap_pm_pi(x);
      //x = wrap_2pi(x);
      //auto res = fmod(x, 2*π_f);
      //if(res < 0)
      //  res += 2*π;
      //if(x > π)
        //x = 2*π - x;
      //return x;
    });*/

    return {fr, ϕp};
  }

  template<typename T>
  ArrayXcf repfreq(const FRat<T> &h, const ArrayXf &fr)
  {
    int n = fr.rows();
    ArrayXcf H(n);
    //auto h2 = FRat<float>::rif(h);
    ArrayXf pl = 2.0f * π_f * fr;
    ArrayXcf dfr = polar(pl);
    for(auto i = 0; i < n; i++)
      H(i) = h.horner(dfr(i));
    return H;
  }

  ArrayXcf repfreq(const ArrayXf &h, const ArrayXf &fr)
  {
    return repfreq(FRat<float>::rif(h), fr);
  }

  /*ArrayXcf repfreq(const ArrayXf &h, const ArrayXf &fr)
  {
    int n = fr.rows();
    ArrayXcf H(n);
    auto h2 = FRat<float>::rif(h);
    ArrayXf pl = 2.0f * π_f * fr;
    ArrayXcf dfr = polar(pl);
    for(auto i = 0; i < n; i++)
      H(i) = h2.horner(dfr(i));
    return H;
  }*/

  template<typename T>
  tuple<ArrayXf, ArrayXf> frmag(const FRat<T> &h, int npts)
  {
    ArrayXf fr = linspace(0, 0.5 - (0.5 / npts), npts);
    return{fr, repfreq(h, fr).abs()};
  }

  template<typename T>
  struct LigneARetard: FiltreGen<T>
  {
    Vecteur<T> fenetre;
    int ri = 0, retard = 0;

    LigneARetard(int retard)
    {
      ri = 0;
      this->retard = retard;
      fenetre.setZero(retard);

      msg("ligne à retard : retard = {} échantillons.", retard);
    }
    FRat<float> get_frat()
    {
      return FRat<float>::z_power(-retard);
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      if(retard == 0)
      {
        y = x;
        return;
      }

      int n = x.rows();
      y.resize(n);

      for(auto i = 0; i < n; i++)
      {
        y(i) = fenetre(ri);
        fenetre(ri) = x(i);
        ri = (ri + 1) % retard;
      }
    }

  };

  template<typename T>
    sptr<FiltreGen<T>> ligne_a_retard(int retard)
  {
    return make_shared<LigneARetard<T>>(retard);
  }

  template<typename T, typename Tc>
  struct FiltreRIF: FiltreGen<T>
  {
    Vecteur<T> wnd;
    int index, K;
    Vecteur<Tc> coefs;

    FiltreRIF(const Vecteur<Tc> &c)
    {
      this->coefs = c;
      index = 0;
      K = coefs.rows();
      wnd.setZero(K);
    }

    FRat<Tc> get_frat()
    {
      return FRat<Tc>::rif(coefs);
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      tsd_assert(K > 0);
      int n = x.rows();
      auto iptr = x.data();
      auto optr = y.data();
      if(iptr != optr)
      {
        y.resize(n);
        optr = y.data();
      }

      for(auto j = 0; j < n; j++)
      {
        auto cptr = coefs.data() + (K - 1);

        T sum = 0;

        // Remplace le plus ancien élément
        wnd(index) = *iptr++;
        index = (index + 1) % K;
        T *wptr = wnd.data() + index;

        int K1, K2;

        // Number of samples at end of delay line
        K1 = K - index;
        // Number at begin of delay line
        K2 = K - K1;

        for(auto i = 0; i < K1; i++)
          sum += *wptr++ * *cptr--;

        wptr = wnd.data(); // Restart at beginning of delay line
        for(auto i = 0; i < K2; i++)
          sum += *wptr++ * *cptr--;

        *optr++ = sum;
      }
    }
  };

  template<typename T>
  struct FiltreId: FiltreGen<T>
  {
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      y = x;
    }
  };

  template<typename T>
    sptr<FiltreGen<T>> filtre_id()
  {
    return make_shared<FiltreId<T>>();
  }


  template<typename T>
  struct Decimateur: FiltreGen<T, T>
  {
    int R = 0, cnt = 0;
    Decimateur(int R)
    {
      this->R = R;
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      int n = x.rows();

      y.resize((n + R - 1 - cnt) / R);

      int i, j = 0;
      for(i = cnt; i < n; i += R)
      {
        if(j >= y.rows())
        {
          echec("j={},y.rows()={},i={},n={},cnt={},R={}",j,y.rows(),i,n,cnt,R);
        }

        y(j++) = x(i);
      }

      i -= R;
      tsd_assert(j == y.rows());

      // Si i = n-R : cnt = 0
      // Si i = n-1 : cnt = R-1
      cnt = i - (n-R);

      // R - nombre d'échantillons non traités
      //cnt = (R - (n - 1 - i)) % R;
      //cnt = R - ((cnt + n) % R);
    }
  };

  template<typename T>
  sptr<FiltreGen<T>> decimateur(int R)
  {
    return std::make_shared<Decimateur<T>>(R);
  }

  template<typename Tc, typename T>
  sptr<FiltreGen<T>> filtre_rif(const Eigen::Ref<const Vecteur<Tc>> c)
  {
    return make_shared<FiltreRIF<T,Tc>>(c);
  }

  template<typename T, typename Tc>
  struct FiltreRII: FiltreGen<T>
  {
    Vecteur<Tc> numer, denom;
    Vecteur<T> wndx, wndy;
    int index = 0, index_y = 0;
    int Kx = 0, Ky = 0;

    FiltreRII(const FRat<Tc> &h)
    {
      // Attention, ici, il faut passer en puissance de z-1
      auto h2 = h.eval_inv_z();
      auto num = h2.numer.vers_coefs().coefs;
      auto den = h2.denom.vers_coefs().coefs;

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
      //cout << "RII : coefs numerateur = " << numer.transpose() << "\n";
      //cout << "RII : coefs denom = " << denom.transpose() << "\n";
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      int n = x.rows();
      auto iptr = x.data();
      auto optr = y.data();
      if(optr != iptr)
      {
        y.resize(n);
        optr = y.data();
      }

      // Forme directe I : numérateur puis dénominateur

      // (1) Filtre non récursif
      for(auto j = 0; j < n; j++)
      {
        auto nptr = numer.data() + (Kx - 1);
        T somme = 0;

        // Remplace le plus ancien élément
        wndx(index) = *iptr++;
        index = (index + 1) % Kx;
        T *wptr = wndx.data() + index;

        // Number of samples at end of delay line
        int KK1 = Kx - index;
        // Number at begin of delay line
        int KK2 = Kx - KK1;

        for(auto i = 0; i < KK1; i++)
          somme += *wptr++ * *nptr--;

        // Redémarre au début de la ligne à retard
        wptr = wndx.data();
        for(auto i = 0; i < KK2; i++)
          somme += *wptr++ * *nptr--;

        *optr++ = somme;
      }

      // (2) Filtre récursif
      optr = y.data();
      for(auto j = 0; j < n; j++)
      {
        auto dptr = denom.data() + 1; // Premier élément ignoré
        T somme = *optr;
        // Remplace le plus ancien élément
        // wndy(index) = *ioptr;
        // index = (index + 1) % Ky;
        T *wptr = wndy.data() + index_y;

        // Number of samples at end of delay line
        int KK1 = Ky - index_y;
        // Number at begin of delay line
        int KK2 = Ky - KK1;

        for(auto i = 0; i < KK1; i++)
          somme -= *wptr++ * *dptr++;

        // Redémarre au début de la ligne à retard
        wptr = wndy.data();
        for(auto i = 0; i < KK2; i++)
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
      return make_shared<FiltreRII<T,Tc>>(h);
    }



  template<typename T = float, int order = 2>
  struct AllpassUps: FiltreGen<T>
  {
    AllpassUps()
    {
      reset();
    }
    void reset()
    {
      xe_acc = xo_acc = xe_old = 0;
    }
    void process(T *y, const T *x, unsigned int n)
    {
      for(auto i = 0u; i < n; i++)
      {
        T xe = *x++;

        /* alpha =  1/2 */
        xe_acc = xe_old + ((xe - xe_acc) / 4.0f);
        /* alpha =  3/2 */
        xo_acc = xe_old + ((3.0f * (xe - xo_acc)) / 4.0f);

        *y++   = xe_acc;
        *y++ = xo_acc;

        xe_old = xe;
      }
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      tsd_assert(y.data() != x.data());
      y.resize(x.rows() * 2);
      process(y.data(), x.data(), x.rows());
    }
  private:
    T xe_old, xe_acc, xo_acc;
  };

  template<typename T>
  sptr<FiltreGen<T>> filtre_allpass_ups()
  {
    return make_shared<AllpassUps<T>>();
  }

  template<typename T>
  inline T decim_etape(const T &xe, const T &xo, T &xe_old, T &xo_old, T &xe_acc, T &xo_acc)
  {
    /* α =  1/2 */
    xe_acc = xe_old + ((xe - xe_acc) / 4.0f);
    /* α =  3/2 */
    xo_acc = xo_old + ((3.0f * (xo - xo_acc)) / 4.0f);

    xe_old = xe;
    xo_old = xo;

    return (xe_acc + xo_acc) / 2.0f;
  }

  /** @brief Décimateur réalisé à partir de deux passe-tout
   *         fonctionnant en parallèle. **/
  template<typename T = float, int order = 2>
  struct AllpassDecim : FiltreGen<T>
  {
    AllpassDecim()
    {
      reset();
    }
    void reset()
    {
      a_residu  = false;
      residu    = 0;
      xe_acc    = xo_acc = xe_old = xo_old = 0;
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      if(y.data() != x.data())
        y.resize((x.rows() + 1) / 2);
      process(y.data(), x.data(), x.rows());
      if(y.data() == x.data())
        y = y.head((x.rows() + 1) / 2);
    }
    int process(T *y, const T *x, unsigned int n)
    {
      int res = 0;
      uint32_t i;
      T xe, xo;

      if(n == 0)
        return 0;

      // Premier calcul avec le résidu
      if(a_residu)
      {
        a_residu = false;
        xe = residu;
        xo = *x++;

        *y++ = decim_etape<T>(xe, xo, xe_old, xo_old, xe_acc, xo_acc);
        res++;
        n--;
      }

      for(i = 0; i < n / 2; i++)
      {
        xe = *x++;
        xo = *x++;

        *y++ = decim_etape<T>(xe, xo, xe_old, xo_old, xe_acc, xo_acc);
      }
      res += n / 2;

      if(n & 1)
      {
        a_residu = true;
        residu = *x++;
        //assert(!(isnan(residu)));
      }

      return res;
    }
  private:
    T xe_old, xe_acc, xo_old, xo_acc;
    T residu; // Pour stocker le résidu si nombre d'échantillons impairs
    bool a_residu;
  };


  template<typename T>
  sptr<FiltreGen<T>> filtre_allpass_decim()
  {
    return make_shared<AllpassDecim<T>>();
  }

  struct SOISConfig
  {
    ArrayXf coefs;
    RIIStructure structure = FormeDirecte2;
  };



  /** @brief Second Order IIR Section */
  template<typename T = float, typename Tcoef = float, typename Ty = float>
  struct SOIS: Filtre<T, T, SOISConfig>
  {
    Ty x1, x2;
    Ty y0, y1, y2;
    Tcoef b0, b1, b2, a0, a1, a2;
    bool first_call;
    RIIStructure structure = RIIStructure::FormeDirecte2;

    SOIS()
    {
      a0 = a1 = a2 = b0 = b1 = b2 = 1;
      first_call = true;
    }

    /**       b0 + b1 z-1 + b2 z-2
     *  y/x = ---------------------
     *        a0 + a1 z-1 + a2 z-2
     */
    void configure(float b0, float b1, float b2, float a0, float a1, float a2)
    {
      tsd_assert_msg(a0 != 0, "Section SOIS : a0 doit être non nul.");

      b0 /= a0;
      b1 /= a0;
      b2 /= a0;
      a1 /= a0;
      a2 /= a0;
      a0 = 1.0;

      this->b0 = b0;
      this->b1 = b1;
      this->b2 = b2;
      this->a0 = 1.0;
      this->a1 = a1;
      this->a2 = a2;
      first_call = true;
      /*printf("Original coefs:      %f %f %f %f %f %f\n",
              b0, b1, b2, a0, a1, a2);
      printf("After approximation: %f %f %f %f %f %f\n",
              (float) this->b0, (float) this->b1, (float) this->b2,
              (float) this->a0, (float) this->a1, (float) this->a2);*/
    }
    int configure_impl(const SOISConfig &config)
    {
      auto coefs = config.coefs;
      structure = config.structure;
      tsd_assert_msg(coefs.rows() == 6,
                "Invalid number of coefficients for SOIC initialization (6 required).");
      configure(coefs(0), coefs(1), coefs(2), coefs(3), coefs(4), coefs(5));
      return 0;
    }
    void configure(const FRat<float> &h)
    {
      tsd_assert_msg(h.numer.coefs.rows() == 3,
                "Invalid number of coefficients for SOIC initialization (3 required).");
      tsd_assert_msg(h.denom.coefs.rows() == 3,
                "Invalid number of coefficients for SOIC initialization (3 required).");
      configure(h.numer.coefs(0), h.numer.coefs(1), h.numer.coefs(2),
            h.denom.coefs(0), h.denom.coefs(1), h.denom.coefs(2));
    }
    static void normalise_coefs(ArrayXf &coefs)
    {
      tsd_assert_msg(coefs(3) != 0, "Cannot normalize this section");
      coefs(0) /= coefs(3);
      coefs(1) /= coefs(3);
      coefs(2) /= coefs(3);
      coefs(4) /= coefs(3);
      coefs(5) /= coefs(3);
      coefs(3) = 1.0;
    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      auto n = x.rows();
      if(n == 0)
        return;

      if(x.data() != y.data())
        y.resize(n);

      Ty x0;
      auto iptr = x.data();
      auto optr = y.data();

      if(first_call)
      {
        y0 = y1 = y2 = x2 = x1 = x(0);
        first_call = false;
      }

      if(structure == RIIStructure::FormeDirecte2)
      {
        for(auto i = 0; i < n; i++)
        {
          x0 = *iptr++;

          /*y0  = b0 * x0; // ===> PB ICI: car b0*x0 est r�duit en entier avant la conversion en flottant
                             // C'est pourquoi le type de x0,x1 et x2 est Ty
          y0 += b1 * x1;
          y0 += b2 * x2;
          y0 -= a1 * y1;
          y0 -= a2 * y2;*/

          auto d2 = x0 - a1 * y1 - a2 * y0;
          y2 = b0 * d2 + b1 * y1 + b2 * y0;

          *optr++ = (T) y2;

          y0 = y1;
          y1 = d2;
        }
      }
      else if(structure == RIIStructure::FormeDirecte1)
      {
        for(auto i = 0; i < n; i++)
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
      else
      {
        echec("SOIS : structure non implémentée.");
      }
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

    int configure(Tcoef b0, Tcoef b1, Tcoef a1)
    {
      this->b0 = b0;
      this->b1 = b1;
      this->a1 = a1;
      this->x1 = 0;
      this->y1 = 0;
      return 0;
    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      auto n = x.rows();
      y.resize(n);
      for(auto i = 0; i < n; i++)
      {
        auto x0 = x(i);
        y1 = -a1 * y1 + b0 * x0 + b1 * x1;
        y(i) = y1;
        x1 = x0;
      }
    }
  };


  template<typename T, typename Tcoef, typename Ty>
  struct ChaineSOIS: FiltreGen<T>
  {
    float gain = 1.0f;
    vector<SOIS<T, Tcoef, Ty>> sections;
    RIIFoS<T, Tcoef> rii1;
    bool avec_rii1 = false;

    ChaineSOIS(const FRat<cfloat> &h, RIIStructure structure)
    {
      //infos("Construction d'une chaine de sections du second ordre.");
      auto z = h.numer.roots(), p = h.denom.roots();
      auto nz = z.rows(), np = p.rows();
      //infos("%d zéros, %d racines", nz, np);

      // pour l'instant
      tsd_assert_msg(nz == np, "SOIS : numerateur et dénominateur doivent avoir la même dimension");

      // Pool qui contient tous les index a priori, on vient puiser dedans
      set<int> pool;
      for(auto i = 0; i < nz; i++)
        pool.insert(i);

      // Groupe les racines et poles conjugés
      int i;
      for(i = 0; i + 1 < nz; i += 2)
      {
        auto k = *(pool.begin());
        pool.erase(pool.begin());

        // (z-(z(i)) * (z-z(i+1))
        // z² -(z(i)+z(i+1))*z + z(i)*z(i+1)

        // 1 -(z(i)+z(i+1))*z^-1 + z(i)*z(i+1)*z^-2

        SOIS<T, Tcoef, Ty> s;

        /**       b0 + b1 z-1 + b2 z-2
         *  y/x = ---------------------
         *        a0 + a1 z-1 + a2 z-2
         */
        // int setup(float b0, float b1, float b2, float a0, float a1, float a2)

        float berr = 1e9;

        cfloat sz, pz, sp, pp;
        auto bj = pool.begin();

        for(auto j = 0; j < nz; j++)
        {
          if(pool.count(j) == 0)
            continue;

          auto sz0 = -(z(j) + z(k)), pz0 = z(j)*z(k);
          auto sp0 = -(p(j) + p(k)), pp0 = p(j)*p(k);

          auto err = abs(sz0.imag()) + abs(pz0.imag()) + abs(sp0.imag()) + abs(pp0.imag());
          if(err < berr)
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
        //infos("Section %d : err = %e", i, berr);

        if(berr > 1e-5)
        {
          msg_erreur("SOIS factorisation : err = {}", berr);
        }

        ArrayXf coefs(6);
        coefs << 1.0f, sz.real(), pz.real(),
                 1.0f, sp.real(), pp.real();

        s.Configurable<SOISConfig>::configure({coefs, structure});

        sections.push_back(s);

      }

      for(; i < nz; i++)
      {
        //msg("Ordre impair : insertion d'une section du premier ordre");

        tsd_assert(pool.size() == 1);

        auto zer = z(*(pool.begin()));
        auto pol = p(*(pool.begin()));

        if(abs(pol.imag()) + abs(zer.imag()) > 1e-5f)
        {
          msg("pol.imag = {:e}, zer.imag = {:e}", pol.imag(), zer.imag());
          msg_erreur("Chaine RII: section du premier ordre restante : nombres imaginaires.");
        }

        /**       b0 + b1 z-1         z + b1/b0
         *  y/x = ----------- = b0 * ----------
         *        1 + a1 z-1            z + a1
         */

        float a1 = - pol.real();
        float b0 = h.numer.mlt.real() / h.denom.mlt.real();
        float b1 = - zer.real() * b0;

        rii1.configure(b0, b1, a1);
        avec_rii1 = true;
      }

      if(!avec_rii1)
      {
        gain = h.numer.mlt.real() / h.denom.mlt.real();
      }


    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      y = x;
      for(auto &s: sections)
        y = s.FiltreGen<T>::step(y);
      if(avec_rii1)
        y = rii1.FiltreGen<T>::step(y);
      else
        y = y * gain;
    }
  };


  /*template<typename Tc, typename T>
  sptr<Filtre<T>> filtre_SOIS(const Vecteur<Tc> &coefs)
  {
    auto res =  make_shared<SOIS<T,Tc,T>>();
    res->setup(coefs);
    return res;
  }*/

  template<typename T>
  sptr<FiltreGen<T>> filtre_sois(const FRat<cfloat> &h, RIIStructure structure)
  {
    return make_shared<ChaineSOIS<T,T,T>>(h, structure);
  }

  template<typename T>
  sptr<FiltreGen<T>> filtre_sois(const FRat<float> &h, RIIStructure structure)
  {
    auto numer = h.numer, denom = h.denom;




    FRat<cfloat> h2;
    h2.numer.mlt = h.numer.mlt;
    h2.numer.mode_racines = h.numer.mode_racines;
    h2.numer.coefs = h.numer.coefs;
    h2.denom.mlt = h.denom.mlt;
    h2.denom.mode_racines = h.denom.mode_racines;
    h2.denom.coefs = h.denom.coefs;

    if(!numer.mode_racines)
      h2.numer = Poly<cfloat>::from_roots(h.numer.roots()) * h.numer.coefs(h.numer.coefs.rows() - 1);
    if(!denom.mode_racines)
      h2.denom = Poly<cfloat>::from_roots(h.denom.roots()) * h.denom.coefs(h.denom.coefs.rows() - 1);


    return make_shared<ChaineSOIS<T,T,T>>(h2, structure);
  }

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  //// Second Order IIR Section
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////









float rii1_coef(float fcut)
{
  return 1.0 - exp(-fcut * 2 * π);
}

float rii1_tc_vers_coef(float τ)
{
  return rii1_coef(1.0 / (2 * π * τ));
}

float rii1_fcoupure(float γ)
{
  return -log(1.0 - γ) / (2 * π);
}

float rii1_coef_vers_tc(float γ)
{
  return 1.0 / (2 * π * rii1_fcoupure(γ));
}

static float bloqueur_dc_alpha(float fc)
{
  return (sqrt(3.0f) - 2 * sin(π * fc)) / (sin(π * fc) + sqrt(3.0f) * cos(π * fc));
}

FRat<float> design_bloqueur_dc(float fc)
{
  auto α = bloqueur_dc_alpha(fc);
  auto z = (FRat<float>::z()).inv();
  return (1 - z) / (1 - α * z);
}

FRat<float> design_rii1(float fc)
{
  // yn = a * xn + (1-a) yn-1
  // h = a / (1-(1-a)z^-1)
  //   = az / (z -(1-a)
  float a = rii1_coef(fc);
  auto z = FRat<float>::z();
  return (z * a) / (z  - (1-a));
}

template<typename T = float>
struct FiltreDC: FiltreGen<T>
{
  float α;
  T xp, yp;

  FiltreDC(float fc)
  {
    α = bloqueur_dc_alpha(fc);
    xp = yp = (T) 0;
  }
  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    auto n = x.rows();

    if((void *) &y != (void *) &x)
      y.resize(n);

    for(auto i = 0; i < n; i++)
    {
      y(i) = x(i) - xp + α * yp;
      xp = x(i);
      yp = y(i);
    }
  }
};


/*************************************************************
 *************************************************************
    MOBILE AVERAGE LOWPASS FILTER (TIME DOMAIN DOOR FILTER)
 *************************************************************
 *************************************************************/
template<typename T = float, typename Tacc = float, typename Tinv = float>
struct MoyenneGlissante: FiltreGen<T>
{
  int index = 0, K = 0;
  Tacc accu;
  Vecteur<T> fenetre;
  // Inverse of length
  Tinv K_inv;

  MoyenneGlissante(int n)
  {
    accu  = 0;
    index = 0;
    K     = n;
    fenetre.setZero(n);
    K_inv = (Tinv) (1.0 / ((double) n));
  }
  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    int n = x.rows();

    tsd_assert(K > 0);

    if((void *) &y != (void *) &x)
      y.resize(n);

    for(auto i = 0; i < n; i++)
    {
      accu += x(i);
      accu -= fenetre(index);
      fenetre(index) = x(i);
      index = (index + 1) % K;  //(first + 1 < len) ? first + 1 : 0;
      y(i) = ((T) accu) * K_inv;
    }
  }
};


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//// Exponential first order LP filter
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

/** @brief First order exponential high-pass filter */
template<typename T = float, typename Tacc = float, typename Tgamma = float>
struct ExpHpFilter
{
  Tacc   acc;
  // Filter coefficient
  Tgamma γ;
  bool   init;

  ExpHpFilter()
  {
    γ = 0;
    reset();
  }

  int setup(float γ)
  {
    this->γ = γ;

    msg("lexp setup(γ = {}): approx = {}", γ, (float) this->γ);

    reset();

    return 0;
  }

  void reset()
  {
    init = false;
    acc  = 0;
  }

  /** @param[in] in   Input matrix.
   *  @param[out] out Output matrix (can be the same object as input if in place computation is needed). */
  void process(Vecteur<T> &out, const Vecteur<T> &in)
  {
    uint32_t i, nsples = in.rows();

    if(nsples == 0)
      return;

    out.resize(nsples);

    if(!init)
      acc = in(0); // Will automatically shift input data of n bits, if acc
                   // is fixed point.

    for(i = 0; i < nsples; i++)
    {
      Tacc x = in(i);
      acc += γ * (x - acc);
      out(i) = (T) (x - acc);
    }

    init = true;
  }

};

/** @brief First order exponential low-pass filter
 *  Canonical form:
 *
 *  y_n+1 = (1-g) y_n + g x_n e.g. y_n+1 = y_n + g(x_n - y_n) */
template<typename T = float, typename Tacc = float, typename Tgamma = float>
struct ExpLpFilter : FiltreGen<T>
{
  Tacc   acc   = 0;
  Tgamma γ = 0;
  bool   init  = false;

  ExpLpFilter(float γ)
  {
    this->γ = γ;
    init = false;
    acc  = 0;
  }

  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    auto n = x.rows();

    if(n == 0)
      return;

    y.resize(n);

    if(!init)
      acc = x(0); // Will automatically shift input data of n bits, if acc
                   // is fixed point.

    // acc_i = xpi + (1-g)*acc_{i-1}


    for(auto i = 0; i < n; i++)
    {
      Tacc xi = x(i);
      acc += γ * (xi - acc);
      y(i) = (T) acc;
    }

    init = true;
  }
};


template<typename T> sptr<FiltreGen<T>> filtre_rii1(float γ)
{
  return make_shared<ExpLpFilter<T,T,T>>(γ);
}

template<typename T>
  sptr<FiltreGen<T>> filtre_dc(float fc)
{
  return make_shared<FiltreDC<T>>(fc);
}

template<typename T, typename Tacc>
  sptr<FiltreGen<T>> filtre_mg(int n)
{
  return make_shared<MoyenneGlissante<T,Tacc,T>>(n);
}


template
  Tableau<float> forme_polyphase<float>(const Eigen::Ref<const Vecteur<float>> x, unsigned int M);

template
  Tableau<cfloat> forme_polyphase<cfloat>(const Eigen::Ref<const Vecteur<cfloat>> x, unsigned int M);

template
  Vecteur<float> iforme_polyphase<float>(const Eigen::Ref<const Tableau<float>> x);

template
  Vecteur<cfloat> iforme_polyphase<cfloat>(const Eigen::Ref<const Tableau<cfloat>> x);


// Instanciations explicites.
template sptr<FiltreGen<float>> filtre_rii1(float γ);
template sptr<FiltreGen<cfloat>> filtre_rii1(float γ);

template sptr<FiltreGen<float>> filtre_dc(float fc);
template sptr<FiltreGen<cfloat>> filtre_dc(float fc);

template sptr<FiltreGen<float>> filtre_rii(const FRat<float> &h);
template sptr<FiltreGen<cfloat>> filtre_rii(const FRat<cfloat> &h);
template sptr<FiltreGen<float>> ligne_a_retard<float>(int retard);
template sptr<FiltreGen<float>> filtre_mg<float,double>(int n);
template sptr<FiltreGen<cfloat>> filtre_mg<cfloat,cdouble>(int n);


template sptr<FiltreGen<float>> filtre_sois<float>(const FRat<cfloat> &coefs, RIIStructure structure);
template sptr<FiltreGen<float>> filtre_sois<float>(const FRat<float> &coefs, RIIStructure structure);
template sptr<FiltreGen<cfloat>> filtre_sois<cfloat>(const FRat<cfloat> &coefs, RIIStructure structure);
template sptr<FiltreGen<cfloat>> filtre_sois<cfloat>(const FRat<float> &coefs, RIIStructure structure);

template tuple<ArrayXf, ArrayXf> frgroup<float>(const FRat<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frgroup<float>(const Vecteur<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frmag<float>(const FRat<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frmag<float>(const Vecteur<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frphase<float>(const FRat<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frphase<float>(const Vecteur<float> &h, int npts);
template tuple<ArrayXf, ArrayXf> frgroup<cfloat>(const FRat<cfloat> &h, int npts);
template tuple<ArrayXf, ArrayXf> frgroup<cfloat>(const Vecteur<cfloat> &h, int npts);
template tuple<ArrayXf, ArrayXf> frmag<cfloat>(const FRat<cfloat> &h, int npts);
template tuple<ArrayXf, ArrayXf> frmag<cfloat>(const Vecteur<cfloat> &h, int npts);
template tuple<ArrayXf, ArrayXf> frphase<cfloat>(const FRat<cfloat> &h, int npts);
template tuple<ArrayXf, ArrayXf> frphase<cfloat>(const Vecteur<cfloat> &h, int npts);


namespace hidden {
auto filtre_id1  = filtre_id<float>;
auto filtre_id2  = filtre_id<cfloat>;

auto decimateur1  = decimateur<float>;
auto decimateur2  = decimateur<cfloat>;

auto filtre_rif1 = filtre_rif<float, float>;
auto filtre_rif2 = filtre_rif<float, cfloat>;
//auto filtre_rif3 = filtre_rif<cfloat, float>;
auto filtre_rif4 = filtre_rif<cfloat, cfloat>;

auto filtre_ad1 = filtre_allpass_decim<float>;
auto filtre_ad2 = filtre_allpass_decim<cfloat>;
auto filtre_up1 = filtre_allpass_ups<float>;
auto filtre_up2 = filtre_allpass_ups<cfloat>;
/*auto frmag_1 = frmag<float>;
auto frphase_1 = frphase<float>;
auto frmag_2 = frmag<cfloat>;
auto frphase_2 = frphase<cfloat>;*/


}

}
