  #include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include "tsd/moniteur-cpu.hpp"


#include <bit>
#include <limits>

const auto FOURIER_MODE_SAFE = false;

template<typename derived>
void checkNaN(const Eigen::ArrayBase<derived> &x, const std::string &s)
{
  if(x.hasNaN())
    echec("{} : has NaN.", s);
}

namespace tsd::fourier
{



  template<typename T>
  auto fft_twiddles_ref(int N, int signe = -1)
  {
    Vecteur<std::complex<T>> res(N);
    for(auto k = 0; k < N; k++)
      res(k) = std::polar<T>(1.0f, signe*2.0*π*k/N);
    return res;
  }

  template<typename T>
    auto fft_twiddles_fast(int N, int signe = -1)
    {
      Vecteur<std::complex<T>> x(N);
      cdouble r = 1, w0 = std::polar<double>(1.0, (signe*2*π)/N);

      // TODO : faire les calculs sur un seul quadrant
      for(auto i = 0; i < N; i++)
      {
        x(i) = r;
        r *= w0;
      }

      return x;
    }

  template<typename T>
    auto fft_twiddles(int N, int signe = -1)
  {
    return fft_twiddles_fast<T>(N, signe);
  }

  template<typename T>
  auto ifft_twiddles(int N)
  {
    return fft_twiddles<T>(N, 1);
  }

  void check_csym(const ArrayXcf &X)
  {
    int n = X.rows();

    ArrayXf verr;
    if((n % 2) == 0)
      verr = (X.segment(1, n/2-1) - X.tail(n/2-1).reverse().conjugate()).abs();
    else
      verr = (X.segment(1, n/2) - X.tail(n/2).reverse().conjugate()).abs();

    int idx;
    float err = verr.maxCoeff(&idx);
    idx++;

    err = std::max(err, std::abs(imag(X(0))));
    if((n % 2) == 0)
      err = std::max(err, std::abs(imag(X(n/2))));

    fmt::print("Csym erreur : {}, imag(0) = {}, imag(n/2) = {}, deseq({}) = {}\n",
        err, imag(X(0)), imag(X(n/2)), idx, err);
  }




  /** @brief Calcul FFT Cooley–Tukey (N doit être une puissance de 2) */
  template<typename T, typename Ttwiddle, bool avant>
  int fft_radix2(Vecteur<std::complex<T>> &X,
                 const Vecteur<std::complex<T>> &x,
                 Vecteur<std::complex<T>> &scratch,
                 const Vecteur<std::complex<Ttwiddle>> &twiddles)
  {
    int32_t N = x.rows();
    auto iteration_paire = ((N & 0x55555555) != 0);

    scratch.resize(N);
    X.resize(N);
    dsp_check_safe(twiddles.rows() == N, "twiddles matrix not init.");

    if(N == 1)
    {
      X(0) = x(0);
      return 0;
    }

    auto E = x.data();
    std::complex<T> *Xp, *Xp2, *Xstart;

    // Commence par calculer les FFT de plus bas niveau
    // puis remonte petit à petit
    for (int n = 1; n < N; n *= 2)
    {
      // Premiére itération : E = x, => Xp, Xs = scratch
      // Deuxième itération : E = scratch, Xs = X
      // ....
      Xstart = iteration_paire ? scratch.data() : X.data();
      int pas = N / (2 * n);
      Xp  = Xstart;
      Xp2 = Xstart + N / 2;
      //auto tptr = twiddles.data();

      /*Eigen::Map<Eigen::MatrixXcf> Mat(E, pas, 2 * n);

      Eigen::VectorXcf e1(2*n), e2(2*n);
      for(auto k = 0; k < n; k++)
      {
        e1(2*k)   = 0.5f;
        e1(2*k+1) = twiddles(k * pas);
        e2(2*k)   = 0.5f;
        e2(2*k+1) = -twiddles(k * pas);
      }*/

      //Eigen::Map<ArrayXcf> map_Xp1(Xp, pas);
      //Eigen::Map<Eigen::VectorXcf> Mat(E, pas, 2 * n);
      //Eigen::Map<ArrayXcf> map_big_Xp(Xp, pas * n), map_big_Xp2(Xp2, pas * n);
      //Eigen::Map<const ArrayXcf> map_big_E(E, 2 * pas * n);
      //map_big_Xp =


      for(int k = 0; k < n; k++)
      {
        auto twid = twiddles(k * pas);

        if constexpr(!avant)
            twid = conj(twid);

#       if 0
        Eigen::Map<ArrayXcf> map_Xp(Xp, pas), map_Xp2(Xp2, pas);
        Eigen::Map<const ArrayXcf> map_E0(E, pas), map_E1(E + pas, pas);

        auto p  = twid * map_E1;
        map_Xp  = map_E0/2 + p;
        map_Xp2 = map_E0/2 - p;
        E += pas;
        Xp += pas;
        Xp2 += pas;
#       else
        const auto tr = twid.real(), ti = twid.imag();
        for(int m = 0; m < pas; m++)
        {
          const auto g = *(E + pas), e = *E;
          const auto gx = g.real(), gy = g.imag();
          const std::complex<T> p(tr * gx - ti * gy, tr * gy + ti * gx);

          //*Xp++  = e/2.0f + p;
          //*Xp2++ = e/2.0f - p;
          *Xp++  = e + p;
          *Xp2++ = e - p;
          E++;
        }
#       endif
        E += pas;
      }
      E = Xstart;
      iteration_paire = !iteration_paire;
    }

    // On a divisé par 2 à chaque itération au lieu de diviser par srqt(2)
    //X *= std::sqrt((float) N);
    X /= std::sqrt((float) N);
    return 0;
  }



  // A SUPPRIMER
  /** @brief Compute the FFT of real samples:
   *  From n real samples, produce n/2 complex value (only positive frequencies) **/
  /** Algorithm: numerical recipes, p512 */
  template<typename T>
  Vecteur<std::complex<T>> rfft_compute(const Vecteur<T> &x)
  {
    int n = x.rows();
    Vecteur<std::complex<T>> x2(n/2), Xt(n/2), X(n/2);

    for(auto i = 0; i < n / 2; i++)
    {
      x2(i).real(x(2*i));
      x2(i).imag(x(2*i+1));
    }

    // Compute n/2 FFT points in X
    Xt = fft(x2);

    const std::complex<T> j2(0,0.25);
    const std::complex<T> r2(0.25,0);

    Vecteur<std::complex<T>> tw(n/2+1);

    for(auto k = 0; k <= n/2; k++)
    {
      tw(k).real(cos(-(2*π*k)/n));
      tw(k).imag(sin(-(2*π*k)/n));
    }

    for(auto i = 0; i <= n / 2; i++)
    {
      std::complex<T> X1, X2, X3;

      if(i == n/2)
        X1 = Xt(0);
      else
        X1 = Xt(i);

      if(i > 0)
        X2 = Xt(n/2-i);
      else
        X2 = Xt(0);

      X3 = r2 * (X1 + conj(X2)) - j2 *  (X1 - conj(X2)) * tw(i);

      if(i < n / 2)
        X(i) = X3;
      else
        X(0).imag(X3.real()); // Note: actually X(0).imag() = 0, but we store here X(n/2)
    }

    return X;
  }


  /** @brief Restore real samples from complex positive frequencies
   *  From n/2 complex values (positive frequencies), produces n real samples **/
  /** Algorithm: numerical recipes, p512 */
  template<typename basetype>
  Vecteur<basetype> irfft_compute(const Vecteur<std::complex<basetype>> &X) // Size n/2
  {
    int n = 2 * X.rows();

    //dsp_check(x.rows() == (int) n, "Output IFFT size must be n (= 2 size of complex FFT)");


    Vecteur<basetype> x(n);
    Vecteur<std::complex<basetype>> X2(n/2);
    std::complex<basetype> J(0,1.0);

    for(auto i = 0; i < n / 2; i++)
    {
      std::complex<basetype> Xi, Xp;

      if(i == 0)
        Xp = X(0).imag();
      else
        Xp = X(n/2-i);

      if(i == 0)
        Xi = X(0).real();
      else
        Xi = X(i);

      std::complex<basetype> Xe = (Xi + conj(Xp)),
                             Xo = (Xi - conj(Xp))
                                * std::polar<basetype>((double) 1.0, (double) (2*π*i)/n);


      X2(i) = Xe + J * Xo;
    }

    Vecteur<std::complex<basetype>> x2(n/2);

    // Compute n/2 FFT points in X
    //fft_compute<basetype, Ttwiddle, false>(x2, X2, scratch, itwiddles);
    x2 = ifft(X2);

    // x2 is the same as x, but viewed as a complex vector of n/2 samples
    for(auto i = 0; i < n / 2; i++)
    {
      x(2*i)   = x2(i).real();
      x(2*i+1) = x2(i).imag();
    }

    return x;
  }




  /*
   *
   *     //if(!std::ispow2(N))
    if((N & (N - 1)) != 0)
    {
      infos("Avertissement : FFT sur N != 2^k");
      int N2 = prochaine_puissance_de_2(N);

    }
   *
   */


  // twidlles  : n2
  // chirp     : 2*n-1
  ArrayXcf fft_czt_impl(const ArrayXcf &x, int n2, const ArrayXcf &twiddles, const ArrayXcf &chirp)
  {
    int n = x.rows();

    // chirp = exp(-2pi * i * [-(n-1)/n, -(n-2)/n, ... 0, ..., (n-1)/n])
    ArrayXcf xp = ArrayXcf::Zero(n2);
    xp.head(n) = x * chirp.tail(n);

    ArrayXcf icp = ArrayXcf::Zero(n2);
    icp.head(2*n-1) = chirp.conjugate();

    ArrayXcf scratch(n2);
    ArrayXcf y, y2;

    //y = ifft(fft(xp) * fft(icp))
    ArrayXcf Xp, Xc;
    fft_radix2<float, float, true>(Xp, xp,      scratch, twiddles);
    fft_radix2<float, float, true>(Xc, icp,     scratch, twiddles);
    fft_radix2<float, float, false>(y2, Xp * Xc, scratch, twiddles);
    return y2.segment(n-1, n) * chirp.tail(n) * (sqrt(n2) / sqrt(n));
  }

  /** Conversion du résultat d'une fft vers une ifft
   *  E.g. Y(n) = X(N-n) */
  ArrayXcf fftversifft(const ArrayXcf &X)
  {
    int n = X.rows();
    ArrayXcf Y(n);
    Y(0) = X(0);

    if((n & 1) == 0)
    {
      Y(n/2) = X(n/2);
      Y.segment(1, n/2-1) = X.tail(n/2-1).reverse();
      Y.tail(n/2-1) = X.segment(1, n/2-1).reverse();
    }
    else
    {
      Y.segment(1, n/2) = X.tail(n/2).reverse();
      Y.tail(n/2) = X.segment(1, n/2).reverse();
    }

    //fmt::print("fftversifft: X = {}, Y = {}\n", X, Y);

    return Y;
  }

  struct FFTPlan;

  template<typename T>
  struct RFFTPlan: FiltreGen<T, std::complex<T>>
  {
    int n;
    sptr<FFTPlan> cplan;

    ArrayXcf twiddles;

    RFFTPlan(int n)
    {
      configure(n);
    }
    int configure(int n)
    {
      this->n = n;

      if(n > 0)
      {
        if((n & 1) == 0)
        {
          cplan     = fftplan_création(n/2);
          twiddles  = fft_twiddles<T>(n);
        }
        else
        {
          // Si n est impair, on fait juste bêtement une FFT.
          cplan = fftplan_création(n);
        }
      }

      return 0;
    }
    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<std::complex<T>> &y)
    {
      if(x.rows() != n)
        configure(x.rows());
      if((n & 1) == 0)
      {
        y.resize(n);
        Vecteur<std::complex<T>> x2(n/2), Xt(n/2);

        for(auto i = 0; i < n / 2; i++)
        {
          x2(i).real(x(2*i));
          x2(i).imag(x(2*i+1));
        }

        // Compute n/2 FFT points in X
        cplan->step(x2, Xt, true);

        const std::complex<T> j2(0,0.5 / sqrt(2.0f));
        const std::complex<T> r2(0.5 / sqrt(2.0f),0);

        for(auto i = 0; i <= n / 2; i++)
        {
          std::complex<T> X1, X2;

          if(i == n/2)
            X1 = Xt(0);
          else
            X1 = Xt(i);

          if(i > 0)
            X2 = Xt(n/2-i);
          else
            X2 = Xt(0);

          y(i) = r2 * (X1 + conj(X2)) - j2 *  (X1 - conj(X2)) * twiddles(i);
        }
        force_csym(y);
      }
      else
      {
        ArrayXcf x2 = x;
        cplan->step(x2, y, true);
      }
    }
  };




  struct FFTPlanDefaut: FiltreGen<cfloat>, FFTPlan
  {
    bool normalize = true;
    bool avant = true;
    int n = 0, n2 = 0;
    ArrayXcf scratch, twiddles, chirp;

    sptr<FFTPlanDefaut> sousplan;

    FFTPlanDefaut(int n = -1, bool avant = true, bool normalize = true)
    {
      configure(n, avant, normalize);
    }

    int configure(int n, bool avant, bool normalize)
    {
      this->n = n;
      this->avant = avant;

      if(n == -1)
        return 0;

      n2 = n;


      if((n & (n - 1)) != 0)
      {
        // n est pair
        if((n & 1) == 0)
        {
          sousplan = std::make_shared<FFTPlanDefaut>(n / 2, avant, normalize);
        }
        else
        {
          //infos("Avertissement : FFT sur n != 2^k (n = %d)", n);
          n2 = prochaine_puissance_de_2(2*n-1);

          ArrayXf t = linspace(-(n-1), n-1, 2*n-1).square() / 2;
          t *= -2*π/n;
          // t = -2pi * [-(n-1)/n, -(n-2)/n, ... 0, ..., (n-1)/n]
          chirp = polar(t);
        }
      }
      //msg("plan fft : n = {}, n2 = {}", n, n2);

      scratch.resize(n2);
      //if(avant || (n2 != n))
        twiddles = fft_twiddles<float>(n2);
      //else
        //twiddles = ifft_twiddles<float>(n2);
      return 0;
    }

    void step(const Eigen::Ref<const ArrayXcf> x, ArrayXcf &y)
    {
      step(x, y, this->avant);
    }

    void step(const Eigen::Ref<const ArrayXcf> x, ArrayXcf &y, bool avant)
    {
      tsd_assert(x.rows() > 0);

      if((int) n != x.rows())
        configure(x.rows(), this->avant, this->normalize);

      if(n2 != n)
      {
        y = fft_czt_impl(x, n2, twiddles, chirp);
        if(!avant)
        {
          y = fftversifft(y);
        }
      }
      else
      {
        // Puissance de 2
        if((n & (n - 1)) == 0)
        {
          y.resize(n);
          if(avant)
            fft_radix2<float, float, true>(y, x, scratch, twiddles);
          else
            fft_radix2<float, float, false>(y, x, scratch, twiddles);
        }
        else
        {
          ArrayXcf xe(n/2), xo(n/2);
          for(auto i = 0; i < n/2; i++)
          {
            xe(i) = x(2*i);
            xo(i) = x(2*i+1);
          }

          ArrayXcf E, O;
          sousplan->step(xe, E, avant);
          sousplan->step(xo, O, avant);
          y.resize(n);

          if(avant)
          {
            y.head(n/2) = E + twiddles.head(n/2) * O;
            y.tail(n/2) = E + twiddles.tail(n/2) * O;
          }
          else
          {
            y.head(n/2) = E + twiddles.head(n/2).conjugate() * O;
            y.tail(n/2) = E + twiddles.tail(n/2).conjugate() * O;
          }
          if(normalize)
            y *= 1.0f / std::sqrt(2.0f);
        }
      }
    }

  };

  std::function<sptr<FFTPlan>()> fftplan_defaut = []()
  {
      return std::make_shared<FFTPlanDefaut>();
  };


  sptr<FFTPlan> fftplan_création(int n, bool avant, bool normalize)
  {
    auto res = fftplan_defaut();
    res->configure(n, avant, normalize);
    return res;
  }

  sptr<FiltreGen<float, cfloat>> rfftplan_création(int n)
  {
    return std::make_shared<RFFTPlan<float>>(n);
  }


  // Produit de corrélation effectué dans le domaine fréquentiel
  template<typename T>
  Vecteur<T> correlation_freq(const Vecteur<T> &X0, const Vecteur<T> &X1)
  {
    int n = X0.rows();
    tsd_assert(n == X1.rows());
    Vecteur<T> Y(n);

    Y(0) = X0(0) * conj(X1(0));
    Y.tail(n-1) = X0.tail(n-1).reverse() * X1.tail(n-1).reverse().conjugate() ;

    Y *= std::sqrt(n);
    return Y;
  }


  struct FFTCorrelateurBloc
  {
    FFTPlanDefaut plan;

    /** @note Should pad x0 and x1 with K zeros if K first lags are examined
     *  (e.g. n zeroes if all lags are used).
     *
     *  Result: y(0): lag = 0
     *          y(1): lag = 1
     *          ....
     *          y(n-2): lag = -2
     *          y(n-1): lag = -1
     */
    ArrayXcf step(const ArrayXcf &x0, const ArrayXcf &x1)
    {
      ArrayXcf x1p = x1.rows() == 0 ? x0 : x1;

      tsd_assert_msg(x0.rows() == x1p.rows(), "the two input vectors should have the dimension {} != {}.", x0.rows(), x1.rows());

      //msg("fft correlateur : ffts...");
      // (1) Compute FFT on x0 and x1
      ArrayXcf X0, X1, X2, y;

      plan.step(x0, X0, true);
      plan.step(x1p, X1, true);

      //msg("fft correlateur : cfreq...");
      X2 = correlation_freq(X0, X1);

      //msg("fft correlateur : ifft...");
      plan.step(X2, y, false);

      // y(0): lag = 0
      // y(1): lag = 1
      // ...
      // y(n-1): lag = n-1 == -1 from the periodicity hypothesis.

      return y.conjugate();
    }
  };


  std::tuple<ArrayXf, ArrayXcf> ccorr(const ArrayXcf &x0, const ArrayXcf &x1)
  {
    FFTCorrelateurBloc fc;
    int m = x0.rows();
    ArrayXf lags = linspace(0, m-1, m);
    return {lags, fc.step(x0, x1) / m};
  }

  std::tuple<ArrayXf, ArrayXcf> xcorrb(const ArrayXcf &x, const ArrayXcf &y, int m)
  {
    FFTCorrelateurBloc fc;

    int n = x.rows();

    if(m < 0)
      m = n;

    ArrayXcf yp = y.rows() == 0 ? x : y;

    ArrayXcf x2 = ArrayXcf::Zero(m + n + m);
    ArrayXcf y2 = ArrayXcf::Zero(m + n + m);

    x2.segment(m, n) = x;
    y2.segment(m, n) = yp;

    ArrayXcf r = fc.step(x2, y2);

    ArrayXcf res(2*m-1);
    res.tail(m) = r.head(m) / n;
    res.head(m-1) = r.tail(m-1) / n;
    ArrayXf lags = linspace(-(m-1), m-1, 2*m-1);

    return {lags, res};
  }

  std::tuple<ArrayXf, ArrayXcf> xcorr(const ArrayXcf &x, const ArrayXcf &y, int m)
  {
    int n = x.rows();

    if(m < 0)
      m = n;

    // Corrélation biaisée
    auto [lags, zb] = xcorrb(x, y, m);


    zb.head(m-1) /= linspace(n-(m-1), n-1, m-1) / n;
    zb.tail(m-1) /= linspace(n-1, n-(m-1), m-1) / n;

    return {lags, zb};
  }


  ArrayXf coherence(const ArrayXcf &x, const ArrayXcf &y)
  {
    ArrayXcf X = fft(x), Y = fft(y);
    return (X * Y.conjugate()).abs() / (X.abs() * Y.abs());
  }


  template<typename T>
  Vecteur<std::complex<T>> delais_fractionnaire_c(const Vecteur<std::complex<T>> &x, float τ)
  {
    int n = 2 * x.rows();

    // (1) Complète par des zéros à gauche et à droite
    Vecteur<std::complex<T>> x2 = Vecteur<std::complex<T>>::Zero(n);
    x2.segment(n/4, n/2) = x;

    // Modulation dans le domaine fréquentiel
    Vecteur<std::complex<T>> X = fft(x2);

    //ArrayXcf rot = sigexp(-τ / n, n);
    //force_csym(rot);

    ArrayXcf rot(n);
    for(auto i = 0; i < n; i++)
    {
      rot(i) = std::polar(1.0f, -2*π_f*i*τ/n + π_f*τ);
    }

    rot = fftshift(rot);

    X *= rot;

    // Retire le padding
    return ifft(X).segment(n/4, n/2);
  }



  // PB : ne marche pas si T = complexe
  template<typename T = float>
  Vecteur<T> delais_fractionnaire_r(const Vecteur<T> &x, float τ)
  {
    uint32_t i, n = 2 * x.rows();

    // (1) Padd x with zeroes at begin and end
    Vecteur<T> x2 = Vecteur<T>::Zero(n);
    x2.segment(n/4, n/2) = x;

    // A delay in time domain is a modulation in frequency domain
    Vecteur<std::complex<T>> X = rfft_compute(x2);

    std::complex<T> rot(1.0, 0), dphi;

    dphi = std::polar<T>((T) 1.0, (T) - τ * 2 * π / n);

    // TODO
    // Attention X(0).real -> freq(0)
    //           X(0).imag -> freq(n/2)
    for(i = 0; i < n / 2; i++)
    {
      if(i == 0)
        X(i).real((X(i).real() * rot).real());
      else
        X(i) *= rot;

      rot  *= dphi;
    }

    X(0).imag((X(0).imag() * rot).real());

    Vecteur<T> y2 = irfft_compute(X);

    // Now depadd
    return y2.segment(n/4, n/2);
  }

  template<typename T = float>
  Vecteur<T> delais_entier(const Vecteur<T> &x, int τ)
  {
    if(τ == 0)
      return x;

    int n = x.rows();
    Vecteur<T> y = Vecteur<T>::Zero(n);

    // Delais positif
    if(τ > 0)
      y.tail(n - τ) = x.head(n - τ);
    // Delais négatif
    else
      y.head(n + τ) = x.tail(n + τ);

    return y;
  }


  template<typename T>
    Vecteur<T> délais(const Vecteur<T> &x, float delais)
  {
    if(std::floor(delais) != delais)
    {
      if constexpr (std::is_same<T, cfloat>::value)
        return delais_fractionnaire_c(x, delais);
      else
        return delais_fractionnaire_r(x, delais);
    }
    else
    {
      return delais_entier<T>(x, (int) delais);
    }
  }




//[C,Nf,Nz] = complexite(fe, fs, N, Ne)

// Complexité, par échantillon d'entrée, en FLOPS
// M  : taille de filtre ou motif
// Ne : taille de bloc d'entrée
void ola_complexité(int M, int Ne, float &C, int &Nf, int &Nz)
{
  Nf = prochaine_puissance_de_2(Ne + M - 1);
  Nz = Nf - Ne;
  C = (1.0f / Ne) * 2 * 5 * Nf * log(1.0f*Nf) / log(2.0f);
}

void ola_complexité_optimise(int M, float &C_, int &Nf_, int &Nz_, int &Ne_)
{
  // Nf doit être au minimum égal à M
  int kmin = ceil(log(M) / log(2));

  for(int k = kmin; k < kmin + 20; k++)
  {
    int Nf, Nz, Ne = (1 << k) - (M - 1);
    float C;
    ola_complexité(M, Ne, C, Nf, Nz);
    if((k == kmin) || (C < C_))
    {
      Nf_   = Nf;
      Nz_   = Nf - Ne;
      Ne_   = Ne;
      C_    = C;
    }
  }
}







template<typename T = cfloat>
struct OLA: Filtre<T, T, FiltreFFTConfig>
{
  ArrayXcf padded, X, last, x2;
  ArrayXcf svg;
  // TODO !
  FFTPlanDefaut plan;
  ArrayXf fenêtre;

  /** Nb ech. TFD */
  int N = 0;

  /** Nb zéros insérés */
  int N_zeros = 0;

  /** Nb ech. d'entrée / sortie */
  int Ne = 0;

  /** Nombre d'échantillons en cours */
  int64_t cnt_ech = 0;

  sptr<SinkGen<cfloat>> tampon;
  bool tampon_vide = true;

  std::vector<ArrayXcf> bo;

  int configure_impl(const FiltreFFTConfig &config)
  {
    Ne      = config.dim_blocs_temporel;

    if(!config.traitement_freq)
    {
      msg_erreur("configuration OLA : traitement fréquentiel non précisé.");
      return -1;
    }

    if(Ne <= 0)
      Ne = 512;

    // Attention, post-affectation
    //Configurable<FiltreFFTConfig>::config.dim_blocs_temporel = Ne;

    N       = prochaine_puissance_de_2(Ne + config.nb_zeros_min);
    N_zeros = N - Ne;

    // Position initiale
    cnt_ech = - ((int) Ne) / 2;

    msg("OLA init: Ne = {} (nb échan entrée), N = {} (dim TFD), Nz = {} (nb zéros).", Ne, N, N_zeros);

    padded.setZero(N);
    last.setZero(Ne);
    svg.setZero(Ne);

    /*if((config.facteur_recouvrement != 1) && (config.facteur_recouvrement != 2))
    {
      erreur("OLA::configure : facteur de recouvrement : doit etre 2.");
      return -1;
    }*/

    if(config.avec_fenetrage)
    {
      fenêtre = tsd::filtrage::fenetre("hn", Ne, false);
      tsd_assert(fenêtre.rows() == Ne);
    }

    //plan  = creation_fft_plan(N, true);
    //iplan = creation_fft_plan(N, false);


    tampon_vide = true;
    tampon = tampon_création<cfloat>(Ne,
        [&](const Vecteur<cfloat> &x)
        {
          ArrayXcf y;
          step_interne(x, y);

          //msg("callback tampon : {} --> {}", x.rows(), y.rows());

          bo.push_back(y);
        });

    return 0;
  }

  void step(IArrayXcf x, ArrayXcf &y)
  {
    if(tampon_vide && (x.rows() == Ne))
    {
      step_interne(x, y);
      return;
    }
    tampon_vide = false;
    tampon->step(x);
    auto n = 0;
    for(auto &b: bo)
      n += b.rows();
    y.resize(n);
    auto i = 0;
    for(auto &b: bo)
    {
      y.segment(i, b.rows()) = b;
      i += b.rows();
    }
    //msg("OLA step : {} échan en entrée, Ne = {}, nb échans sortie : {}.", x.rows(), Ne, y.rows());
    bo.clear();
  }



  void step_interne(IArrayXcf x, ArrayXcf &y)
  {
    auto &config = Configurable<FiltreFFTConfig>::config;

    tsd_assert(x.rows() == Ne);

    if(!config.avec_fenetrage)
    {
      // OLA simple (pas de recouvrement sur les tampons d'entrée,
      // recouvrement sur les tampons de sorties)

      //msg("OLA/STI : Ne = {}, padded.rows()={}", Ne, padded.rows());

      // Ne = 16k, Padded = 32k
      padded.tail(Ne) = x;


      // Conso CPU des plans : ~ 20 % chacun

      plan.step(padded, X);

      if constexpr(FOURIER_MODE_SAFE)
      {
        checkNaN(x, "ola: x");
        checkNaN(X, "ola: X (avant traitement)");
      }

      config.traitement_freq(X);


      plan.step(X, x2, false);
      //x2.setZero(2*Ne);

      // Overlap-add
      // svg de taille Ne
      // si N_zeros = Ne, c'est tout le bloc qui est mis à jour
      svg.tail(N_zeros) += x2.head(N_zeros);
      y = svg;
      svg = x2.tail(Ne);
      cnt_ech += Ne;


      if constexpr(FOURIER_MODE_SAFE)
      {
        checkNaN(X, "ola: X (après traitement)");
        checkNaN(x2, "ola: x2");
        checkNaN(svg, "ola: svg");
      }
    }
    // OLA recouvrement facteur 2, et fenêtrage
    else
    {
      // 1) prev + nv
      padded.tail(Ne / 2) = x.head(Ne / 2);
      padded.tail(Ne) *= fenêtre;

      plan.step(padded, X);
      config.traitement_freq(X);
      plan.step(X, x2, false);

      // Ne + Nz = N
      // Overlap-add
      svg.segment(Ne - N_zeros, N_zeros) += x2.head(N_zeros);

      {
        last.tail(Ne/2)  += svg.head(Ne/2) / 2.0f;
        if(cnt_ech >= 0)
        {
          y = last;
        }
        else
          y.resize(0);
          //sortie.dispatch(last);

        last.head(Ne/2) = svg.tail(Ne/2) / 2.0f;
        last.tail(Ne/2).setZero();
      }
      tsd_assert(x2.rows() == N);
      svg = x2.tail(Ne);//x2.segment(N_zeros, Ne);

      cnt_ech += Ne / 2;

      // 2) nv seul
      padded.tail(Ne) = x * fenêtre;//.segment(N_zeros, Ne) = x;

      //energie = padded.abs2().sum();

      plan.step(padded, X);
      config.traitement_freq(X);
      plan.step(X, x2, false);

      // Overlap-add
      //svg.segment(Ne - N_zeros, N_zeros) += x2.head(N_zeros);
      svg.tail(N_zeros) += x2.head(N_zeros);
      last += svg / 2.0f;
      svg = x2.segment(N_zeros, Ne);

      cnt_ech += Ne / 2;

      // SVG pour prochain calcul
      //padded.segment(N_zeros, Ne/2)
      //            = padded.segment(N_zeros+Ne/2, Ne/2);
      padded.segment(N_zeros, Ne/2) = x.tail(Ne/2);

      //energie_bloc = energie;
    }
  }
};


std::tuple<sptr<Filtre<cfloat, cfloat, FiltreFFTConfig>>, int> filtre_fft(const FiltreFFTConfig &config)
{
  auto res = std::make_shared<OLA<cfloat>>();
  res->configure(config);
  return {res, res->N};
}


/*std::tuple<sptr<Filtre<float, float, FiltreFFTConfig>>, int> filtre_rfft(const FiltreFFTConfig &config)
{
  auto res = std::make_shared<OLA<float>>();
  res->configure(config);
  return {res, res->N};
}*/






template<typename T>
struct FiltreFFTRIF: FiltreGen<T>
{
  OLA<cfloat> ola;
  ArrayXcf H;


  FiltreFFTRIF(const ArrayXf &h)
  {
    FiltreFFTConfig ola_config;
    ola_config.nb_zeros_min           = h.rows();
    ola_config.traitement_freq    = [&](ArrayXcf &X)
    {
      X *= H;
    };
    ola.configure(ola_config);

    ArrayXf h2 = ArrayXf::Zero(ola.N);
    h2.tail(h.rows()) = h;
    //h2.head(h.rows()) = h;
    H = fft(h2);
    H *= std::sqrt(ola.N);

    if constexpr(FOURIER_MODE_SAFE)
    {
      checkNaN(h, "RIFFFT: h");
      checkNaN(h2, "RIFFFT: h2");
      checkNaN(H, "RIFFFT: H");
    }
  }
  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    ArrayXcf y1, x1;
    x1 = x;
    ola.step(x1, y1);
    if constexpr(est_complexe<T>())
      y = y1;
    else
      y = y1.real();
  }
};


}

namespace tsd::filtrage {

// Filtre OLA
template<typename T>
sptr<FiltreGen<T>> filtre_rif_fft(const ArrayXf &h)
{
  return std::make_shared<tsd::fourier::FiltreFFTRIF<T>>(h);
}

}

namespace tsd::fourier {




// Alignement de deux signaux suivant un délais variable
struct AlignementSignal::Impl
{
  ArrayXXcf residu[2];
  int N = 0;
  int delais_en_cours = 0;

  // Tampon de sortie synchronisé (matrice complexe à deux colonnes)
  sptr<SinkGen<cfloat>> sortie[2];

  Impl()
  {
    configure(1);
  }
  int configure(int N)
  {
    this->N = N;
    delais_en_cours = 0;
    sortie[0] = tampon_création<cfloat>(N, [&](const ArrayXcf &)
        {

        });
    sortie[1] = tampon_création<cfloat>(N, [&](const ArrayXcf &)
            {

            });
    residu[0].resize(0,1);
    residu[1].resize(0,1);
    return 0;
  }
  void step(const ArrayXXcf &x, const ArrayXXcf &y, int delais)
  {
    // Principe : on a un grand tampon de sortie synchronisé
    // Deux nouveaux blocs en entrée + un résidu sur une des entrées
    int n0 = x.rows();
    tsd_assert(y.rows() == n0);

    // Si pas de changement de délais,
    // copie du résidu + nv données vers buffer de sortie + stockage résidu

    // Combien on peut envoyer ?
    //  - n0 sur celui qui n'a pas de résidu
    //  - n0-r sur l'autre
    // n0 >= r => ok
    // n0 < r => on ne peut pas envoyer tout le résidu
    // celui qui n'a pas de résidu
    //unsigned int nb;
    //if(residu[0].rows() > 0)
    //  nb = y.rows();


    // Ou alors plus simple :
    // On pose les deux vers les résidus
    // Et on vide les résidus


    // n[0] : nombre de données à prendre sur x
    // n[1] : nombre de données à prendre sur y
    int n[2] = {n0, n0};

    // Cas nominal : prends autant sur x et y ?

    // Consomme un peu moins sur la première entrée
    if(delais > delais_en_cours)
    {
      n[1] -= (delais - delais_en_cours);
      if(n[1] < 0)
      {
        n[1] = 0;
        delais_en_cours += n0;
      }
      else
        delais_en_cours = delais;
    }
    // Sinon sur la deuxième entrée
    else if(delais < delais_en_cours)
    {
      n[0] -= (delais_en_cours - delais);
      if(n[0] < 0)
      {
        n[0] = 0;
        delais_en_cours -= n0;
      }
      else
        delais_en_cours = delais;
    }

    // n[1] = n0
    // n[0] = n0-7

    const ArrayXXcf *xin[2] = {&x, &y};



    for(auto i = 0u; i < 2; i++)
    {
      // printf("xin[%d] : %d elts, %f, %f\n", i, xin[i]->rows(), (*xin[i])(0,0).real(), (*xin[i])(1,0).real());
      // Transfert n[i] à partir de xin[i] vers résisu
      ArrayXXcf c(n[i] + residu[i].rows(), 1);
      //if(residu[i].rows() > 0)
      c << residu[i], xin[i]->block(n0 - n[i], 0, n[i], 1);
      //if(n[i] > 0)
      residu[i] = c;
      //printf("residu[%d]: %d elts: %f\n", i, residu[i].rows(), residu[i](1).real());
    }

    // Vidange des résidus
    int nmin = std::min(residu[0].rows(), residu[1].rows());
    ArrayXXcf s(nmin, 2);

    for(auto i = 0u; i < 2; i++)
    {
      for(auto j = 0; j < nmin; j++)
      {
        s(j,i) = residu[i](j);
      }
      residu[i] = residu[i].block(nmin, 0, residu[i].rows() - nmin, 1);
    }

    //sortie.on_event(s);
    sortie[0]->step(s.col(0));
    sortie[1]->step(s.col(1));
  }
};



// Alignement de deux signaux suivant un délais variable
/*sptr<ProcesseurConfigurable<cfloat, cfloat, int>> creation_aligneur()
{
  return std::make_shared<AlignementSignal>();
}*/

AlignementSignal::AlignementSignal()
{
  impl = std::make_shared<Impl>();
}

int AlignementSignal::configure(int N)
{
  return impl->configure(N);
}

void AlignementSignal::step(const ArrayXXcf &x, const ArrayXXcf &y, int delais)
{
  impl->step(x, y, delais);
}


int SpectrumConfig::Nf() const
{
  return BS / nsubs;
}

int SpectrumConfig::Ns() const
{
  if(sweep.active)
    return Nf() + (nsubs-1) * sweep.step;
  return Nf();
}

struct Spectrum: Filtre<cfloat,float,SpectrumConfig>
{
  sptr<FFTPlan> plan;
  // Fenêtre
  ArrayXf f;
  ArrayXf mag_moy, mag_cnt;
  ArrayXf masque;
  int cntmag = 0;

  // Dim FFT
  int Nf;
  // Dim spectre total
  int Ns;


  int configure_impl(const SpectrumConfig &config)
  {
    cntmag  = 0;
    Nf      = config.BS / config.nsubs;
    Ns      = Nf;

    masque = ArrayXf::Ones(Nf);
    if(config.sweep.masque_hf > 0)
    {
      masque.head(config.sweep.masque_hf).setZero();
      masque.tail(config.sweep.masque_hf).setZero();
    }
    if(config.sweep.masque_bf > 0)
      masque.segment(Nf/2 - config.sweep.masque_bf, 2 * config.sweep.masque_bf).setZero();

    if(config.sweep.active)
    {
      Ns = Nf + (config.nsubs-1) * config.sweep.step;
      mag_cnt = ArrayXf::Zero(Ns);
      for(auto i = 0; i < config.nsubs; i++)
        mag_cnt.segment(i * config.sweep.step, Nf) += masque;
      // pour éviter une division par zéro si le saut configuré est trop important
      mag_cnt = mag_cnt.cwiseMax(1.0f).eval();
    }

    mag_moy = ArrayXf::Zero(Ns);
    f       = tsd::filtrage::fenetre(config.fenetre, Nf, false);


    // Normalise les coefficients de la fenêtre, afini que l'énergie totale ne soit pas modifiée
    // (du moins, si le signal n'est pas  corrélé avec la fenêtre)
    f = std::sqrt(Nf / f.abs2().sum()) * f;

    // Note : si fenêtre rectangulaire, f est inchangée.

    auto n2 = f.abs2().sum();
    auto err_rel = std::abs(n2 - Nf) / Nf;
    if(err_rel > 1e-3)
    {
      msg("Après normalisation de la fenêtre : nrm={}, devrait être égal à Nf={}", n2, Nf);
      msg_avert("Erreur trop importante : {} %", err_rel * 100);
    }

    /*{
      tsd::vue::Figures fig;
      fig.subplot().plot(f);
      fig.subplot().plot_psd(f);
      fig.enregistrer("./build/spfen.png");
    }*/

    //msg("Fenêtre spectre: {}", (int) config.fenetre);
    plan = config.plan ? config.plan : fftplan_defaut();
    plan->configure(Nf, true, false);
    return 0;
  }

  ArrayXcf yt;
  void step(const Eigen::Ref<const Vecteur<cfloat>> x, ArrayXf &y)
  {
    auto &config = Configurable<SpectrumConfig>::config;
    // TODO (optionally): overlap 1/2
    tsd_assert_msg(x.rows() == config.BS, "Spectrum : dimension invalide.");

    moniteur_spectrum.commence_op();
    //msg("spec : début.");

    //msg("nsubs: {}", config.nsubs);

    // Mode multi-threadé
    if(config.nsubs > 1)
    {
      ArrayXf yft[config.nsubs];
#     if LIBTSD_USE_OMP
#     pragma omp parallel for
#     endif
      for(auto i = 0; i < config.nsubs; i++)
      {
        ArrayXcf tmp;
        //double en1 = x.segment(i * Nf, Nf).
        plan->step(x.segment(i * Nf, Nf) * f, tmp, true);


        /*{
          float E1 = x.segment(i * Nf, Nf).abs2().sum();
          float E2 = tmp.abs2().sum() / Nf;
          msg("E1 = {}, E2 = {}", E1, E2);
        }*/


        yft[i] = fftshift(tmp.abs2());
      }
      // Si balayage, il faut faire un mag_moy un peu différent...
      // On suppose dans tous les cas que le paquet reçu (dim = BS)
      // contient l'ensemble du balayage


      if(config.sweep.active)
      {
        for(auto i = 0; i < config.nsubs; i++)
          mag_moy.segment(i * config.sweep.step, Nf) += yft[i] * masque;
      }
      else
      {
        for(auto i = 0; i < config.nsubs; i++)
          mag_moy += yft[i];
      }

    }
    else
    {
      plan->step(x * f, yt, true);
      mag_moy += fftshift(yt.abs2());
    }
    cntmag++;
    if(cntmag == config.nmeans)
    {
      //mag_moy /= (config.nmeans * config.nsubs * std::sqrt((float) Nf));
      // Carré de sqrt(Nf) car ici on est en puissance
      // Facteur 1/Nf pour avoir une FFT normalisée
      mag_moy /= (config.nmeans * config.nsubs * Nf);
      if(config.sweep.active)
        mag_moy /= mag_cnt;
      y = 10 * log10(mag_moy + std::numeric_limits<float>::min());


#     if 0
      // Attention :
      //  - En temporel : puissance = énergie / durée
      //  - En fréquentiel :
      //     - Découpage en N sous blocs (N = nmeans * nsub)
      //     - Calcul de la TF normalisée sur chaque sous bloc


      auto puissance_temporel = [](const ArrayXcf &x, float fs, float R, float g)
      {
        float ts = 1 / fs;
        ArrayXf xe = x.abs2() / (g * g);
        float energie_totale = ts * xe.sum() / R;
        float T = ts * xe.rows();
        float puissance = energie_totale / T;
        float p_dBm = 10*log10(puissance) + 30;
        return p_dBm;
      };

      auto puissance_frequentiel = [](const ArrayXf &X, float fs, float R, float g)
      {
        ArrayXf P_Watt_par_bin = X / (X.rows() * g * g * R);
        float puissance_totale = P_Watt_par_bin.sum();
        float p_dBm = 10*log10(puissance_totale) + 30;
        return p_dBm;
      };



      {
        float fs = 20e6;
        float R = 50;
        float g = 21.8; // En lsb / Volt

        msg("Puissance instantanée - temporel :    {:.1f} dBm.", puissance_temporel(x, fs, R, g));
        msg("Puissance instantanée - fréquentiel : {:.1f} dBm.", puissance_frequentiel(mag_moy, fs, R, g));

        //msg("nsubs = {}, nmeans = {}", config.nsubs, config.nmeans);
      }
#     endif

      mag_moy.setZero();
      cntmag = 0;
    }
    else
      y.resize(0);

    //msg("spec : fin.");
    moniteur_spectrum.fin_op();
  }
};

sptr<Filtre<cfloat,float,SpectrumConfig>> rt_spectrum(const SpectrumConfig &config)
{
  auto res = std::make_shared<Spectrum>();
  res->configure(config);
  return res;
}



// Adapté depuis le script Scilab czt.sci
ArrayXcf czt(IArrayXcf x, int m, cfloat W, cfloat z0)
{
   int n = x.rows();
   int nm = std::max(n,m);

   // create sequence h(n)=[w*exp(-j*ϕ)]**(-n*n/2)
   ArrayXcf h(2*nm-1);

   for(auto i = 0; i < nm; i++)
     h(i) =  std::pow(W, -0.5f * i * i);
   for(auto i = nm; i < 2*nm-1; i++)
     h(i) = h(nm - 1 - (i - nm));
   // h(nm+1:2*nm-1) = h(nm-1:-1:1);

   tsd_assert(!x.hasNaN());

   ArrayXcf g(n);
   for(auto i = 0; i < n; i++)
   {
     tsd_assert(std::abs(h(nm + i - 1)) != 0);
     g(i) = x(i) * std::pow(z0, -i) / h(nm + i - 1);
   }


   tsd_assert(!g.hasNaN());

   // convolve h(n) and g(n)
   ArrayXcf hc(m + n - 1);
   hc.head(m)   = h.segment(nm-1, m);
   hc.tail(n-1) = h.segment(nm-n,n-1);

   ArrayXcf gc = ArrayXcf::Zero(m + m - 1);
   gc.head(n) = g;

   tsd_assert(!hc.hasNaN());
   tsd_assert(!gc.hasNaN());

   auto hcg = ifft(fft(hc) * fft(gc));

   tsd_assert(!hcg.hasNaN());

   // preserve m points and divide by h(n)
   return hcg.head(m) / h.segment(nm-1,m);
}

template<typename T>
  Vecteur<T> reechan_freq(const Vecteur<T> &x, float lom)
{
  if(lom == 1)
    return x;

  int n  = x.rows();
  int n2 = std::round(n * lom);
  auto X = fft(x);

  if(lom > 1)
  {
    decltype(X) X2;
    X2.setZero(n2);
    X2.head(n/2) = X.head(n/2);
    X2.tail(n/2) = X.tail(n/2);
    decltype(X) x = ifft(X2) * std::sqrt(lom);
    if constexpr(est_complexe<T>())
      return x;
    else
      return x.real();
  }
  else
  {
    decltype(X) X2;
    X2.setZero(n2);
    X2.head(n2/2) = X.head(n2/2);
    X2.tail(n2/2) = X.tail(n2/2);
    //auto x = ifft(X2) * std::sqrt(lom);
    if constexpr(est_complexe<T>())
      return ifft(X2) * std::sqrt(lom);
    else
      return (ifft(X2) * std::sqrt(lom)).real();
  }
}


template
  Vecteur<float> reechan_freq(const Vecteur<float> &x, float lom);

template
  Vecteur<cfloat> reechan_freq(const Vecteur<cfloat> &x, float lom);


template
Vecteur<float> délais<float>(const Vecteur<float> &x, float delay);

template
Vecteur<cfloat> délais<cfloat>(const Vecteur<cfloat> &x, float delay);






}

namespace tsd::filtrage {
template
sptr<FiltreGen<float>> filtre_rif_fft<float>(const ArrayXf &h);

template
sptr<FiltreGen<cfloat>> filtre_rif_fft<cfloat>(const ArrayXf &h);
}

namespace tsd::tf {
ArrayXXf periodogramme_tfd(const ArrayXcf &x, int N)
{
  std::vector<ArrayXf> mags;

  //int N = 512;
  //int N = 4096;
  // 48k/50 = 1000
  //int N = 4096;

  //int N = 4096;

  tsd::fourier::FiltreFFTConfig ola_config;
  ola_config.avec_fenetrage     = true;
  ola_config.dim_blocs_temporel = N;
  ola_config.nb_zeros_min           = 0;
  ola_config.traitement_freq    = [&](ArrayXcf &X)
  {
    ArrayXf xa = 10 * (X.abs2() + 1e-20f).log10().head(X.rows()/2);
    mags.push_back(xa);
  };

  auto [ola, N2] = tsd::fourier::filtre_fft(ola_config);

  ArrayXcf x2 = x;

  ola->step(x2);

  ArrayXXf M(N2/2, mags.size());
  for(auto i = 0u; i < mags.size(); i++)
    M.col(i) = mags[i];

  return M.transpose();
}



////////////////////////////////////////////////////////////




}


