#include "tsd/tsd-all.hpp"
#include "tsd/moniteur-cpu.hpp"


#include <bit>
#include <limits>


using namespace std;

const auto FOURIER_MODE_SAFE = non;

template<typename T>
void checkNaN(const Vecteur<T> &x, cstring s)
{
  si(x.hasNaN())
    échec("{} : has NaN.", s);
}

namespace tsd::fourier
{

template<typename T>
auto tfr_rotation_ref(entier N, entier signe = -1)
{
  Vecteur<complex<T>> y(N);
  pour(auto k = 0; k < N; k++)
    y(k) = std::polar<T>(1.0f, signe*2.0*π*k/N);
  retourne y;
}

template<typename T>
auto tfr_rotation_rapide(entier n, entier signe = -1)
{
  Vecteur<complex<T>> x(n);
  cdouble r = 1, w0 = std::polar<double>(1.0, (signe*2*π)/n);

  // TODO : faire les calculs sur un seul quadrant
  pour(auto i = 0; i < n; i++)
  {
    x(i) = r;
    r *= w0;
  }

  retourne x;
}

template<typename T>
auto tfr_rotation(entier n, entier signe = -1)
{
retourne tfr_rotation_rapide<T>(n, signe);
}

template<typename T>
auto itfr_rotation(entier n)
{
  retourne tfr_rotation<T>(n, 1);
}


/** @brief Calcul TFR Cooley–Tukey (N doit être une puissance de 2) */
template<typename T, typename Trotation, bouléen avant>
void tfr_radix2(Vecteur<complex<T>> &X,
                const Vecteur<complex<T>> &x,
                Vecteur<complex<T>> &scratch,
                const Vecteur<complex<Trotation>> &rotations)
{
  soit N               = x.rows();
  soit iteration_paire = ((N & 0x55555555) != 0);

  scratch.resize(N);
  X.resize(N);
  assertion_safe(rotations.rows() == N, "La matrice de twiddles n'est pas initialisée.");

  si(N == 1)
  {
    X(0) = x(0);
    retourne;
  }

  soit E = x.data();
  complex<T> *Xp, *Xp2, *Xstart;

  // Commence par calculer les TFR de plus bas niveau
  // puis remonte petit à petit
  pour(entier n = 1; n < N; n *= 2)
  {
    // Premiére itération : E = x, => Xp, Xs = scratch
    // Deuxième itération : E = scratch, Xs = X
    // ....
    Xstart = iteration_paire ? scratch.data() : X.data();
    soit pas = N / (2 * n);
    Xp  = Xstart;
    Xp2 = Xstart + N / 2;

    pour(entier k = 0; k < n; k++)
    {
      soit rot = rotations(k * pas);

      si constexpr(!avant)
          rot = conj(rot);

      const soit tr = rot.real(), ti = rot.imag();
      pour(entier m = 0; m < pas; m++)
      {
        const soit g = *(E + pas), e = *E;
        const soit gx = g.real(), gy = g.imag();
        const complex<T> p(tr * gx - ti * gy, tr * gy + ti * gx);
        *Xp++  = e + p;
        *Xp2++ = e - p;
        E++;
      }
      E += pas;
    }
    E = Xstart;
    iteration_paire = !iteration_paire;
  }

  // On a divisé par 2 à chaque itération au lieu de diviser par srqt(2)
  X /= sqrt((float) N);
}



// A SUPPRIMER
/** Calcule la TFR d'un vecteur réel :
 *  n échantillons réels -> n/2 échantillons complexes (fréquences positives) **/
/** Algorithme : numerical recipes, p512 */
template<typename T>
Vecteur<complex<T>> rtfr_calcule(const Vecteur<T> &x)
{
  using cT = complex<T>;

  soit n = x.rows();
  Vecteur<cT> x2(n/2), Xt(n/2), X(n/2);

  pour(auto i = 0; i < n / 2; i++)
  {
    x2(i).real(x(2*i));
    x2(i).imag(x(2*i+1));
  }

  // Compute n/2 FFT points in X
  Xt = fft(x2);

  const cT j2(0,0.25), r2(0.25,0);

  Vecteur<cT> rot(n/2+1);

  pour(auto k = 0; k <= n/2; k++)
  {
    rot(k).real(cos(-(2*π*k)/n));
    rot(k).imag(sin(-(2*π*k)/n));
  }

  pour(auto i = 0; i <= n / 2; i++)
  {
    cT X1, X2, X3;

    si(i == n/2)
      X1 = Xt(0);
    sinon
      X1 = Xt(i);

    si(i > 0)
      X2 = Xt(n/2-i);
    sinon
      X2 = Xt(0);

    X3 = r2 * (X1 + conj(X2)) - j2 *  (X1 - conj(X2)) * rot(i);

    si(i < n / 2)
      X(i) = X3;
    sinon
      X(0).imag(X3.real()); // Note: actually X(0).imag() = 0, but we store here X(n/2)
  }

  retourne X;
}


/** @brief Restore real samples from complex positive frequencies
 *  From n/2 complex values (positive frequencies), produces n real samples **/
/** Algorithm: numerical recipes, p512 */
template<typename T>
Vecteur<T> irtfr_calcule(const Vecteur<complex<T>> &X) // Size n/2
{
  using cT = complex<T>;

  soit n = 2 * X.rows();

  //assertion(x.rows() == (entier) n, "Output IFFT size must be n (= 2 size of complex FFT)");

  Vecteur<T> x(n);
  Vecteur<cT> X2(n/2);
  cT J(0,1.0);

  pour(auto i = 0; i < n / 2; i++)
  {
    cT Xi, Xp;

    si(i == 0)
      Xp = X(0).imag();
    sinon
      Xp = X(n/2-i);

    si(i == 0)
      Xi = X(0).real();
    sinon
      Xi = X(i);

    cT Xe = (Xi + conj(Xp)),
       Xo = (Xi - conj(Xp))
          * polar<T>(1.0, (2*π*i)/n);


    X2(i) = Xe + J * Xo;
  }

  // Compute n/2 FFT points in X
  soit x2 = ifft(X2);

  // x2 is the same as x, but viewed as a complex vector of n/2 samples
  pour(auto i = 0; i < n / 2; i++)
  {
    x(2*i)   = x2(i).real();
    x(2*i+1) = x2(i).imag();
  }

  retourne x;
}



// twidlles  : n2
// chirp     : 2*n-1
Veccf tfr_czt_impl(const Veccf &x, entier n2, const Veccf &twiddles, const Veccf &chirp)
{
  soit n = x.rows();

  // chirp = exp(-2pi * i * [-(n-1)/n, -(n-2)/n, ... 0, ..., (n-1)/n])
  soit xp = Veccf::zeros(n2);
  xp.head(n) = x * chirp.tail(n);

  soit icp = Veccf::zeros(n2);
  icp.head(2*n-1) = chirp.conjugate();

  Veccf scratch(n2), y, y2, Xp, Xc;

  //y = ifft(fft(xp) * fft(icp))
  tfr_radix2<float, float, oui>(Xp, xp,       scratch, twiddles);
  tfr_radix2<float, float, oui>(Xc, icp,      scratch, twiddles);
  tfr_radix2<float, float, non>(y2, Xp * Xc,  scratch, twiddles);
  retourne y2.segment(n-1, n) * chirp.tail(n) * (sqrt((float)n2) / sqrt((float)n));
}

/** Conversion du résultat d'une tfr vers une itfr
 *  E.g. Y(n) = X(N-n) */
Veccf tfr2itfr(const Veccf &X)
{
  soit n = X.rows();
  Veccf Y(n);
  Y(0) = X(0);

  si((n & 1) == 0)
  {
    Y(n/2)              = X(n/2);
    Y.segment(1, n/2-1) = X.tail(n/2-1).reverse();
    Y.tail(n/2-1)       = X.segment(1, n/2-1).reverse();
  }
  sinon
  {
    Y.segment(1, n/2)   = X.tail(n/2).reverse();
    Y.tail(n/2)         = X.segment(1, n/2).reverse();
  }

  retourne Y;
}

template<typename T>
struct RTFRPlan: FiltreGen<T, complex<T>>
{
  using cT = complex<T>;

  entier n;
  sptr<FFTPlan> cplan;
  Veccf rotations;

  RTFRPlan(entier n)
  {
    configure(n);
  }
  void configure(entier n)
  {
    this->n = n;

    si(n > 0)
    {
      si((n & 1) == 0)
      {
        cplan      = tfrplan_création(n/2);
        rotations  = tfr_rotation<T>(n);
      }
      sinon
      {
        // si n est impair, on fait juste bêtement une FFT.
        cplan = tfrplan_création(n);
      }
    }
  }
  void step(const Vecteur<T> &x, Vecteur<cT> &y)
  {
    si(x.rows() != n)
      configure(x.rows());
    si((n & 1) == 0)
    {
      y.resize(n);
      Vecteur<cT> x2(n/2);

      pour(auto i = 0; i < n / 2; i++)
      {
        x2(i).real(x(2*i));
        x2(i).imag(x(2*i+1));
      }

      // Compute n/2 FFT points in X
      soit Xt = cplan->step(x2);

      const cT j2(0, 0.5 / sqrt(2)), r2(0.5 / sqrt(2), 0);

      pour(auto i = 0; i <= n / 2; i++)
      {
        cT X1, X2;

        si(i == n/2)
          X1 = Xt(0);
        sinon
          X1 = Xt(i);

        si(i > 0)
          X2 = Xt(n/2-i);
        sinon
          X2 = Xt(0);

        y(i) = r2 * (X1 + conj(X2)) - j2 *  (X1 - conj(X2)) * rotations(i);
      }
      csym_forçage(y);
    }
    sinon
    {
      Vecteur<cT> y1 = x.as_complex();
      cplan->step(y1, y);
    }
  }
};




struct TFRPlanDefaut: FiltreGen<cfloat>, FFTPlan
{
  bouléen normaliser  = oui, avant = oui;
  entier n = 0, n2 = 0;
  Veccf scratch, rotations, chirp;
  sptr<FFTPlan> sousplan;

  TFRPlanDefaut(entier n = -1, bouléen avant = oui, bouléen normalize = oui)
  {
    configure(n, avant, normalize);
  }

  void configure(entier n, bouléen avant, bouléen normalize)
  {
    this->n = n;
    this->avant = avant;

    si(n == -1)
      retourne;

    n2 = n;

    // Pas une puissance de 2 ?
    si((n & (n - 1)) != 0)
    {
      // n est pair
      si((n & 1) == 0)
      {
        // Décompose tant que pair
        sousplan = make_shared<TFRPlanDefaut>(n / 2, avant, normalize);
      }
      sinon
      {
        //msg("Avertissement : FFT sur n != 2^k (n = %d)", n);
        n2 = prochaine_puissance_de_2(2*n-1);
        // CZT
        soit t = square(linspace(-(n-1), n-1, 2*n-1)) / 2;
        t *= -2*π/n;
        // t = -2pi * [-(n-1)/n, -(n-2)/n, ... 0, ..., (n-1)/n]
        chirp = polar(t);
      }
    }
    //msg("plan fft : n = {}, n2 = {}", n, n2);
    scratch.resize(n2);
    rotations = tfr_rotation<float>(n2);
  }

  void step(const Veccf &x, Veccf &y)
  {
    step(x, y, this->avant);
  }

  void step(const Veccf &x, Veccf &y, bouléen avant)
  {
    assertion(x.rows() > 0);

    si((entier) n != x.rows())
      configure(x.rows(), this->avant, this->normaliser);

    si(n2 != n)
    {
      y = tfr_czt_impl(x, n2, rotations, chirp);
      si(!avant)
      {
        y = tfr2itfr(y);
      }
    }
    sinon
    {
      // Puissance de 2
      si((n & (n - 1)) == 0)
      {
        y.resize(n);
        si(avant)
          tfr_radix2<float, float, oui>(y, x, scratch, rotations);
        sinon
          tfr_radix2<float, float, non>(y, x, scratch, rotations);
      }
      sinon
      {
        Veccf xe(n/2), xo(n/2);
        pour(auto i = 0; i < n/2; i++)
        {
          xe(i) = x(2*i);
          xo(i) = x(2*i+1);
        }

        soit E = sousplan->step(xe, avant),
             O = sousplan->step(xo, avant);
        y.resize(n);

        si(avant)
        {
          y.head(n/2) = E + rotations.head(n/2) * O;
          y.tail(n/2) = E + rotations.tail(n/2) * O;
        }
        sinon
        {
          y.head(n/2) = E + rotations.head(n/2).conjugate() * O;
          y.tail(n/2) = E + rotations.tail(n/2).conjugate() * O;
        }
        si(normaliser)
          y *= 1 / sqrt(2.0f);
      }
    }
  }

};

fonction<sptr<FFTPlan>()> fftplan_defaut = []()
{
    retourne make_shared<TFRPlanDefaut>();
};


sptr<FFTPlan> tfrplan_création(entier n, bouléen avant, bouléen normalize)
{
  soit res = fftplan_defaut();
  si(n >= 0)
    res->configure(n, avant, normalize);
  retourne res;
}

sptr<FiltreGen<float, cfloat>> rtfrplan_création(entier n)
{
  retourne make_shared<RTFRPlan<float>>(n);
}


// Produit de corrélation effectué dans le domaine fréquentiel
template<typename T>
Vecteur<T> correlation_freq(const Vecteur<T> &X0, const Vecteur<T> &X1)
{
  soit n = X0.rows();
  assertion(n == X1.rows());
  Vecteur<T> Y(n);

  Y(0) = X0(0) * conj(X1(0));
  si(n > 1)
    Y.tail(n-1) = X0.tail(n-1).reverse() * X1.tail(n-1).reverse().conjugate();

  Y *= sqrt((float) n);
  retourne Y;
}


struct TFRCorrelateurBloc
{
  sptr<FFTPlan> plan = fftplan_defaut();

  /** @note Should pad x0 and x1 with K zeros si K first lags are examined
   *  (e.g. n zeroes si all lags are used).
   *
   *  Result: y(0): lag = 0
   *          y(1): lag = 1
   *          ....
   *          y(n-2): lag = -2
   *          y(n-1): lag = -1
   */
  Veccf step(const Veccf &x0, const Veccf &x1)
  {
    soit &x1p = (x1.rows() == 0) ? x0 : x1;

    assertion_msg(x0.rows() == x1p.rows(),
        "the two input vectors should have the dimension {} != {}.",
        x0.rows(), x1p.rows());

    // (1) Compute FFT on x0 and x1
    soit X0 = plan->step(x0, oui),
         X1 = plan->step(x1p, oui);

    soit X2 = correlation_freq(X0, X1);
    retourne plan->step(X2, non);

    // y(0): lag = 0
    // y(1): lag = 1
    // ...
    // y(n-1): lag = n-1 == -1 from the periodicity hypothesis.

    //retourne y;//.conjugate();
  }
};


tuple<Vecf, Veccf> ccorr(const Veccf &x0, const Veccf &x1)
{
  TFRCorrelateurBloc fc;
  soit m    = x0.rows();
  soit lags = linspace(0, m-1, m);
  retourne {lags, fc.step(x0, x1) / m};
}

tuple<Vecf, Veccf> xcorrb(const Veccf &x, const Veccf &y, entier m)
{
  TFRCorrelateurBloc fc;

  soit n = x.rows();

  si(m < 0)
    m = n;

  soit yp = (y.rows() == 0) ? x : y;

  soit x2 = Veccf::zeros(m + n + m),
       y2 = x2;

  x2.segment(m, n) = x;
  y2.segment(m, n) = yp;

  soit r = fc.step(x2, y2);

  Veccf res(2*m-1);
  res.tail(m)   = r.head(m) / n;
  res.head(m-1) = r.tail(m-1) / n;

  soit lags = linspace(-(m-1), m-1, 2*m-1);

  retourne {lags, res};
}

tuple<Vecf, Veccf> xcorr(const Veccf &x, const Veccf &y, entier m)
{
  soit n = x.rows();

  si(m < 0)
    m = n;

  // Corrélation biaisée
  soit [lags, zb] = xcorrb(x, y, m);

  si(m > 1)
  {
    zb.head(m-1) /= (linspace(n-(m-1), n-1, m-1) / n).as<cfloat>();
    zb.tail(m-1) /= (linspace(n-1, n-(m-1), m-1) / n).as<cfloat>();
  }

  retourne {lags, zb};
}


Vecf coherence(const Veccf &x, const Veccf &y)
{
  soit X = fft(x), Y = fft(y);
  retourne abs(X * Y.conjugate()) / (abs(X) * abs(Y));
}


template<typename T>
Vecteur<complex<T>> delais_fractionnaire_c(const Vecteur<complex<T>> &x, float τ)
{
  soit n = 2 * x.rows();

  // (1) Complète par des zéros à gauche et à droite
  soit x2 = Vecteur<complex<T>>::zeros(n);
  x2.segment(n/4, n/2) = x;

  // Modulation dans le domaine fréquentiel
  soit X = fft(x2);

  soit rot = Veccf::int_expr(n,
    IMAP(polar(1.0f, -2*π_f*i*τ/n + π_f*τ)));

  X *= fftshift(rot);

  // Retire le padding
  retourne ifft(X).segment(n/4, n/2);
}



// PB : ne marche pas si T = complexe
template<typename T = float>
Vecteur<T> delais_fractionnaire_r(const Vecteur<T> &x, float τ)
{
  soit n = 2 * x.rows();

  // (1) Padd x with zeroes at begin and end
  soit x2 = Vecteur<T>::zeros(n);
  x2.segment(n/4, n/2) = x;

  // A delay in time domain is a modulation in frequency domain
  soit X = rtfr_calcule(x2);

  complex<T> rot(1.0, 0), dphi;

  dphi = polar<T>((T) 1.0, (T) - τ * 2 * π / n);

  // TODO
  // Attention X(0).real -> freq(0)
  //           X(0).imag -> freq(n/2)
  pour(auto i = 0; i < n / 2; i++)
  {
    si(i == 0)
      X(i).real((X(i).real() * rot).real());
    sinon
      X(i) *= rot;

    rot  *= dphi;
  }

  X(0).imag((X(0).imag() * rot).real());

  // With depading
  retourne irtfr_calcule(X).segment(n/4, n/2).clone();
}

template<typename T = float>
Vecteur<T> delais_entier(const Vecteur<T> &x, entier τ)
{
  si(τ == 0)
    retourne x;

  soit n = x.rows();
  soit y = Vecteur<T>::zeros(n);

  // Delais positif
  si(τ > 0)
    y.tail(n - τ) = x.head(n - τ);
  // Delais négatif
  sinon
    y.head(n + τ) = x.tail(n + τ);

  retourne y;
}


template<typename T>
  Vecteur<T> délais(const Vecteur<T> &x, float delais)
{
  si(floor(delais) != delais)
  {
    si constexpr (est_complexe<T>())
      retourne delais_fractionnaire_c(x, delais);
    sinon
      retourne delais_fractionnaire_r(x, delais);
  }
  sinon
    retourne delais_entier<T>(x, (entier) delais);
}




//[C,Nf,Nz] = complexite(fe, fs, N, Ne)

// Complexité, par échantillon d'entrée, en FLOPS
// M  : taille de filtre ou motif
// Ne : taille de bloc d'entrée
void ola_complexité(entier M, entier Ne, float &C, entier &Nf, entier &Nz)
{
  Nf = prochaine_puissance_de_2(Ne + M - 1);
  Nz = Nf - Ne;
  C = (1.0f / Ne) * 2 * 5 * Nf * log(1.0f*Nf) / log(2.0f);
}

void ola_complexité_optimise(entier M, float &C_, entier &Nf_, entier &Nz_, entier &Ne_)
{
  // Nf doit être au minimum égal à M
  soit kmin = (entier) ceil(log(M) / log(2));

  pour(auto k = kmin; (k < kmin + 20) && (k < 31); k++)
  {
    entier Nf, Nz, Ne = (1 << k) - (M - 1);

    float C;
    ola_complexité(M, Ne, C, Nf, Nz);
    //msg("k = {} -> Ne = {}, C = {}", k, Ne, C);
    si((k == kmin) || (C < C_))
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
  Veccf padded, X, last, x2, svg;
  // TODO !
  TFRPlanDefaut plan;
  Vecf fenêtre;

  /** Nb ech. TFD */
  entier N = 0,

  /** Nb zéros insérés */
         N_zeros = 0,

  /** Nb ech. d'entrée / sortie */
         Ne = 0;

  /** Nombre d'échantillons en cours */
  int64_t cnt_ech = 0;

  sptr<SinkGen<cfloat>> tampon;
  bouléen tampon_vide = oui;

  vector<Veccf> bo;

  void configure_impl(const FiltreFFTConfig &config)
  {
    Ne      = config.dim_blocs_temporel;

    si(!config.traitement_freq)
      échec("configuration OLA : traitement fréquentiel non précisé.");

    si(Ne <= 0)
      Ne = 512;

    // Attention, post-affectation
    //Configurable<FiltreFFTConfig>::config.dim_blocs_temporel = Ne;

    N       = prochaine_puissance_de_2(Ne + config.nb_zeros_min);
    N_zeros = N - Ne;

    // Position initiale
    cnt_ech = - ((entier) Ne) / 2;

    msg("OLA init: Ne = {} (nb échan entrée), N = {} (dim TFD), Nz = {} (nb zéros).", Ne, N, N_zeros);

    padded.setZero(N);
    last.setZero(Ne);
    svg.setZero(Ne);

    /*si((config.facteur_recouvrement != 1) && (config.facteur_recouvrement != 2))
    {
      msg_erreur("OLA::configure : facteur de recouvrement : doit etre 2.");
      retourne -1;
    }*/

    si(config.avec_fenetrage)
    {
      fenêtre = tsd::filtrage::fenêtre("hn", Ne, non);
      assertion(fenêtre.rows() == Ne);
    }

    //plan  = creation_fft_plan(N, oui);
    //iplan = creation_fft_plan(N, non);


    tampon_vide = oui;
    tampon = tampon_création<cfloat>(Ne,
        [&](const Vecteur<cfloat> &x)
        {
          Veccf y;
          step_interne(x, y);
          bo.push_back(y);
        });
  }

  void step(const Veccf &x, Veccf &y)
  {
    si(tampon_vide && (x.rows() == Ne))
    {
      step_interne(x, y);
      retourne;
    }
    tampon_vide = non;
    tampon->step(x);
    soit n = 0;
    pour(auto &b: bo)
      n += b.rows();
    y.resize(n);
    soit i = 0;
    pour(auto &b: bo)
    {
      y.segment(i, b.rows()) = b;
      i += b.rows();
    }
    bo.clear();
  }



  void step_interne(const Veccf &x, Veccf &y)
  {
    soit &config = Configurable<FiltreFFTConfig>::config;

    assertion(x.rows() == Ne);

    si(!config.avec_fenetrage)
    {
      // OLA simple (pas de recouvrement sur les tampons d'entrée,
      // recouvrement sur les tampons de sorties)

      //msg("OLA/STI : Ne = {}, padded.rows()={}", Ne, padded.rows());

      padded.tail(Ne) = x;


      // Conso CPU des plans : ~ 20 % chacun

      plan.step(padded, X);

      si constexpr(FOURIER_MODE_SAFE)
      {
        checkNaN(x, "ola: x");
        checkNaN(X, "ola: X (avant traitement)");
      }

      config.traitement_freq(X);

      plan.step(X, x2, non);

      // Overlap-add
      // svg de taille Ne
      // si N_zeros = Ne, c'est tout le bloc qui est mis à jour
      svg.tail(N_zeros) += x2.head(N_zeros);
      y = svg.clone();
      svg.copie(x2.tail(Ne));
      cnt_ech += Ne;


      si constexpr(FOURIER_MODE_SAFE)
      {
        checkNaN(X, "ola: X (après traitement)");
        checkNaN(x2, "ola: x2");
        checkNaN(svg, "ola: svg");
      }
    }
    // OLA recouvrement facteur 2, et fenêtrage
    sinon
    {
      // 1) prev + nv
      padded.tail(Ne / 2) = x.head(Ne / 2);
      padded.tail(Ne) *= fenêtre.as_complex();

      plan.step(padded, X);
      config.traitement_freq(X);
      plan.step(X, x2, non);

      // Ne + Nz = N
      // Overlap-add
      svg.segment(Ne - N_zeros, N_zeros) += x2.head(N_zeros);

      {
        last.tail(Ne/2)  += svg.head(Ne/2) / 2.0f;
        si(cnt_ech >= 0)
          y = last.clone();
        sinon
          y.resize(0);

        last.head(Ne/2) = svg.tail(Ne/2) / 2.0f;
        last.tail(Ne/2).setZero();
      }
      assertion(x2.rows() == N);
      svg = x2.tail(Ne).clone();

      cnt_ech += Ne / 2;

      // 2) nv seul
      padded.tail(Ne) = x * fenêtre.as_complex();

      plan.step(padded, X);
      config.traitement_freq(X);
      plan.step(X, x2, non);

      // Overlap-add
      svg.tail(N_zeros) += x2.head(N_zeros);
      last += svg / 2.0f;
      svg = x2.segment(N_zeros, Ne);

      cnt_ech += Ne / 2;

      // SVG pour prochain calcul
      padded.segment(N_zeros, Ne/2) = x.tail(Ne/2);

    }
  }
};


tuple<sptr<Filtre<cfloat, cfloat, FiltreFFTConfig>>, entier> filtre_fft(const FiltreFFTConfig &config)
{
  soit res = make_shared<OLA<cfloat>>();
  res->configure(config);
  retourne {res, res->N};
}





template<typename T>
struct FiltreFFTRIF: FiltreGen<T>
{
  OLA<cfloat> ola;
  Veccf H;

  FiltreFFTRIF(const Vecf &h)
  {
    FiltreFFTConfig ola_config;
    ola_config.nb_zeros_min       = h.rows(); // -1 ?
    ola_config.traitement_freq    = [&](Veccf &X)
    {
      X *= H;
    };
    ola.configure(ola_config);

    soit h2 = Vecf::zeros(ola.N);
    h2.tail(h.rows()) = h;
    H = fft(h2);
    H *= sqrt(ola.N);

    si constexpr(FOURIER_MODE_SAFE)
    {
      checkNaN(h,  "RIFFFT: h");
      checkNaN(h2, "RIFFFT: h2");
      checkNaN(H,  "RIFFFT: H");
    }
  }
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    y = real(ola.Filtre<cfloat, cfloat, FiltreFFTConfig>::step(x.as_complex()));
  }
};


}

namespace tsd::filtrage {

// Filtre OLA
template<typename T>
sptr<FiltreGen<T>> filtre_rif_fft(const Vecf &h)
{
  retourne make_shared<tsd::fourier::FiltreFFTRIF<T>>(h);
}

}

namespace tsd::fourier {




// Alignement de deux signaux suivant un délais variable
struct AlignementSignal::Impl
{
  Veccf residu[2];
  entier N = 0;
  entier delais_en_cours = 0;

  // Tampon de sortie synchronisé (matrice complexe à deux colonnes)
  sptr<SinkGen<cfloat>> sortie[2];

  Impl()
  {
    configure(1);
  }
  void configure(entier N)
  {
    this->N = N;
    delais_en_cours = 0;
    sortie[0] = tampon_création<cfloat>(N, [&](const Veccf &)
        {

        });
    sortie[1] = tampon_création<cfloat>(N, [&](const Veccf &)
            {

            });
    residu[0].resize(0);
    residu[1].resize(0);
  }
  void step(const Veccf &x, const Veccf &y, entier delais)
  {
    // Principe : on a un grand tampon de sortie synchronisé
    // Deux nouveaux blocs en entrée + un résidu sur une des entrées
    entier n0 = x.rows();
    assertion(y.rows() == n0);

    // si pas de changement de délais,
    // copie du résidu + nv données vers buffer de sortie + stockage résidu

    // Combien on peut envoyer ?
    //  - n0 sur celui qui n'a pas de résidu
    //  - n0-r sur l'autre
    // n0 >= r => ok
    // n0 < r => on ne peut pas envoyer tout le résidu
    // celui qui n'a pas de résidu
    //unsigned int nb;
    //si(residu[0].rows() > 0)
    //  nb = y.rows();


    // Ou alors plus simple :
    // On pose les deux vers les résidus
    // Et on vide les résidus


    // n[0] : nombre de données à prendre sur x
    // n[1] : nombre de données à prendre sur y
    entier n[2] = {n0, n0};

    // Cas nominal : prends autant sur x et y ?

    // Consomme un peu moins sur la première entrée
    si(delais > delais_en_cours)
    {
      n[1] -= (delais - delais_en_cours);
      si(n[1] < 0)
      {
        n[1] = 0;
        delais_en_cours += n0;
      }
      sinon
        delais_en_cours = delais;
    }
    // Sinon sur la deuxième entrée
    sinon si(delais < delais_en_cours)
    {
      n[0] -= (delais_en_cours - delais);
      si(n[0] < 0)
      {
        n[0] = 0;
        delais_en_cours -= n0;
      }
      sinon
        delais_en_cours = delais;
    }

    // n[1] = n0
    // n[0] = n0-7

    const Veccf *xin[2] = {&x, &y};


    pour(auto i = 0u; i < 2; i++)
    {
      // printf("xin[%d] : %d elts, %f, %f\n", i, xin[i]->rows(), (*xin[i])(0,0).real(), (*xin[i])(1,0).real());
      // Transfert n[i] à partir de xin[i] vers résisu
      //Tabcf c(n[i] + residu[i].rows(), 1);
      //si(residu[i].rows() > 0)

      //c << residu[i], xin[i]->block(n0 - n[i], 0, n[i], 1);
      //si(n[i] > 0)
      //residu[i] = c;


      residu[i] = vconcat(residu[i], xin[i]->segment(n0 - n[i], n[i]));

      //printf("residu[%d]: %d elts: %f\n", i, residu[i].rows(), residu[i](1).real());
    }

    // Vidange des résidus
    soit nmin = min(residu[0].rows(), residu[1].rows());
    Tabcf s(nmin, 2);

    pour(auto i = 0; i < 2; i++)
    {
      pour(auto j = 0; j < nmin; j++)
      {
        s(j,i) = residu[i](j);
      }
      residu[i] = residu[i].segment(nmin, residu[i].rows() - nmin);
    }

    sortie[0]->step(s.col(0));
    sortie[1]->step(s.col(1));
  }
};



// Alignement de deux signaux suivant un délais variable
/*sptr<ProcesseurConfigurable<cfloat, cfloat, entier>> creation_aligneur()
{
  retourne make_shared<AlignementSignal>();
}*/

AlignementSignal::AlignementSignal()
{
  impl = make_shared<Impl>();
}

void AlignementSignal::configure(entier N)
{
  impl->configure(N);
}

void AlignementSignal::step(const Tabcf &x, const Tabcf &y, entier delais)
{
  impl->step(x, y, delais);
}


entier SpectrumConfig::Nf() const
{
  retourne BS / nsubs;
}

entier SpectrumConfig::Ns() const
{
  si(sweep.active)
    retourne Nf() + (nsubs-1) * sweep.step;
  retourne Nf();
}

struct Spectrum: Filtre<cfloat,float,SpectrumConfig>
{
  sptr<FFTPlan> plan;
  // Fenêtre
  Vecf f, mag_moy, mag_cnt, masque;
  entier cntmag = 0,
  // Dim FFT
         Nf = 0,
  // Dim spectre total
         Ns = 0;


  void configure_impl(const SpectrumConfig &config)
  {
    cntmag  = 0;
    Nf      = config.BS / config.nsubs;
    Ns      = Nf;

    masque = Vecf::ones(Nf);
    si(config.sweep.masque_hf > 0)
    {
      masque.head(config.sweep.masque_hf).setZero();
      masque.tail(config.sweep.masque_hf).setZero();
    }
    si(config.sweep.masque_bf > 0)
      masque.segment(Nf/2 - config.sweep.masque_bf, 2 * config.sweep.masque_bf).setZero();

    si(config.sweep.active)
    {
      Ns = Nf + (config.nsubs-1) * config.sweep.step;
      mag_cnt = Vecf::zeros(Ns);
      pour(auto i = 0; i < config.nsubs; i++)
        mag_cnt.segment(i * config.sweep.step, Nf) += masque;
      // pour éviter une division par zéro si le saut configuré est trop important
      mag_cnt = mag_cnt.cwiseMax(1.0f);
    }

    mag_moy = Vecf::zeros(Ns);
    f       = tsd::filtrage::fenêtre(config.fenetre, Nf, non);


    // Normalise les coefficients de la fenêtre, afini que l'énergie totale ne soit pas modifiée
    // (du moins, si le signal n'est pas  corrélé avec la fenêtre)
    f = sqrt(Nf / abs2(f).somme()) * f;

    // Note : si fenêtre rectangulaire, f est inchangée.

    soit n2 = abs2(f).somme();
    soit err_rel = abs(n2 - Nf) / Nf;
    si(err_rel > 1e-3)
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

    //msg("Fenêtre spectre: {}", (entier) config.fenetre);
    plan = config.plan ? config.plan : fftplan_defaut();
    plan->configure(Nf, oui, non);
  }

  Veccf yt;
  void step(const Vecteur<cfloat> &x, Vecf &y)
  {
    soit &config = Configurable<SpectrumConfig>::config;
    // TODO (optionally): overlap 1/2
    assertion_msg(x.rows() == config.BS, "Spectrum : dimension invalide.");

    moniteur_spectrum.commence_op();
    //msg("spec : début.");

    //msg("nsubs: {}", config.nsubs);

    // Mode multi-threadé
    si(config.nsubs > 1)
    {
      Vecf yft[config.nsubs];
#     if LIBTSD_USE_OMP
#     pragma omp parallel pour
#     endif
      pour(auto i = 0; i < config.nsubs; i++)
      {
        soit tmp = plan->step(x.segment(i * Nf, Nf) * f);
        yft[i] = fftshift(abs2(tmp));
      }
      // si balayage, il faut faire un mag_moy un peu différent...
      // On suppose dans tous les cas que le paquet reçu (dim = BS)
      // contient l'ensemble du balayage

      si(config.sweep.active)
      {
        pour(auto i = 0; i < config.nsubs; i++)
          mag_moy.segment(i * config.sweep.step, Nf) += yft[i] * masque;
      }
      sinon
      {
        pour(auto i = 0; i < config.nsubs; i++)
          mag_moy += yft[i];
      }

    }
    sinon
    {
      plan->step(x * f, yt, oui);
      mag_moy += fftshift(abs2(yt));
    }
    cntmag++;
    si(cntmag == config.nmeans)
    {
      // Carré de sqrt(Nf) car ici on est en puissance
      // Facteur 1/Nf pour avoir une FFT normalisée
      mag_moy /= (config.nmeans * config.nsubs * Nf);
      si(config.sweep.active)
        mag_moy /= mag_cnt;
      y = pow2db(mag_moy + numeric_limits<float>::min());


#     if 0
      // Attention :
      //  - En temporel : puissance = énergie / durée
      //  - En fréquentiel :
      //     - Découpage en N sous blocs (N = nmeans * nsub)
      //     - Calcul de la TF normalisée sur chaque sous bloc


      soit puissance_temporel = [](const ArrayXcf &x, float fs, float R, float g)
      {
        float ts = 1 / fs;
        ArrayXf xe = x.abs2() / (g * g);
        float energie_totale = ts * xe.sum() / R;
        float T = ts * xe.rows();
        float puissance = energie_totale / T;
        float p_dBm = 10*log10(puissance) + 30;
        retourne p_dBm;
      };

      soit puissance_frequentiel = [](const ArrayXf &X, float fs, float R, float g)
      {
        ArrayXf P_Watt_par_bin = X / (X.rows() * g * g * R);
        float puissance_totale = P_Watt_par_bin.sum();
        float p_dBm = 10*log10(puissance_totale) + 30;
        retourne p_dBm;
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
    sinon
      y.resize(0);

    //msg("spec : fin.");
    moniteur_spectrum.fin_op();
  }
};

sptr<Filtre<cfloat,float,SpectrumConfig>> rt_spectrum(const SpectrumConfig &config)
{
  soit res = make_shared<Spectrum>();
  res->configure(config);
  retourne res;
}



// Adapté depuis le script Scilab czt.sci
Veccf czt(const Veccf &x, entier m, cfloat W, cfloat z0)
{
   soit n  = x.rows(),
        nm = max(n,m);

   // create sequence h(n)=[w*exp(-j*ϕ)]**(-n*n/2)
   Veccf h(2*nm-1);

   pour(auto i = 0; i < nm; i++)
     h(i) =  pow(W, -0.5f * i * i);
   pour(auto i = nm; i < 2*nm-1; i++)
     h(i) = h(nm - 1 - (i - nm));

   assertion(!x.hasNaN());

   Veccf g(n);
   pour(auto i = 0; i < n; i++)
   {
     assertion(abs(h(nm + i - 1)) != 0);
     g(i) = x(i) * pow(z0, -i) / h(nm + i - 1);
   }


   assertion(!g.hasNaN());

   // convolve h(n) and g(n)
   Veccf hc(m + n - 1);
   hc.head(m)   = h.segment(nm-1, m);
   hc.tail(n-1) = h.segment(nm-n,n-1);

   soit gc = Veccf::zeros(m + m - 1);
   gc.head(n) = g;

   assertion(!hc.hasNaN());
   assertion(!gc.hasNaN());

   soit hcg = ifft(fft(hc) * fft(gc));

   assertion(!hcg.hasNaN());

   // preserve m points and divide by h(n)
   retourne hcg.head(m) / h.segment(nm-1, m);
}

template<typename T>
  Vecteur<T> rééchan_freq(const Vecteur<T> &x, float lom)
{
  si(lom == 1)
    retourne x;

  soit n  = x.rows(),
       n2 = (entier) round(n * lom);

  soit X = fft(x);

  si(lom > 1)
  {
    decltype(X) X2;
    X2.setZero(n2);
    X2.head(n/2) = X.head(n/2);
    X2.tail(n/2) = X.tail(n/2);
    decltype(X) x = ifft(X2) * sqrt(lom);
    retourne real(x);
  }
  sinon
  {
    decltype(X) X2;
    X2.setZero(n2);
    X2.head(n2/2) = X.head(n2/2);
    X2.tail(n2/2) = X.tail(n2/2);
    retourne real(ifft(X2) * sqrt(lom));
  }
}


template
  Vecteur<float> rééchan_freq(const Vecteur<float> &x, float lom);

template
  Vecteur<cfloat> rééchan_freq(const Vecteur<cfloat> &x, float lom);


template
Vecteur<float> délais<float>(const Vecteur<float> &x, float delay);

template
Vecteur<cfloat> délais<cfloat>(const Vecteur<cfloat> &x, float delay);






}

namespace tsd::filtrage {
template
sptr<FiltreGen<float>> filtre_rif_fft<float>(const Vecf &h);

template
sptr<FiltreGen<cfloat>> filtre_rif_fft<cfloat>(const Vecf &h);
}

namespace tsd::tf {
Tabf periodogramme_tfd(const Veccf &x, entier N)
{
  vector<Vecf> mags;

  //entier N = 512;
  //entier N = 4096;
  // 48k/50 = 1000
  //entier N = 4096;

  //entier N = 4096;

  tsd::fourier::FiltreFFTConfig ola_config;
  ola_config.avec_fenetrage     = oui;
  ola_config.dim_blocs_temporel = N;
  ola_config.nb_zeros_min       = 0;
  ola_config.traitement_freq    = [&](Veccf &X)
  {
    soit xa = 10 * (log10(abs2(X) + 1e-20f)).head(X.rows()/2);
    mags.push_back(xa);
  };

  soit [ola, N2] = tsd::fourier::filtre_fft(ola_config);

  ola->step(x);

  Tabf M(N2/2, mags.size());
  pour(auto i = 0u; i < mags.size(); i++)
    M.col(i) = mags[i];

  retourne M.transpose();
}



}


