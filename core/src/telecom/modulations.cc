#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <iostream>
#include <cassert>

using namespace tsd::vue;
using namespace std::complex_literals;

namespace tsd::telecom {


/*template<typename T>
  void clean(T &t)
{
    t = (t.abs() > 1e-8f).select(t, 0.0f);
}*/


static inline int decode_QPSK(const cfloat &x)
{
  if(x.real() >= 0)
  {
    if(x.imag() >= 0)
      return 0;
    else
      return 3;
  }
  else if(x.imag() >= 0)
    return 1;
  else
    return 2;
}

static inline int decode_BPSK(const cfloat &x)
{
  //return (x.real() >= 0) ? 0 : 1;
  return (x.real() >= 0) ? 1 : 0;
}

static ArrayXcf psk_constellation(int M)
{
  ArrayXcf x(M);

  if(M == 2)
  {
    x << -1, 1;
    return x;
  }

  float dec = (M == 4) ? π / 4 : 0;
  for(auto i = 0; i < M; i++)
    x(i) = std::polar(1.0, dec + (2 * π * i) / M);
  return x;
}

static ArrayXf ask_constellation(int M, float K1, float K2)
{
  return K1 + linspace(0, M-1, M) * (K2 / (M-1));
}

template<typename T>
  void clean(Vecteur<std::complex<T>> &t)
{
  auto tr = (t.real().abs() > 1e-7f).select(t.real(), 0.0f);
  auto ti = (t.imag().abs() > 1e-7f).select(t.imag(), 0.0f);
  t = tr + 1.0if * ti;
}


template<typename D>
 auto modulo(const Eigen::ArrayBase<D> &x, unsigned int M)
{
  return x - M * (x * (1.0f / M)).floor();
}



// Train binaire -> index
ArrayXi symmap_binaire(const BitStream &x, int k)
{
  int n = x.lon();

  if(k <= 0)
    echec("symmap_binaire : k <= 0.");

  int nsymbs = (n + k-1)/k;
  ArrayXi y(nsymbs);
  for(auto i = 0; i < nsymbs; i++)
  {
    int symb = 0;
    int e = 1;
    for(auto j = 0; (j < k) && (i*k+j < n); j++)
    {
      symb = symb + e * x[i*k+j];
      e   *= 2;
    }
    y(i) = symb;
  }
  return y;
}

// Index -> train binaire
void symdemap_binaire(BitStream &bs, const ArrayXi &x, int k)
{
  int n = x.rows();
  for(auto i = 0; i < n; i++)
    for(auto j = 0; j < k; j++)
      bs.push((((int) x(i)) >> j) & 1);
}

ArrayXcf FormeOnde::génère_symboles(const BitStream &bs)
{
  auto x = symmap_binaire(bs, infos.k);
  auto n = x.rows();
  ArrayXcf y(n);
  for(auto i = 0; i < n; i++)
  {
    cnt = i;
    y(i) = lis_symbole(x(i));
  }
  return y;
}

struct CtxGenImpl: FormeOnde::CtxGen
{
  FormeOnde *parent;
  sptr<FiltreGen<cfloat>> f;
  int ncoefs, osf;

  CtxGenImpl(FormeOnde *parent, int ncoefs, int osf)
  {
    this->parent = parent;
    f = parent->filtre.filtre_mise_en_forme(ncoefs, osf);
    this->osf     = osf;
    this->ncoefs  = parent->filtre.get_coefs(ncoefs, osf).rows();
  }
  void reset()
  {

  }
  ArrayXcf step(const BitStream &bs)
  {
    return ArrayXcf();
    /*auto &mconfig = config.format.modulation;

    auto osf = mconfig.fe / mconfig.fsymb;

    int   d_ech           = (int) ceil(mod->delais());
    int   nbits_par_symb  = mconfig.wf->k;
    float nbits_par_echan = ((float) nbits_par_symb) / osf;
    int   d_bit           = (int) ceil(d_ech * nbits_par_echan);

    // Ici, non, il faut concaténer deux modulateurs
    return mod->step(bs);*/
  }
};


sptr<FormeOnde::CtxGen> FormeOnde::get_contexte_tx(int ncoefs, int osf)
{
  return std::make_shared<CtxGenImpl>(this, ncoefs, osf);
}


// TODO : stocker dans le (un) contexte les différents objet de filtrage
ArrayXcf FormeOnde::génère_échantillons(const BitStream &bs, int ncoefs,
    int osf, float &retard)
{
  ArrayXcf symbs = génère_symboles(bs);
  auto f = filtre.filtre_mise_en_forme(ncoefs, osf);
  ncoefs = filtre.get_coefs(ncoefs, osf).rows();
  int nsymbs_flush = ncoefs;
  symbs = vconcat(symbs, ArrayXcf::Zero(nsymbs_flush));

  // Temps vers le milieu du premier bit transmis
  retard = tsd::filtrage::filtre_rif_ups_délais(ncoefs, osf);

  ArrayXcf y = f->step(symbs);

  if(infos.est_fsk)
  {
    // df = 0.5 * h * fsymb
    // df sur un symbole = 0.5 * h / osf
    // 2 π / osf <=> h = 2
    // => θ = h * 2 * pi / osf / 2 = h * pi / osf

    float Ω_max = (π_f * infos.index) / osf;
    // h = 2 -> Omega_max = 2 * pi / osf

    msg("Ω max = {} degrés.", rad2deg(Ω_max));

    // h = 2 fd / fsymb = excursion / fsymb
    // => omega_max = 2%pi*fd = %pi * h * fsymb
    // => Omega_max = omega_max / fs
    //              = %pi * h / ovs;
    // Conversion phase -> IQ


    // normalisation entre [-fmax,fmax]
    ArrayXf vfreqs = y.real() * (Ω_max / y.real().abs().maxCoeff());

    // TODO : échantillon précédent !
    y = polar(cumsum(vfreqs));
  }

  return y;
}


struct CtxSansMemoire: FormeOnde::Ctx
{
  const FormeOnde *fo;
  CtxSansMemoire(const FormeOnde *fo)
  {
    this->fo = fo;
  }
  // Pas de contexte
  void reset(){}


  std::tuple<int, cfloat> step(cfloat x)
  {
    int i = fo->symbole_plus_proche(x);
    return {i, fo->lis_symbole(i)};
  }

  /*void avance(int symbole){}
  cfloat lis_symbole(int i)
  {
    return fo->lis_symbole(i);
  }
  int symbole_plus_proche(cfloat x)
  {
    return fo->symbole_plus_proche(x);
  }*/
};

sptr<FormeOnde::Ctx> FormeOnde::get_ctx(int OSF) const
{
  return std::make_shared<CtxSansMemoire>(this);
}

// Optimisable
void  FormeOnde::decode_symboles(BitStream &bs, const ArrayXcf &x)
{
  auto n = x.rows();
  ArrayXi x2(n);
  for(auto i = 0; i < n; i++)
  {
    cnt = i;
    x2(i) = symbole_plus_proche(x(i));
  }
  symdemap_binaire(bs, x2, infos.k);
}


#if 0
ArrayXf symdemap(const ArrayXcf &xc, unsigned int k, const char enc = 'b')
{
// Maps a symbol sequence to a binary sequence with binary or gray decoding
//
// Calling Sequence
// b = symdemap(x,k[,enc])
//
// Parameters
// x: Input symbol vector (real or complex, depending on the encoding)
// k: Number of bits / symbol
// enc: Specify binary ('b', default), qam ('q'), or phase ('p') encoding
// b: Output binary sequence (a vector of 0 and 1)
//
// Description
// Maps a symbol sequence to a binary sequence (sequence of 0 or 1), using one of the following possible encodings:
// <itemizedlist>
//   <listitem><emphasis role="bold">Binary encoding (PAM / NRZ)</emphasis> / <varname>enc = 'b'</varname>: produces a vector of real symbols, uniformly distribued between -1 and 1 (if k = 2, only -1 and 1 are produced, e.g. NRZ encoding).</listitem>
// <listitem><emphasis role="bold">Phase encoding (PSK)</emphasis> / <varname>enc = 'p'</varname>: produces complex symbols of unit magnitude and uniformly distribued phase between 0 and <latex>$2\pi$</latex> radians (for k = 1 bit/symbol, this is BPSK, which symbols are identical to NRZ, for k = 2 bits/symbol, this is QPSK, ...).</listitem>
// <listitem><emphasis role="bold">Quadrature amplitude Encoding (QAM)</emphasis> / <varname>enc = 'q'</varname>: produces complex symbols where real and imaginary parts of the signal are independently modulated using each one half of the bitstream.</listitem>
// </itemizedlist>
//
// <refsection><title>Example</title></refsection>
// <programlisting>
//   b = symdemap([-1 1 1 -1], 1, 'b') // NRZ decoding
// </programlisting>
//
// See also
//  symmap
//
// Authors
//  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/sct">http://www.tsdconseil.fr/log/sct</ulink>

  ArrayXf x, b;

  if(xc.rows() == 0)
    return x;

  auto M = (unsigned int) std::pow(2, k);

  if(enc == 'g')
  {
    msg_erreur("symdemap: TODO (gray encoding)");
    return ArrayXf();
  }
  if((enc != 'b') && (enc != 'p') & (enc != 'q'))
  {
    msg_erreur("symdemap: encoding must be ''b'', ''g'' or ''pp''.");
    return ArrayXf();
  }

  if(enc == 'p')
  {
    x = xc.arg();
    //x.flo
    x = (x * (M / (2 * pi))).round();
    x += (float) M;
    x = modulo(x, M);
  }
  else if(enc == 'q')
  {
    if((k & 1) != 0)
    {
      msg_erreur("symdemap / QAM mode: k should be multiple of 2.");
      return ArrayXf();
    }
    ArrayXf b1 = symdemap(xc.real(), k/2, 'b');
    ArrayXf b2 = symdemap(xc.imag(), k/2, 'b');
    /*b.resize(2 * b1.rows());
    // TODO: à vectoriser
    for(auto i = 0; i < b1.rows(); i++)
    {
      b(2*i)   = b1(i);
      b(2*i+1) = b2(i);
    }*/
    auto n = b1.rows();
    b.resize(2 * n);
    for(auto i = 0; i < n; i++)
    {
      b(2*i)   = b1(i);
      b(2*i+1) = b2(i);
    }
    //b(Eigen::seq(0, 2 * n - 2, 2)) = b1;
    //b(Eigen::seq(1, 2 * n - 1, 2)) = b2;
    return b;
  }
  else
  {
    auto K1 = -1;
    auto K2 = 2 / (M-1);
    x = (xc.real() - K1) / K2;
  }

  // Map to binary sequence
  auto n = x.rows();
  b = ArrayXf::Zero(n*k);
  x = modulo((M+x.real().round()), M);
  for(auto i = 0; i < n; i++)
  {
    for(auto j = 0u; j < k; j++)
    {
      if(((int) x(i)) bitand (1 << j))
        b(i*k+j) = 1;
    }
  }

  return b;
}
ArrayXcf symmap(const ArrayXf &bs, int k, char mode = 'b')
{
  ArrayXf b = bs;
  auto n = b.rows();

  ArrayXf b2 = ArrayXf::Zero(n / k) ;
  ArrayXcf x = ArrayXcf::Zero(n / k) ;

  // Maps a binary sequence to a symbol sequence
  //
  // Calling Sequence
  // x = symmap(b,k[,enc])
  //
  // Parameters
  // b: Input binary sequence (must be a vector of 0.0 and 1.0)
  // k: Number of bits / symbol
  // enc: Specify binary ('b', default), phase ('p'), or qam ('q') encoding
  // x: Output symbol vector (real or complex)
  //
  // Description
  // Maps a binary sequence (sequence of 0 or 1) to a sequence of symbols, using one of the following possible encodings:
  // <itemizedlist>
  //   <listitem><emphasis role="bold">Binary encoding (PAM / NRZ)</emphasis> / <varname>enc = 'b'</varname>: produces a vector of real symbols, uniformly distribued between -1 and 1 (if k = 2, only -1 and 1 are produced, e.g. NRZ encoding).</listitem>
  // <listitem><emphasis role="bold">Phase encoding (PSK)</emphasis> / <varname>enc = 'p'</varname>: produces complex symbols of unit magnitude and uniformly distribued phase between 0 and <latex>$2\pi$</latex> radians (for k = 1 bit/symbol, this is BPSK, which symbols are identical to NRZ, for k = 2 bits/symbol, this is QPSK, ...).</listitem>
  // <listitem><emphasis role="bold">Quadrature amplitude Encoding (QAM)</emphasis> / <varname>enc = 'q'</varname>: produces complex symbols where real and imaginary parts of the signal are independently modulated using each one half of the bitstream.</listitem>
  // </itemizedlist>
  //
  // <refsection><title>Example</title></refsection>
  // <programlisting>
  //   x1 = symmap([0 1 0 0], 1, 'b') // NRZ encoding
  //   x2 = symmap([0 1 0 0], 2, 'p') // QPSK encoding
  // </programlisting>
  //
  // See also
  //  symdemap
  //
  // Authors
  //  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/sct">http://www.tsdconseil.fr/log/sct</ulink>


  if((k == 1) && (mode == 'b'))
      // Plus rapide !
      return 2 * (b - 0.5);

      //if(enc == 'g')
      //    error("symmap: TODO (gray encoding)");
      //end;
  if((mode != 'b') && (mode != 'p') && (mode != 'q'))
  {
    msg_erreur("symmap: encoding must be ''b'', ''g'' or ''pp''.");
    return ArrayXcf();
  }

  if(mode == 'q')
  {
    if(k & 1)
    {
      msg_erreur("symmap / QAM mode: k should be multiple of 2.");
      return ArrayXcf();
    }
    if(n & 1)
      n--;

    //auto b1 = bs(Eigen::seq(0, n-2, 2));
    //auto b2 = bs(Eigen::seq(1, n-1, 2));
    //auto b1 = sousech(bs, 2);
    //auto b2 = Eigen::Map<const ArrayXf,0,Eigen::InnerStride<>>(
    //    x.data()+1, bs.rows() / 2, Eigen::InnerStride<>(2));
        //sousech(bs, 2);

    ArrayXf b1(n/2), b2(n/2);
    for(auto i = 0; i < n/2; i++)
    {
      b1(i) = bs(2*i);
      b2(i) = bs(2*i+1);
    }


    auto s1 = symmap(b1,k/2,'b').real();
    auto s2 = symmap(b2,k/2,'b').real();

    x.resize(s1.rows());
    x.real() = s1;
    x.imag() = s2;

    return x;
  }

  msg("Symmap : k = {} bits / symbole.", k);

  for(auto i = 0; i < n/k; i++)
  {
    auto symb = 0;
    auto e = 1;
    for(auto j = 0; j < k; j++)
    {
      symb = symb + e * b(i*k+j);
      e = e * 2;
    }
    b2(i) = symb;
  }



  auto M = std::pow(2.0f, (float) k);

  if(mode == 'p')
  {
    ArrayXf phase = 2 * pi * b2 / M;

    //std::cout << "b2 = " << b2.transpose() << "\n";
    //std::cout << "phase = " << phase.transpose() << "\n";

    if(M == 4)
      phase += pi / 4;
    x = polar(phase);
    clean(x);
    if(M == 4)
      x *= sqrt(2.0f);
  }
  else
  {
    auto K1 = -1.0f;
    auto K2 = 2.0f / (M-1); // eg. K2 = 2 si k = 1 bit/symbole
    x = K1 + K2 * b2;
  }
  return x;
}
#endif

float FormeOnde::excursion() const
{
  return 1.0f;
}

int FormeOnde::symbole_plus_proche(const cfloat &point) const
{
  float bdist = 1e100;
  int res = 0;
  for(auto i = 0; i < infos.M; i++)
  {
    auto d2 = std::norm(point - lis_symbole(i));
    if(d2 < bdist)
    {
      bdist = d2;
      res = i;
    }
  }
  tsd_assert((res >= 0) && (res < (int) infos.M));
  return res;
}


struct FormeOndeASK: FormeOnde
{
  ArrayXcf symbs;
  float K1 = -1, K2 = 2;

  int symbole_plus_proche(const cfloat &x) const
  {
    return FormeOnde::symbole_plus_proche(x);
  }

  std::string desc_courte() const
  {
    return fmt::format("{}-ASK({},{})", infos.M, K1, K2);
  }

  std::string desc() const
  {
    return fmt::format("{}-ASK({},{})", infos.M, K1, K2);
  }


  FormeOndeASK(int M, float K1, float K2, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre  = filtre;
    infos.M       = M;
    this->K1      = K1;
    this->K2      = K2;
    infos.est_lineaire  = true;
    infos.est_ask       = true;
    infos.k             = std::log2(M);
    symbs         = ask_constellation(M, K1, K2);
  }

  ArrayXcf constellation() const
  {
    return symbs;
  }

  float ber(float EbN0)
  {
    // Equation 5.2.46 de Proakis
    float EbN0_lin = db2pow(EbN0);
    return ((infos.M-1.0f) / infos.M)
        * erfc(sqrt(3.0f * log2((float) infos.M) * EbN0_lin / (infos.M*infos.M - 1.0f)))/infos.k;
  }

  cfloat lis_symbole(unsigned int index) const
  {
    tsd_assert((int) index < infos.M);
    return symbs(index);
  }
};


struct FormeOndePSK: FormeOnde
{
  ArrayXcf symbs;

  int symbole_plus_proche(const cfloat &x) const
  {
    if(infos.M == 2)
      return x.real() > 0 ? 1 : 0;
      //return x.real() > 0 ? 0 : 1;
    return FormeOnde::symbole_plus_proche(x);
  }

  std::string desc_courte() const
  {
    if(infos.M == 2)
      return "BPSK";
    else if(infos.M == 4)
      return "QPSK";
    else
      return fmt::format("{}PSK", infos.M);
  }

  std::string desc() const
  {
    return desc_courte() + fmt::format(", {}", filtre);
  }


  FormeOndePSK(unsigned int M, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire = true;
    infos.M = M;
    infos.est_psk = true;
    infos.k = std::log2(M);
    symbs = psk_constellation(M);
  }

  ArrayXcf constellation() const
  {
    return symbs;
  }

  float ber(float EbN0)
  {
    float EbN0_lin = db2pow(EbN0);
    //eq 5.2.61 from proakis


    // Equation générale (valable sauf en BPSK)
    float ber = erfc(std::sqrt(infos.k*EbN0_lin)*sin(π/infos.M))/infos.k;

    // Si k = 1, M = 2 (BPSK) :
    // sin(pi/2) = 1
    // -> erfc(std::sqrt(ebnoD))

    // Si k = 2, M = 4 (QPSK) :
    // sin(pi/4) = sqrt(2)/2
    // erfc(std::sqrt(ebnoD) * sqrt(2) * sqrt(2) / 2) = erfc(std::sqrt(ebnoD)) / 2

    if(infos.M == 2)
      ber /= 2;
    return ber;
    //if(opt=='i')
      //ber=2*ber; // Differential modulation
    //end
  }

  cfloat lis_symbole(unsigned int index) const
  {
    tsd_assert((int) index < infos.M);
    return symbs(index);
  }
};


struct FormeOnde_π4QPSK: FormeOnde
{
  ArrayXcf symbs[2];

  int symbole_plus_proche(const cfloat &x) const
  {
    return FormeOnde::symbole_plus_proche(x);
  }

  std::string desc_courte() const
  {
    return "pi4-QPSK";
  }

  std::string desc() const
  {
    return fmt::format("pi4-QPSK, {}", filtre);
  }

  struct Ctxπ4QPSK: FormeOnde::Ctx
  {
    const FormeOnde_π4QPSK *parent;
    int cnt = 0;
    const cfloat rot = std::polar(1.0f, -π_f/4);

    Ctxπ4QPSK(const FormeOnde_π4QPSK *parent)
    {
      this->parent = parent;
    }

    void reset()
    {
      cnt = 0;
    }

    /** Index = -1 si pas d'échantillon à sortir */
    std::tuple<int, cfloat> step(cfloat x)
    {
      if(cnt++ & 1)
        x *= rot;
      int r = decode_QPSK(x);
      return {r, parent->symbs[(cnt+1) & 1](r)};
    }
  };

  sptr<Ctx> get_ctx(int OSF) const
  {
    return std::make_shared<Ctxπ4QPSK>(this);
  }


  FormeOnde_π4QPSK(const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire = true;
    infos.M = 4;
    infos.est_psk = true;
    infos.k = 2;
    symbs[0] = psk_constellation(4);
    symbs[1] = symbs[0] * std::polar(1.0f, π_f / 4);
  }

  ArrayXcf constellation() const
  {
    return vconcat(symbs[0], symbs[1]);
  }

  float ber(float EbN0)
  {
    // Ber idem QPSK
    float EbN0_lin = db2pow(EbN0);
    return erfc(std::sqrt(infos.k*EbN0_lin)*sin(π/infos.M))/infos.k;
  }

  cfloat lis_symbole(unsigned int index) const
  {
    tsd_assert((int) index < infos.M);
    return symbs[cnt & 1](index);
  }
};


struct FormeOndeQAM: FormeOnde
{
  ArrayXcf symbs;

  std::string desc_courte() const
  {
    return fmt::format("QAM{}", infos.M);
  }

  std::string desc() const
  {
    return fmt::format("QAM{}, {}", infos.M, filtre);
  }

  FormeOndeQAM(unsigned int M, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire = true;
    infos.M = M;
    infos.est_qam = true;
    infos.k = std::log2(M);
    //nom = fmt::format("QAM{}", M);

    auto M2 = (unsigned int) std::sqrt(M);

    if(M2 * M2 != M)
    {
      msg_erreur("FormeOnde QAM : M devrait être un carré (M = {}).", M);
      return;
    }

    symbs.resize(M);

    auto map = Eigen::Map<ArrayXXcf>(symbs.data(), M2, M2);
    for(auto x = 0u; x < M2; x++)
    {
      for(auto y = 0u; y < M2; y++)
      {
        float xf = (float) x, yf = (float) y;

        xf /= (M2 - 1);
        yf /= (M2 - 1);

        // xf, yf compris entre 0 1
        // 0 -> -1
        // 1 ->  1
        map(x, y) = cfloat(xf * 2 - 1, yf * 2 - 1);
      }
    }
  }

  ArrayXcf constellation() const
  {
    return symbs;
  }

  /*ArrayXcf gene_symboles(const ArrayXf &bs)
  {
    return symmap(bs, k, 'q');
  }
  ArrayXf demodulation(const ArrayXcf &x)
  {
    return symdemap(x, k, 'q');
  }*/

  cfloat lis_symbole(unsigned int index) const
  {
    tsd_assert((int) index < infos.M);
    return symbs(index);
  }
  float ber(float EbN0)
  {
    float ebnoD = std::pow(10.0f, EbN0/10);
    // eq 5.2.79 from proakis
    return (2.0f/infos.k)*(1-(1/sqrt(infos.M)))*erfc(sqrt(3.0f*infos.k*ebnoD/(2.0f*(infos.M-1))));
    //return (1-carre(1-tmp))/k;
  }
};






struct FormeOndeFSK: FormeOnde
{
  ArrayXf symbs;

  struct CtxFSK: FormeOnde::Ctx
  {
    int OSF; // Comment le configurer ?
    const FormeOndeFSK *fo;

    // Etat de phase initial du symbole en cours
    cfloat etat_phase = 0;
    // Nombre de pas depuis le début du symbole en cours
    int    cnt = 0;

    ArrayXcf symbs;

    // Tableau d'erreur pour chaque trajectoire de phase
    ArrayXf   err;
    ArrayXcf  cands;

    CtxFSK(const FormeOndeFSK *fo, int OSF)
    {
      this->OSF = OSF;
      this->fo = fo;

      int M = fo->infos.M;

      symbs.resize(M);
      if(M != 2)
        echec("TODO: CtxFSK : M != 2");

      // Si index = 2, il faut un OSF > 1
      symbs(0) = std::polar(1.0, - 2 * π * fo->infos.index / (2 * OSF));
      symbs(1) = std::polar(1.0, + 2 * π * fo->infos.index / (2 * OSF));

      err.resize(M);
      cands.resize(M);
      reset();
    }

    void reset()
    {
      etat_phase = 0;
      cnt        = 0;
      err.setZero();
      cands.setConstant(1.0f);
    }

    std::tuple<int, cfloat> step(cfloat x)
    {
      // Pour chaque symbole, on fait M hypothèses
      for(auto i = 0; i < (int) fo->infos.M; i++)
      {
        cands(i) = (cands(i) * symbs(i)) / std::abs(cands(i));
        err(i)   += std::norm(cands(i) - x); // std::norm ABS2
      }
      cnt++;
      if(cnt == OSF)
      {
        int i;
        cnt = 0;
        err.minCoeff(&i);
        err.setZero();
        cands.setConstant(cands(i));
        return {i, cands(i)};
      }
      else
      {
        return {-1, 0}; // Pas de symbole à sortir pour l'instant
      }
    }

    /*int step(cfloat obs)
    {

      return 0;
    }*/



    /*void avance(int i)
    {
      etat *= lis_symbole(i);
    }
    cfloat lis_symbole(int i) const
    {
      return etat_phase * std::pow(symbs(i), cnt);
      // return fo->lis_symbole(i);
    }*/
  };


  /*sptr<Ctx> get_ctx(int OSF) const
  {
    return std::make_shared<CtxFSK>(this, OSF);
  }*/

  std::string desc_courte() const
  {
    std::string nom = format("{}{}SK",
        filtre.type == SpecFiltreMiseEnForme::Type::GAUSSIEN ? "G" : "",
            infos.index == 0.5 ? "M" : "F");
    if(infos.M != 2)
      nom = fmt::format("{}{}", infos.M, nom);
    return nom;
  }

  std::string desc() const
  {
    return desc_courte() + fmt::format(", index={}", infos.index);
  }

  float excursion() const
  {
    return infos.index;
  }

  FormeOndeFSK(unsigned int M, float index, const SpecFiltreMiseEnForme &filtre)
  {
    infos.est_lineaire = false;
    infos.M = M;
    infos.index = index;
    this->filtre = filtre;
    infos.est_fsk = true;
    infos.k = std::log2(M);

    //if(filt == 'g')
//      filtre = SpecFiltreMiseEnForme::gaussien(BT, 31);
//    else
//      filtre = SpecFiltreMiseEnForme::nrz();

    symbs.resize(M);

    for(auto x = 0u; x < M; x++)
    {
      float xf = (float) x;

      xf /= (M - 1);

      // xf, yf compris entre 0 1
      // 0 -> -1
      // 1 ->  1
      symbs(x) = xf * 2 - 1;
    }
  }

  ArrayXcf constellation() const
  {
    return ArrayXcf();
  }
  float ber(float EbN0)
  {
    // Détecteur FSK cohérent et orthogonal (?)
    /*float F = std::pow(2.0f,k-1.0f)/(std::pow(2.0f,(float)k)-1);
    function [y]=fskcoherorth(q)
      A=sqrt(2*k*ebnoD(i)),
      B=q*q/2,
      C=0.5*erfc((-q-A)/sqrt(2)),
      y=(C^(M-1))/sqrt(2*(%pi))*exp(-B),
    endfunction
    Ps(i)=1-intg(-100,100,fskcoherorth);
    Pb(i)=F*Ps(i);*/
    // Détecteur FSK cohérent et non orthogonal
    float ebnoD = std::pow(10.0f, EbN0/10);
    cfloat Rho = 0.0f; // ?
    float A = ebnoD*(1-real(Rho));
    return 0.5 * erfc(sqrt(A/2));
  }
  // En supposant une discrimination polaire faite avant...
  cfloat lis_symbole(unsigned int index) const
  {
    //return 0.0f;
    tsd_assert((int) index < infos.M);
    return symbs(index);
  }
  float etat_phase = 0;
  cfloat lis_symbole(int index) const
  {
    //return 0.0f;
    tsd_assert(index < infos.M);
    return symbs(index);
  }
};



//  wf = wf_init('fsk'[,M=2,index=0.5,filt='n',BT=0.8])
//  wf = wf_init('gfsk'[,index=0.5])
//  wf = wf_init('gmsk')
//  wf = wf_init('msk')

sptr<FormeOnde> forme_onde_fsk(unsigned int M, float index, const SpecFiltreMiseEnForme &filtre)
{
  return std::make_shared<FormeOndeFSK>(M, index, filtre);
}


sptr<FormeOnde> forme_onde_psk(unsigned int M, const SpecFiltreMiseEnForme &filtre)
{
  return std::make_shared<FormeOndePSK>(M, filtre);
}

sptr<FormeOnde> forme_onde_ask(int M, float M1, float M2, const SpecFiltreMiseEnForme &filtre)
{
  return std::make_shared<FormeOndeASK>(M, M1, M2, filtre);
}

sptr<FormeOnde> forme_onde_qam(unsigned int M, const SpecFiltreMiseEnForme &filtre)
{
  return std::make_shared<FormeOndeQAM>(M, filtre);
}

sptr<FormeOnde> forme_onde_bpsk(const SpecFiltreMiseEnForme &filtre)
{
  return forme_onde_psk(2, filtre);
}

sptr<FormeOnde> forme_onde_qpsk(const SpecFiltreMiseEnForme &filtre)
{
  return forme_onde_psk(4, filtre);
}





sptr<FormeOnde> forme_onde_π4_qpsk(const SpecFiltreMiseEnForme &filtre)
{
  return std::make_shared<FormeOnde_π4QPSK>(filtre);
}




ArrayXf SpecFiltreMiseEnForme::get_coefs(int ncoefs, int osf) const
{
  tsd_assert_msg(osf > 0, "SpecFiltreMiseEnForme : osf invalide ({})", osf);

  ArrayXf coefs = ArrayXf::Zero(osf);
  coefs(0) = 1;

  if(osf == 1)
  {
    return coefs; // Filtre neutre
    //echec("Spec filtre mise en forme::get_coefs() : osf = 1.");
  }

  if(ncoefs == 0)
  {
    ncoefs = 5 * osf + 1;
    if((ncoefs & 1) == 0)
      ncoefs++;
  }

  // Construction du filtre
  if(type == NRZ)
  {
    return ArrayXf::Ones(osf) / osf; // Ntaps = osf = 4
  }
  else if(type == AUCUN)
  {
    return ArrayXf::Ones(1);
  }
  else if(type == GAUSSIEN)
  {
    tsd_assert(ncoefs > 0);
    //msg("gaussien : bt = {}, osf = {}, ncoefs = {}", BT, osf, ntaps);
    // Porte de largeur osf + filtre gaussien
    coefs = tsd::filtrage::design_rif_gaussien_telecom(ncoefs, BT, osf);
  }
  else if(type == RCS)
  {
    // Ntaps = 3 * osf = 12
    tsd_assert(ncoefs > 0);
    //msg("rcs : roff = {}, osf = {}, ncoefs = {}", α, osf, ncoefs);
    coefs = tsd::filtrage::design_rif_rcs1(ncoefs, β, osf);
  }
  else
    echec("Type de filtre de mise en forme inconnu ({})", (int) type);

  /*{
    Figure fg("Filtre de mise en forme");
    fg.plot(coefs, "b|r");
    stdo << fg;
  }*/

  //msg("Filtre de mise en forme : ncoefs = {} ({} prog, osf = {}).", coefs.rows(), ncoefs, osf);

  tsd_assert(coefs.rows() > 0);
  if(coefs.hasNaN())
    echec("{} : Valeurs NaN dans les coefficients.", *this);

  return coefs;
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_adapté(int ntaps, int osf) const
{
  if(osf == 1)
    return tsd::filtrage::filtre_id<cfloat>();
  ArrayXf h = get_coefs(ntaps, osf);
  return tsd::filtrage::filtre_rif<float, cfloat>(h);
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_adapté_décimation(int ntaps, int osf) const
{
  if(osf == 1)
    return tsd::filtrage::filtre_id<cfloat>();
  ArrayXf h = get_coefs(ntaps, osf);
  return tsd::filtrage::filtre_rif_decim<float, cfloat>(h, osf);
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_mise_en_forme(int ntaps, int R) const
{
  if(R == 1)
    return tsd::filtrage::filtre_id<cfloat>();
  ArrayXf h = get_coefs(ntaps, R);
  // Normalisation en énergie

  float en = h.square().sum();

  // msg("filtre_mise_en_forme: en={}, R={}", en, R);

  // Facteur 1/R à cause de normalisation par R faite dans filtre_rif_ups (d'ailleurs à vérifier)
  h *= (std::sqrt(R) / std::sqrt(en)) / R;

  return tsd::filtrage::filtre_rif_ups<float, cfloat>(h, R);
}



ArrayXf FormeOnde::ber(const ArrayXf &EbN0_dB)
{
  int n = EbN0_dB.rows();
  ArrayXf y(n);
  for(auto i = 0; i < n; i++)
    y(i) = ber(EbN0_dB(i));
  return y;
}

std::ostream& operator<<(std::ostream &ss, const SpecFiltreMiseEnForme &t)
{
  //ss << "SpecFiltreMiseEnForme : ";
  if(t.type == SpecFiltreMiseEnForme::NRZ)
    ss << "NRZ";
  else if(t.type == SpecFiltreMiseEnForme::AUCUN)
    ss << "AUCUN";
  else if(t.type == SpecFiltreMiseEnForme::RCS)
    ss << fmt::format("RCS-dep={}", (int) (t.β * 100.0f));
  else if(t.type == SpecFiltreMiseEnForme::GAUSSIEN)
    ss << fmt::format("Gaussien-BT={}", t.BT);
  else
    ss << "?";
  return ss;
}

std::ostream& operator<<(std::ostream &ss, const FormeOnde &t)
{
  ss << t.desc();
  return ss;
}

SpecFiltreMiseEnForme SpecFiltreMiseEnForme::aucun()
{
  return SpecFiltreMiseEnForme{.type = Type::AUCUN};
}

SpecFiltreMiseEnForme SpecFiltreMiseEnForme::gaussien(float BT)
{
  return SpecFiltreMiseEnForme{.type = Type::GAUSSIEN, .BT = BT};
}
SpecFiltreMiseEnForme SpecFiltreMiseEnForme::nrz()
{
  return SpecFiltreMiseEnForme{.type = Type::NRZ};
}
SpecFiltreMiseEnForme SpecFiltreMiseEnForme::rcs(float β)
{
  return SpecFiltreMiseEnForme{.type = Type::RCS, .β = β};
}

// Ajoute si besoin des zéros avant et après
static ArrayXf sympad(const ArrayXf &v, int nc)
{
  tsd_assert(v.rows() <= nc);
  auto r = v;
  if(v.rows() < nc)
  {
    int d = (nc - v.rows()) / 2;
    r = vconcat(ArrayXf::Zero(d), vconcat(v, ArrayXf::Zero(nc - d - v.rows())));
  }
  tsd_assert(r.rows() == nc);
  return r;
}

SpecFiltreMiseEnForme::Analyse SpecFiltreMiseEnForme::analyse(int ncoefs, int osf) const
{
  SpecFiltreMiseEnForme::Analyse res;
  Figures f;

  ArrayXf h = get_coefs(ncoefs, osf);
  h = sympad(h, ncoefs);

  int nc2 = 8*1024-1, osf2 = std::max(20, 2 * osf);
  ArrayXf ht = get_coefs(nc2, osf2);
  ht = sympad(ht, nc2);

  // Supprime les échantillons nuls, pour clarifier l'affichage
  int i = 0;
  while((i < nc2/2) && (std::abs(ht(i)) < 1e-4) && ((nc2 - 2 * i) / ((1.0 * osf2) / osf) > ncoefs))// && (<std::abs(ht(i)) < h(0))<)
    i++;
  nc2 -= 2 * i;
  ht = get_coefs(nc2, osf2);
  ht = sympad(ht, nc2);


  int m = ncoefs/2, m2 = nc2/2;
  ArrayXf t  = linspace(-(1.0f*m)/osf,(1.0f*m)/osf,ncoefs);
  ArrayXf t2 = linspace(-(1.0f*m2)/osf2,(1.0f*m2)/osf2,nc2);

  ArrayXf hte = ht * ((osf2 * 1.0f) / osf);

  f.subplot();
  f.gcf().plot(t, h,  "|bo", "Réel");
  f.gcf().plot(t2, hte, "-g", "Théorique");
  f.gcf().titres(fmt::format("{} - {} coefs, OSF={}", *this, ncoefs, osf), "Symboles");
  f.subplot();
  auto [fr, mag] = tsd::filtrage::frmag(h);
  auto [fr2, mag2] = tsd::filtrage::frmag(ht);
  f.gcf().plot(fr, mag,  "-b", "Réel");

  float r = ((osf2 * 1.0f) / osf);
  int ns = mag2.rows() / r;
  ArrayXf fr3(ns);
  for(auto i = 0; i < ns; i++)
    fr3(i) = fr2((int)floor(i * r));
  f.gcf().plot(fr3, mag2.head(ns), "-g", "Théorique");
  f.gcf().titre("Réponse fréquentielle");


  f.subplot();
  f.gcf().plot(fr, 20 * log10(mag + 1e-10),  "-b", "Réel");
  f.gcf().plot(fr3, 20 * log10(mag2.head(ns) + 1e-10), "-g", "Théorique");
  f.gcf().titre("Réponse fréquentielle (log)");


  res.fig = f;


  return res;
}

}

