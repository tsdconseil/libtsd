#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/telecom.hpp"
#include "tsd/telecom/modulations.hpp"
#include "tsd/figure.hpp"
#include <Eigen/QR>

using namespace tsd::vue;

namespace tsd::telecom {


#if 0
//template<typename T, typename Tc>
//Eigen::Ref<Vecteur<T>> convol
//  (const Eigen::Ref<const Vecteur<Tc>> h, const Eigen::Ref<const Vecteur<T>> x)

template<typename T, typename Tc>
  auto convol(const Eigen::ArrayBase<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  auto f = tsd::filtrage::filtre_rif<T, Tc>(h);
  return f->step(x);
}
#endif


auto convol(IArrayXf h, IArrayXf x)
{
  auto f = tsd::filtrage::filtre_rif<float, float>(h);
  return f->step(x);
}

auto convol(IArrayXf h, IArrayXcf x)
{
  auto f = tsd::filtrage::filtre_rif<float, cfloat>(h);
  return f->step(x);
}



struct Egaliseur: FiltreGen<cfloat>
{
  sptr<FormeOnde> wf;
  float K = 1, gain = 0.01, g2 = 1;
  int N1 = 11, N2 = 11, cnt = 0;

  ArrayXf h, h_avant, h_arriere;
  ArrayXcf wnd, wnd_deci;

  bool init_ok = false;

  enum Structure
  {
    DDE, DFE
  } structure = DDE;

  enum Errf
  {
    SLICER, CMA
  } errf = SLICER;

  /*std::string lis_structure() const
  {
    return structure ==
  }*/

  Egaliseur(sptr<FormeOnde> wf, const std::string &s_structure, const std::string &s_errf,
      float osf, float gain, int N1, int N2)
  {
    init_ok = false;
    this->K = osf;
    this->gain = gain;
    this->N1 = N1;
    this->N2 = N2;
    // Creation of an equalizer object (LMS, CMA, ...)
    //
    // Calling Sequence
    // eq = equalizer_init(wf); // Default equalizer (DDE, K = 1)
    // eq = equalizer_init(wf,K,errf,struct,gain,N1[,N2]]);
    //
    // Parameters
    // wf: Waveform specification (see <link linkend="wf_init">wf_init</link>)
    // K: Input over-sampling ratio (number of samples per symbol)
    // errf: Error function (either 'slicer' or 'cma')
    // struct: Equalizer structure (either 'dde' or 'dfe')
    // gain: Update gain (determine the convergence speed)
    // N1: Number of taps of the direct FIR equalization filter
    // N2: Number of taps of the feedback FIR equalization filter (only for DFE equalization)
    // eq: Resulting equalizer object (can be used with <link linkend="equalizer_process">equalizer_process</link>)
    //
    //
    // Description
    // Initialize an equalizer object, ever symbol sampled (K = 1), or
    // fractionnaly sampled (K > 1). Supported equalization structures are the
    // following:
    // <itemizedlist>
    //   <listitem><para><emphasis role="bold">Decision Directed Equalization (DDE)</emphasis> <varname>(strucuture = 'dde'):</varname> A direct FIR filter is tuned (every symbol period) to minimize the output error.</para></listitem>
    //   <listitem><para><emphasis role="bold">Decision Feedback Equalization (DFE)</emphasis> <varname>(strucuture = 'dfe'):</varname> Both a direct FIR filter (working at <latex>$K\cdot f_{symb}$</latex>) and a feedback filter (working at symbol frequency) on the decision (slicer) outputs are used.</para></listitem>
    // </itemizedlist>
    //
    // Possible cost function are the following:
    // <itemizedlist>
    //   <listitem><para><emphasis role="bold">Slicer</emphasis> <varname>(errf = 'slicer'):</varname> <latex>$E=(d-y)^2$</latex>. Note: with this error function, the algorithm is better known as <emphasis role="bold">LMS (Least Mean Square).</emphasis></para></listitem>.
    //   <listitem><para><emphasis role="bold">Constant modulus algorithm (CMA)</emphasis> <varname>(errf = 'cma'):</varname> <latex>$E=\left(R-|y|^2\right)^2$</latex></para></listitem>
    // </itemizedlist>
    // For a complete example, see the <link linkend="equalizer_process">equalizer_process</link> function.
    //
    // See also
    //  equalizer_process
    //  wf_init
    //
    // Authors
    //  J. Arzi, full documentation available on <ulink url="http://www.tsdconseil.fr/log/sct">http://www.tsdconseil.fr/log/sct</ulink>

        this->wf  = wf;


        if((s_errf != "cma") && (s_errf != "slicer"))
        {
          msg_erreur("init égaliseur: supported error function = \"slicer\" or \"cma\" (requested : \"{}\").", s_errf);
          return;
        }

        if((s_structure != "dde") && (s_structure != "dfe"))
        {
          msg_erreur("init égaliseur: supported structure = \"dde\" or \"dfe\" (requested : \"{}\").", s_structure);
          return;
        }

        structure = s_structure == "dde" ? DDE : DFE;
        errf = s_errf == "slicer" ? SLICER : CMA;


        if(structure == DFE)
        {
          h_avant    = vconcat(ArrayXf::Zero(N1-1), ArrayXf::Ones(1));//[zeros(N1-1,1) ; 1];
          h_arriere  = ArrayXf::Zero(N2);
          wnd        = ArrayXf::Zero(N1);
          wnd_deci   = ArrayXf::Zero(N2);
          cnt        = 0;
        }
        else if(structure == DDE)
        {
          h     = vconcat(ArrayXf::Zero(N1-1), ArrayXf::Ones(1));
          wnd   = ArrayXf::Zero(N1);
        }
        else
        {
          msg_erreur("equalizer_init: supported structures = 'dde' or 'dfe'.");
          return;
        }

        //if((eq.errf == 'cma') & (eq.structure ~= 'dde'))
          //  error("equalizer_init: CMA (Constant Modulus Algorithm) works only with DDE (Decision Directed Eq.).");
        //end;

        if(errf == CMA)
        {
          g2 = 1; // For PSK, FSK
          if(wf->est_qam)
            msg_erreur("TODO: scale factor for QAM and CMA.");
        }

        init_ok = true;
  }

  void step(IArrayXcf x, ArrayXcf &y)
  {
    if(init_ok)
    {
      ArrayXf err;
      step_interne(x, y, err);
    }
  }
  void step_interne(IArrayXcf x, ArrayXcf &y, ArrayXf &err)
  {
    // Equalization of an input signal
    //
    // Calling Sequence
    // [eq,y] = equalizer_process(eq,x);
    // [eq,y,err] = equalizer_process(eq,x);
    //
    // Parameters
    //  x: input signal to be equalized
    //  y: output signal
    //  err: Error for each symbol
    //
    // Description
    //  Perform CMA or LMS equalization, using DDE (Decision directed) or DFE (Decision feedback) structure.
    //
    // <refsection><title>Example</title></refsection>
    // <programlisting>
    //   fs = 4; fi = 0; fsymb = 1;
    //   wf = wf_init('qpsk');
    //   eq = equalizer_init(wf, fs/fsymb, 'cma', 'dfe', 0.01,21);
    //   h = [1 0.1 -.1 0.15 0.05];
    //   mod = mod_init(wf, fs,fi,fsymb);
    //  [mod,x] = mod_process(mod,prbs(2000));
    //   y = convol(h,x);
    //   y = awgn(y, 0.02, 'c');
    //   [eq,z] = equalizer_process(eq, y);
    // </programlisting>
    // <mediaobject><imageobject><imagedata fileref="ex-dfe-cma-const.png" format="PNG"/></imageobject><caption><para>Example: DFE / CMA equalization</para></caption></mediaobject>
    //
    // See also
    //  equalizer_init
    //
    // Authors
    //  J. Arzi, full documentation available on <ulink url="http://www.tsdconseil.fr/log/sct">http://www.tsdconseil.fr/log/sct</ulink>

    auto n = x.rows();
    y = ArrayXf::Zero(n);
    err = ArrayXf::Zero(n);

    if(structure == DDE)
    {
      for(auto i = 0; i < n; i++)
      {
        // wnd = [wnd(2:N1) ; x(i)];
        wnd.head(N1-1) = wnd.tail(N1-1).eval();
        wnd(N1-1) = x(i);

        auto egalise  = (h * wnd).sum();
        auto d = wf->lis_symbole(wf->symbole_plus_proche(egalise));
        float e;

        if(errf == CMA)
        {
          e = g2 - std::abs(egalise) * std::abs(egalise);
          h = h + gain * e * (wnd * std::conj(egalise)).real();
        }
        else
        {
          auto ec = d - egalise;
          e = std::abs(ec);
          h = h + gain * (std::conj(ec) * wnd).real();
        }
        y(i)   = egalise;
        err(i) = e;
      }
    }
    else if(structure == DFE)
    {
      auto j = 1;
      for(auto i = 0; i < n; i++)
      {
        wnd.head(N1-1) = wnd.tail(N1-1).eval();
        wnd(N1-1) = x(i);
        cnt++;
        if(cnt != K)
          continue;

        cnt = 0;
        auto retroaction = (h_arriere * wnd_deci).sum();
        auto egalise  = (h_avant * wnd).sum() + retroaction;

        auto d = wf->lis_symbole(wf->symbole_plus_proche(egalise));
        float e;
        if(errf == CMA)
        {
          e = g2 - std::abs(egalise * egalise);
          h_avant   += gain * e * (wnd * std::conj(egalise)).real();
          h_arriere += gain * e * (wnd_deci * std::conj(egalise)).real();
        }
        else
        {
          auto ec = d - egalise;
          e = std::abs(ec);
          h_avant   += gain * (std::conj(ec) * wnd).real();
          h_arriere += gain * (std::conj(ec) * wnd_deci).real();
        }

        err(j) = e;
        y(j++) = egalise;

        wnd_deci.head(N2-1) = wnd_deci.tail(N2-1).eval();
        wnd_deci(N2-1) = d;

      } // for i

      y.resize(j);
      err.resize(j);
    } // structure DFE
  }

};


sptr<FiltreGen<cfloat>> égaliseur_création(sptr<FormeOnde> wf, const std::string &structure, const std::string &errf,
    float osf, float gain, int N1, int N2)
{
  return std::make_shared<Egaliseur>(wf, structure, errf, osf, gain, N1, N2);
}


ArrayXf égaliseur_zfe(IArrayXf h, int n)
{
  int m = h.rows();

  if(n < m)
  {
    msg_erreur("egaliseur_zfe(h[{}],n={}): n doit être supérieur ou égal à la longueur de la réponse du canal.", m, n);
    return ArrayXf();
  }

  // (1) Calcule un délais
  // Suppose que h et g sont centrées, delai = m/2 + n/2
  auto d = (int) round(m/2 + n/2);

  // (2) Construit la matrice de filtrage : n+m-1 lignes, n colonnes
  Eigen::MatrixXf F = Eigen::MatrixXf::Zero(n+m-1,n);

  // Ligne 0 à m-1
  for(auto i = 0; i < m; i++)
    F.block(i,0,1,i+1) = h.segment(0, i+1).reverse().transpose();//h(i+1:-1:1)';

  // Ligne m à n-1
  for(auto i = m; i <= n-1; i++)
  {
    // i = m : col = [1...1+m-1]
    // i = n-1: col = [n-m...n-1]
    F.block(i,i+1-m,1,m) = h.reverse().transpose();//h($:-1:1)';
  }

  // Ligne n à n+m-2
  for(auto i = n; i < n+m-1; i++)
    for(auto j = 1; j <= m-1-(i-n); j++)
      F(i,n-j) = h(j+(i-n));

  Eigen::VectorXf δ = Eigen::VectorXf::Zero(n+m-1);
  δ(d, 0) = 1.0f;

  //g = F \ δ;
  ArrayXf g = F.colPivHouseholderQr().solve(δ).array();

  return g;
}






}


