#include "tsd/filtrage.hpp"
#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"

using namespace std;
using namespace tsd;
using namespace tsd::filtrage;

namespace tsd::telecom {

  // Initialization of a baseband down-converter.
  // Calling Sequence
  // dn = downconvert_init(nu)
  // dn = downconvert_init(nu, opt)
  //
  // Parameters
  // nu: Normalized RF or IF frequency (between -0.5 and 0.5)
  // opt: 'r' (default) for real input signal, 'c' for complex input signal
  //
  // Description
  // If the input signal is complex, no image filtering is down.
  // If the input signal is real, image is removed with a low-pass filter
  // with cut-off frequency set to <latex>$\nu$</latex> if <latex>$\nu < \frac{1}{4}$</latex>, and <latex>$\frac{1}{2}-\nu$</latex> otherwise
  // (in either case, half-way between baseband signal and image signal).

  // Down-conversion to baseband.
  //
  // Calling Sequence
  // [dn,y] = downconvert_process(dn,x)
  //
  // Parameters
  // x: input signal (real or complex)
  // dn: downconversion object (created with <link linkend="downconvert_init">downconvert_init</link>)
  // y: output baseband signal (complex)
  //
  // Description
  // If the input signal is complex, no image filtering is down.
  // If the input signal is real, image is removed with a low-pass filter
  // with cut-off frequency set to <latex>$\nu$</latex> if <latex>$\nu < \frac{1}{4}$</latex>, and <latex>$\frac{1}{2}-\nu$</latex>  otherwise,
  // where <latex>$\nu$</latex> is the normalized intermediate frequency
  // (in either case, half-way between baseband signal and image signal).
  //
  // <refsection><title>Example 1: downconversion of a BPSK complex signal (without image filtering)</title></refsection>
  // <programlisting>
  // // Load RF test signal (BPSK, I/Q)
  // [x,fs,fi,fsymb] = sct_test_signal('c');
  // ;
  // // Down-conversion to baseband
  // dn = downconvert_init(fi/fs,'c');
  // [dn,y] = downconvert_process(dn,x);
  // </programlisting>
  // Scatter plot view:
  // <imageobject><imagedata fileref="ex_downconvert1_const.png" format="PNG"/></imageobject>
  // Frequency view:
  // <imageobject><imagedata fileref="ex_downconvert1_psd.png" format="PNG"/></imageobject>
  //
  // <refsection><title>Example 2: downconversion of a BPSK real signal (with image filtering)</title></refsection>
  // <programlisting>
  // // Load RF test signal (BPSK, I/Q)
  // [x,fs,fi,fsymb] = sct_test_signal();
  // ;
  // // Down-conversion to baseband
  // dn = downconvert_init(fi/fs,'r');
  // [dn,y] = downconvert_process(dn,x);
  // </programlisting>
  // Scatter plot view:
  // <imageobject><imagedata fileref="ex_downconvert2_const.png" format="PNG"/></imageobject>
  // Frequency view:
  // <imageobject><imagedata fileref="ex_downconvert2_psd.png" format="PNG"/></imageobject>

  // A FAIRE : intégration d'une décimation polyphase.
  //

template<typename Te>
struct TranspoBB: Filtre<Te, cfloat, TranspoBBConfig>
{
  float fc = 0;
  char mode = 'r';
  sptr<SourceGen<cfloat>> ol;
  sptr<FiltreGen<cfloat>> filtre_image;
  //TranspoBBConfig config;

  int configure_impl(const TranspoBBConfig &config)
  {
    if(mode == 'r')
    {
      // Image filtering
      if(config.fi < 0.25)
          fc = config.fi;
      else
          fc = 0.5 - config.fi;

      auto ntaps = 63;
      ArrayXf h = design_rif_fen(ntaps, "lp", fc);
      filtre_image = tsd::filtrage::filtre_rif<float,cfloat>(h);
    }
    ol = source_ohc(-config.fi);
    return 0;
  }

  TranspoBB(const TranspoBBConfig &config, char mode)
  {
    this->mode    = mode;
    Configurable<TranspoBBConfig>::configure(config);
  }

  void step(const Eigen::Ref<const Vecteur<Te>> x, ArrayXcf &y)
  {
    y = x * ol->step(x.rows());

    if constexpr(!est_complexe<Te>())
      y = filtre_image->step(y);
  }
};

template<typename T>
  sptr<Filtre<T,cfloat,TranspoBBConfig>> transpo_bb(const TranspoBBConfig &config)
  {
    if constexpr(est_complexe<T>())
      return make_shared<TranspoBB<cfloat>>(config, 'c');
    else
      return make_shared<TranspoBB<float>>(config, 'r');
  }

auto tbb1 = transpo_bb<float>;
auto tbb2 = transpo_bb<cfloat>;


} // namespace

