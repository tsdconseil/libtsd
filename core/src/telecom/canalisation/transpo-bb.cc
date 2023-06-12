#include "tsd/filtrage.hpp"
#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"

using namespace std;
using namespace tsd;
using namespace tsd::filtrage;

namespace tsd::telecom {

// A FAIRE : intégration d'une décimation polyphase.
//

template<typename Te>
struct TranspoBB: Filtre<Te, cfloat, TranspoBBConfig>
{
  float fc = 0;
  char mode = 'r';
  sptr<SourceGen<cfloat>> ol;
  sptr<FiltreGen<cfloat>> filtre_image;

  void configure_impl(const TranspoBBConfig &config)
  {
    si(mode == 'r')
    {
      // Image filtering
      si(config.fi < 0.25)
        fc = config.fi;
      sinon
        fc = 0.5 - config.fi;

      soit ntaps = 63;
      soit h = design_rif_fen(ntaps, "lp", fc);
      filtre_image = tsd::filtrage::filtre_rif<float,cfloat>(h);
    }
    ol = source_ohc(-config.fi);
  }

  TranspoBB(const TranspoBBConfig &config, char mode)
  {
    this->mode    = mode;
    Configurable<TranspoBBConfig>::configure(config);
  }

  void step(const Vecteur<Te> &x, Veccf &y)
  {
    y = x * ol->step(x.rows());

    si constexpr(!est_complexe<Te>())
      y = filtre_image->step(y);
  }
};

template<typename T>
  sptr<Filtre<T,cfloat,TranspoBBConfig>> transpo_bb(const TranspoBBConfig &config)
  {
    si constexpr(est_complexe<T>())
      retourne make_shared<TranspoBB<cfloat>>(config, 'c');
    sinon
      retourne make_shared<TranspoBB<float>>(config, 'r');
  }

soit tbb1 = transpo_bb<float>;
soit tbb2 = transpo_bb<cfloat>;


} // namespace

