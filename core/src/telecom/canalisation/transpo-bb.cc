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

