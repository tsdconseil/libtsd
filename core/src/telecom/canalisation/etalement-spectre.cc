#include "tsd/telecom.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"


namespace tsd::telecom {






struct FHSSModulation: Filtre<cfloat,cfloat,FHSSConfig>
{
  sptr<FiltreGen<cfloat>> ra;
  int index = 0, nfreqs = 0;
  std::vector<sptr<SourceGen<cfloat>>> ol;

  FHSSModulation(const FHSSConfig &config)
  {
    configure(config);
  }

  int configure_impl(const FHSSConfig &c)
  {
    auto &config = Configurable<FHSSConfig>::config;
    config = c;
    ra = tsd::filtrage::filtre_reechan<cfloat>((float) config.osf_out / config.osf_in);

    nfreqs = config.seq.rows();
    // Pré-calcul des n fréquences
    /*rotations.resize(nfreqs);
    for(auto i = 0; i < nfreqs; i++)
      rotations[i] = std::polar((2 * π * i) / nfreqs);*/

    msg("configure FHSS : nfreqs = {}, durée slot = {}, osf in = {}, osd_out = {}",
        nfreqs, c.duree_slot, config.osf_in, config.osf_out);

    ol.resize(nfreqs);
    for(auto i = 0; i < nfreqs; i++)
      ol[i] = source_ohc(((float) i) / nfreqs);

    index = 0;

    return 0;
  }
  void step(IArrayXcf x, ArrayXcf &y)
  {
    msg("Calcul FHSS...");
    auto &config = Configurable<FHSSConfig>::config;
    y = ra->step(x);
    msg("x.rows = {}, y.rows = {}.", x.rows(), y.rows());
    // Chaque slot est multiplié par une exponentielle différente
    for(auto i = 0; i + config.duree_slot <= y.rows(); i += config.duree_slot)
    {
      //ol.configure(rotations[index++ % nfreqs]);
      ArrayXcf osc = ol[i]->step(config.duree_slot);
      y.segment(i, config.duree_slot) *= osc;
    }
    // TODO : gestion si paquet pas multiple de ...
  }
};

sptr<Filtre<cfloat,cfloat,FHSSConfig>> fhss_modulation(const FHSSConfig &config)
{
  return std::make_shared<FHSSModulation>(config);
}

void fhss_exemple()
{
  //////// MODULATION
  /// -> Trivial

  //////// DEMODULATION
  /// Comment se synchroniser ?
  // Exemple avec méthode de Harris ?
  // on a disons une matrice NxM
  // On peut chercher en parallèle ?

  // Lien entre séquence PRBS (binaire) et choix de la fréquence (nombre entier) ?
  // Déterminisme entre fréquence n et fréquence n+1, ou bien état ?

  // Détecteur RAKE ??? --> récepteur multi-trajet

  // DS-CDMA
  // MC-CDMA
  // FH-CDMA

  // "Slow frequency hopping" : Thop > Tsymb
  // "Fast frequency hopping" : Thop < Tsymb (En pratique, Tsymb = k*Thop)
}



struct DSSSModulation: Filtre<cfloat,cfloat,DSSSConfig>
{
  int osf = 0;
  int nbits = 0, index = 0;

  ArrayXf seq;
  sptr<FiltreGen<cfloat>> ra;


  DSSSModulation(const DSSSConfig &c)
  {
    configure(c);
  }

  int configure_impl(const DSSSConfig &c)
  {
    auto &config = Configurable<DSSSConfig>::config;
    index   = 0;
    config  = c;
    nbits   = c.chips.rows();
    osf     = nbits / c.osf_in;
    if((nbits % c.osf_in) != 0)
    {
      msg_erreur("DSSS : nbits doit être un multiple de osf in.");
      return -1;
    }
    ra = tsd::filtrage::filtre_reechan<cfloat>(((float) nbits) / c.osf_in);

    return 0;
  }

  void step(const Eigen::Ref<const ArrayXcf> x, ArrayXcf &y)
  {
    auto &config = Configurable<DSSSConfig>::config;
    // Pour chaque symbole d'entrée, consomme un bit de la séquence
    // Mais, un symbole d'entrée = plusieurs échantillons
    tsd_assert((x.rows() % config.osf_in) == 0);
    //auto nsymb = x.rows() / config.osf_in;

    y = ra->step(x);

    // Maintenant, y est échantillonné à la fréquence bit (chips)
    for(auto j = 0; j < y.rows(); j++)
      y(j) *= config.chips(index++ % nbits);


    // Exemple avec OSF IN = 4
    // nbits = 4
    //      ....++++....++++
    // nbits = 8
    //      ....++++....++++
    //      ........++++++++........++++++++
  }
};


sptr<Filtre<cfloat,cfloat,DSSSConfig>> dsss_modulation(const DSSSConfig &config)
{
  return std::make_shared<DSSSModulation>(config);
}








}


using namespace tsd;
using namespace tsd::telecom;


namespace tsd{
  template<typename T>
    ArrayXf range(const Vecteur<T> &x)
  {
    return linspace(0, x.rows() - 1, x.rows());
  }
}









