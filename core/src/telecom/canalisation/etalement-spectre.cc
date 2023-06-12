#include "tsd/telecom.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"


namespace tsd::telecom {






struct FHSSModulation: Filtre<cfloat,cfloat,FHSSConfig>
{
  sptr<FiltreGen<cfloat>> ra;
  entier index = 0, nfreqs = 0;
  vector<sptr<SourceGen<cfloat>>> ol;

  FHSSModulation(const FHSSConfig &config)
  {
    configure(config);
  }

  void configure_impl(const FHSSConfig &c)
  {
    soit &config = Configurable<FHSSConfig>::config;
    config = c;
    ra = tsd::filtrage::filtre_reechan<cfloat>((float) config.osf_out / config.osf_in);

    nfreqs = config.seq.rows();
    // Pré-calcul des n fréquences
    /*rotations.resize(nfreqs);
    pour(auto i = 0; i < nfreqs; i++)
      rotations[i] = std::polar((2 * π * i) / nfreqs);*/

    msg("configure FHSS : nfreqs = {}, durée slot = {}, osf in = {}, osd_out = {}",
        nfreqs, c.duree_slot, config.osf_in, config.osf_out);

    ol.resize(nfreqs);
    pour(auto i = 0; i < nfreqs; i++)
      ol[i] = source_ohc(((float) i) / nfreqs);

    index = 0;
  }
  void step(const Veccf &x, Veccf &y)
  {
    msg("Calcul FHSS...");
    soit &config = Configurable<FHSSConfig>::config;
    y = ra->step(x);
    msg("x.rows = {}, y.rows = {}.", x.rows(), y.rows());
    // Chaque slot est multiplié par une exponentielle différente
    pour(auto i = 0; i + config.duree_slot <= y.rows(); i += config.duree_slot)
    {
      //ol.configure(rotations[index++ % nfreqs]);
      soit osc = ol[i]->step(config.duree_slot);
      y.segment(i, config.duree_slot) *= osc;
    }
    // TODO : gestion si paquet pas multiple de ...
  }
};

sptr<Filtre<cfloat,cfloat,FHSSConfig>> fhss_modulation(const FHSSConfig &config)
{
  retourne std::make_shared<FHSSModulation>(config);
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
  entier osf = 0;
  entier nbits = 0, index = 0;

  Vecf seq;
  sptr<FiltreGen<cfloat>> ra;


  DSSSModulation(const DSSSConfig &c)
  {
    configure(c);
  }

  void configure_impl(const DSSSConfig &c)
  {
    soit &config = Configurable<DSSSConfig>::config;
    index   = 0;
    config  = c;
    nbits   = c.chips.rows();
    osf     = nbits / c.osf_in;
    si((nbits % c.osf_in) != 0)
      échec("DSSS : nbits doit être un multiple de osf in.");
    ra = tsd::filtrage::filtre_reechan<cfloat>(((float) nbits) / c.osf_in);
  }

  void step(const Veccf &x, Veccf &y)
  {
    soit &config = Configurable<DSSSConfig>::config;
    // pour chaque symbole d'entrée, consomme un bit de la séquence
    // Mais, un symbole d'entrée = plusieurs échantillons
    assertion((x.rows() % config.osf_in) == 0);
    //soit nsymb = x.rows() / config.osf_in;

    y = ra->step(x);

    // Maintenant, y est échantillonné à la fréquence bit (chips)
    pour(auto j = 0; j < y.rows(); j++)
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
  retourne std::make_shared<DSSSModulation>(config);
}


}











