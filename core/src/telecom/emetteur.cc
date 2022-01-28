#include "tsd/telecom.hpp"

using namespace std;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

#define VERB(AA)
//AA



namespace tsd::telecom
{


struct ÉmetteurImpl: Émetteur
{
  ÉmetteurConfig config;
  sptr<Modulateur> mod;

  float retard() const
  {
    // Temps vers le milieu du premier bit transmis
    return tsd::filtrage::filtre_rif_ups_délais(
        config.format.modulation.ncoefs_filtre_mise_en_forme,
        config.format.modulation.fe / config.format.modulation.fsymb);
  }

  ÉmetteurImpl(const ÉmetteurConfig &config)
  {
    configure(config);
  }

  int configure(const ÉmetteurConfig &config)
  {
    this->config = config;
    mod = modulateur_création(config.format.modulation);

    if(config.format.fo_entete)
    {
      msg_majeur("Émetteur: fo entete = {}", *config.format.fo_entete);
      msg_majeur("Émetteur: fo données = {}", *config.format.modulation.wf);
    }

    return 0;
  }

  ArrayXcf step(const BitStream &bs)
  {
    BitStream bs2;

    bool fo_entete_specifique = config.format.fo_entete ? true : false;

    if(fo_entete_specifique)
      bs2 = bs;
    else
      bs2 = config.format.entete + bs;

    // Ajoute des zéros pour flusher
    auto &mconfig = config.format.modulation;
    auto osf = mconfig.fe / mconfig.fsymb;

    int   d_ech           = (int) ceil(mod->delais());
    int   nbits_par_symb  = mconfig.wf->k;
    //float nbits_par_echan = ((float) nbits_par_symb) / osf;
    //int   d_bit           = (int) ceil(d_ech * nbits_par_echan);


    //msg("Délais = {} échantillons, k = {} bits/symb, {} bits/échan => flush {} bits.",
    //    mod->delais(), config.forme_onde->k,
    //    nbits_par_echan, d_bit);

    // Pour le flush du filtre
    //for(auto i = 0; i < d_bit; i++)
      //bs2.push(0);
    // PB : ceci termine avec une séquence I/Q constante non nulle...

    ArrayXcf x;
    if(fo_entete_specifique)
    {
      nbits_par_symb  = config.format.fo_entete->k;
      mod->def_fo(config.format.fo_entete);

      BitStream tmp = config.format.entete;
      tmp.pad_mult(nbits_par_symb);

      ArrayXcf x1 = mod->step(tmp);

      nbits_par_symb  = mconfig.wf->k;
      bs2.pad_mult(nbits_par_symb);
      mod->def_fo(config.format.modulation.wf);
      ArrayXcf x2 = mod->step(bs2);
      x = vconcat(x1, x2);
    }
    else
    {
      bs2.pad_mult(nbits_par_symb);
      x = mod->step(bs2);
    }

    // Flush du modulateur
    x = vconcat(x, mod->flush(d_ech));

    // TODEL
    // Flush sur 4 symboles
    //x = vconcat(x, mod->flush(4.0f / nbits_par_echan));

    return x;
  }

  MoniteursStats moniteurs()
  {
    return {};
  }
};

sptr<Émetteur> émetteur_création(const ÉmetteurConfig &ec)
{
  return std::make_shared<ÉmetteurImpl>(ec);
}



}
