#include "tsd/telecom.hpp"

using namespace std;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

#define VERB(AA)

namespace tsd::telecom
{


struct ÉmetteurImpl: Émetteur
{
  ÉmetteurConfig config;
  sptr<Modulateur> mod;

  float retard() const
  {
    auto &mod = config.format.modulation;
    // Temps vers le milieu du premier bit transmis
    retourne tsd::filtrage::filtre_rif_ups_délais(
        mod.ncoefs_filtre_mise_en_forme,
        mod.fe / mod.fsymb);
  }

  ÉmetteurImpl(const ÉmetteurConfig &config)
  {
    configure(config);
  }

  entier configure(const ÉmetteurConfig &config)
  {
    this->config = config;
    mod = modulateur_création(config.format.modulation);

    si(config.format.fo_entete)
    {
      msg_majeur("Émetteur: fo entete = {}", *config.format.fo_entete);
      msg_majeur("Émetteur: fo données = {}", *config.format.modulation.forme_onde);
    }

    retourne 0;
  }

  Veccf step(const BitStream &bs)
  {
    BitStream bs2;
    soit &mconfig = config.format.modulation;
    soit fo_entete_specifique = config.format.fo_entete ? oui : non;

    si(fo_entete_specifique)
      bs2 = bs;
    sinon
    {
      BitStream et2 = config.format.entete;
      et2.pad_mult(mconfig.forme_onde->infos.k);
      bs2 = et2 + bs;
    }

    // Ajoute des zéros pour flusher

    //soit osf = mconfig.fe / mconfig.fsymb;

    soit  d_ech           = (entier) ceil(mod->delais());
    soit  nbits_par_symb  = mconfig.forme_onde->infos.k;
    //float nbits_par_echan = ((float) nbits_par_symb) / osf;
    //entier   d_bit           = (entier) ceil(d_ech * nbits_par_echan);


    //msg("Délais = {} échantillons, k = {} bits/symb, {} bits/échan => flush {} bits.",
    //    mod->delais(), config.forme_onde->k,
    //    nbits_par_echan, d_bit);

    // pour le flush du filtre
    //pour(soit i = 0; i < d_bit; i++)
      //bs2.push(0);
    // PB : ceci termine avec une séquence I/Q constante non nulle...

    Veccf x;
    si(fo_entete_specifique)
    {
      nbits_par_symb  = config.format.fo_entete->infos.k;
      mod->def_forme_onde(config.format.fo_entete);

      BitStream tmp = config.format.entete;
      tmp.pad_mult(nbits_par_symb);

      soit x1 = mod->step(tmp);

      nbits_par_symb  = mconfig.forme_onde->infos.k;
      bs2.pad_mult(nbits_par_symb);
      mod->def_forme_onde(config.format.modulation.forme_onde);
      soit x2 = mod->step(bs2);
      x = vconcat(x1, x2);
    }
    sinon
    {
      bs2.pad_mult(nbits_par_symb);
      x = mod->step(bs2);
    }

    // Flush du modulateur
    x = vconcat(x, mod->flush(d_ech));

    // TODEL
    // Flush sur 4 symboles
    //x = vconcat(x, mod->flush(4.0f / nbits_par_echan));

    retourne x;
  }

  MoniteursStats moniteurs()
  {
    retourne {};
  }
};

sptr<Émetteur> émetteur_création(const ÉmetteurConfig &ec)
{
  retourne make_shared<ÉmetteurImpl>(ec);
}



}
