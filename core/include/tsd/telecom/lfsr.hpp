#pragma once

#include "tsd/tsd.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"

#include <vector>

namespace tsd::telecom {




/** @addtogroup telecom-ber
 *  @{
 */






/** @brief Structure de configuration générateur LFSR */
struct LFSRConfig
{
  int reglen = 0;
  uint32_t pol    = 0;
  uint32_t p0     = 0;
  uint32_t pol_sortie = 0;
  // Sinon, sortie = poids fort
  // bool sortie_est_poids_faible = false;
  // bool sortie_est_pol = false;

  enum
  {
    POIDS_FAIBLE,
    POIDS_FORT,
    POL
  } sortie = POIDS_FAIBLE;

};

/** @brief PRBS generator context */
class LFSRGenerateur
{
public:

  LFSRGenerateur(unsigned int reglen = 16);

  /** @brief PRBS generator context initialization
   *  @param reglen: length of the PRBS register
   *  @note Sequence length will be (2^reglen)-1
   *  Requires 2 <= reglen <= PRBS_MAX_REGLEN */
  int configure(unsigned int reglen);

  int configure(const LFSRConfig &config);

  /** @brief Generate PRBS data */
  void step(BitStream &bs, unsigned int nbits);

private:
  LFSRConfig config;
  //uint16_t reglen;
  uint32_t reg;//, pol;
};


/** @brief LFSR decoder context */
class LFSRRecepteur
{
public:
  /** @brief Constructeur */
  LFSRRecepteur();

  /** @brief PRBS receiver context initialization
   *  @param reglen: length of the PRBS register
   *  @note Sequence length will be (2^reglen)-1
   *  Requires 2 <= reglen <= PRBS_MAX_REGLEN
   *  @param  nb_bits_to_ignore Les premiers bits ne sont pas pris en compte pour le calcul du BER. */
  int configure(uint16_t reglen, int nb_bits_to_ignore = 0);

  void step(const BitStream &bs);

  /** @brief Reset the PRBS receiver in unlocked state */
  void reset();

  /** @param[out] is_locked
   *  @param[out] ber        Current bit error rate, or 0.5 if not locked. */
  void lis_etat(bool &is_locked, float &ber) const;

  void affiche_etat() const;

private:
  struct Impl;
  sptr<Impl> impl;
};

/** @} */

}

