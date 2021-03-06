#pragma once

#include "tsd/telecom.hpp"
#include "tsd/telecom/bitstream.hpp"

namespace tsd::telecom {

/** @addtogroup telecom-protos
 *  @{
 */

struct ADSBDecodeurConfig
{
  float fe;
  unsigned int Ne;
};

struct ADSBTrame
{
  BitStream bs;
  float score;
  std::string texte;
};

struct POCSAGDecodeurConfig
{
  bool debug_actif = false;
  float fe = 1, fi = 0;

  // Normalement, 512, 1200 ou 2400 bauds
  // -1 = détermination automatique
  int debit = -1;

};

struct POCSAGMessage
{
  uint32_t ric;
  uint16_t function;
  std::string texte;
};


/** @brief Création d'un démodulateur ADSB */
extern sptr<ProtocoleDemodulateur<ADSBDecodeurConfig, ADSBTrame>> demodulateur_adsb();

/** @brief Création d'un démodulateur POCSAG */
extern sptr<ProtocoleDemodulateur<POCSAGDecodeurConfig, POCSAGMessage>> demodulateur_pocsag();

/** @} */


}
