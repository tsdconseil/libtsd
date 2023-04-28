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
  string texte;
};

struct POCSAGDecodeurConfig
{
  bouléen debug_actif = non;
  float fe = 1, fi = 0;

  // Normalement, 512, 1200 ou 2400 bauds
  // -1 = détermination automatique
  entier debit = -1;

};

struct POCSAGMessage
{
  uint32_t ric;
  uint16_t function;
  string texte;
};


/** @brief Création d'un démodulateur ADSB */
extern sptr<ProtocoleDemodulateur<ADSBDecodeurConfig, ADSBTrame>> demodulateur_adsb();

/** @brief Création d'un démodulateur POCSAG */
extern sptr<ProtocoleDemodulateur<POCSAGDecodeurConfig, POCSAGMessage>> demodulateur_pocsag();

/** @} */


}
