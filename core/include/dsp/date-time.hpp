#pragma once

#include "tsd/date-heure.hpp"

namespace dsp::time {

/** @addtogroup tsd-heure
 *  @{
 */


/** @brief Represents a time interval, as number of micro-seconds.
 *
 */
using TimeSpan = tsd::temps::IntervalleTemps;

using Calendar      = tsd::temps::Calendrier;

using DateComposite = tsd::temps::DateComposite;

using TimeComposite = tsd::temps::HeureComposite;

/** @brief Date time representation, with conversions to / from sidereal time, UTC, ... */
using DateTime = tsd::temps::DateHeure;

// Non, reprendre la classe, toutes ses méthodes, etc.


/** @} */

}
