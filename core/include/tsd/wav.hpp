#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"


namespace tsd::audio {

  /** @addtogroup audio
   *  @{
   */

#define ERREUR_WAV_FICHIER_INACCESSIBLE 0x2000
#define ERREUR_WAV_EXTENSION            0x2001

struct WavConfig
{
  entier ncanaux = 1;
  double freq_ech = 48000;
  enum
  {
    PCM_16 = 0,
    PCM_FLOAT
  } format = PCM_16;
};

/** @brief Chargement d'un fichier audio monovoie (format .wav)
 *  @returns Données, et fréquence d'échantillonnage
 *  @~english @brief Loading of a single channel .wav audio file*/
extern tuple<Vecf, float> wav_charge(cstring chemin);

/** @~french  @brief Chargement d'un fichier audio stéréo (format .wav)
 *  @returns Données, et fréquence d'échantillonnage
 *  @~english @brief Loading of a two channels (stereo) .wav audio file */
extern tuple<Veccf, float> wav_charge_stereo(cstring schemin);

/** @~french @brief Enregistre un signal réel en tant que fichier .wav monovoie
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern void wav_enregistre(cstring chemin, float fech, const Vecf &x);

/** @brief Save a complex signal as a two channels  (stereo) .wav file
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern void wav_enregistre_stereo(cstring chemin, float fech, const Veccf &x);

/** @cond undoc */
const auto wav_save = wav_enregistre;
/** @endcond */


/**  @~english @brief Save a real signal as a single channel .wav file
  *  @~english @param path Output file path,
  *  @~english @param fs Sampling frequency,
  *  @~english @param x Signal to be saved. */
inline void wave_save(cstring path, float fs, const Vecf &x)
{
  wav_enregistre(path, fs, x);
}


const auto wav_load = wav_charge;

struct WavLecteur
{
  WavLecteur();

  void charge(cstring filename);

  /** @brief Duration in samples */
  entier lis_nechantillons() const;

  /** @brief Duration in seconds */
  float lis_duree() const;

  WavConfig lis_config() const;

  /** @returns nb échantillons lus */
  entier lis_donnees(int16_t *res, entier n);

  Vecf lis_flottant(entier n);

  void redemarre();
  void ferme();

  bouléen eof() const;

  entier lis_position() const;

  float lis_position_secondes() const;

  void recherche(float temps);

  string lis_infos() const;
  _PIMPL_
};



struct WavEcrivain
{
  WavEcrivain();
  void init(cstring chemin, const WavConfig &config);

  void ecris(int16_t *donnees, entier n);
  void ecris(const Veccf &x);
  void ecris(const Vecf &x);

  void ferme();
  float lis_position_secondes() const;
  string lis_infos() const;
  _PIMPL_
};



/** @} */

}



