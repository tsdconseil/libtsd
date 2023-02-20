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
  unsigned int ncanaux = 1;
  double freq_ech = 48000;
  enum
  {
    PCM_16 = 0,
    PCM_FLOAT
  } format = PCM_16;
};

/** @~french  @brief Chargement d'un fichier audio monovoie (format .wav)
 *  @returns Données, et fréquence d'échantillonnage
 *  @~english @brief Loading of a single channel .wav audio file*/
extern std::tuple<Vecf, float> wav_charge(const std::string &chemin);

/** @~french  @brief Chargement d'un fichier audio stéréo (format .wav)
 *  @returns Données, et fréquence d'échantillonnage
 *  @~english @brief Loading of a two channels (stereo) .wav audio file */
extern std::tuple<Veccf, float> wav_charge_stereo(const std::string &chemin);

/** @~french @brief Enregistre un signal réel en tant que fichier .wav monovoie
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern entier wav_enregistre(const std::string &chemin, float fech, const Vecf &x);

/** @brief Save a complex signal as a two channels  (stereo) .wav file
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern entier wav_enregistre_stereo(const std::string &chemin, float fech, const Veccf &x);

/** @cond undoc */
const auto wav_save = wav_enregistre;
/** @endcond */


/**  @~english @brief Save a real signal as a single channel .wav file
  *  @~english @param path Output file path,
  *  @~english @param fs Sampling frequency,
  *  @~english @param x Signal to be saved. */
inline entier wave_save(const std::string &path, float fs, const Vecf &x)
{
  return wav_enregistre(path, fs, x);
}


const auto wav_load = wav_charge;

struct WavLecteur
{
  WavLecteur();

  entier charge(const std::string &filename);

  /** @brief Duration in samples */
  uint32_t lis_nechantillons() const;

  /** @brief Duration in seconds */
  float lis_duree() const;

  entier lis_config(WavConfig &res);
  entier lis_donnees(int16_t *res, uint32_t n);
  Vecf lis_flottant(uint32_t n);

  //ArrayXXcf lis_cplx();

  entier redemarre();
  entier ferme();
  bouléen eof() const;
  uint32_t lis_position() const;
  float lis_position_secondes() const;

  entier recherche(float temps);

  std::string lis_infos() const;



private:
  struct Impl;
  sptr<Impl> impl;
};



struct WavEcrivain
{
  WavEcrivain();
  entier init(const std::string &chemin, const WavConfig &config);

  entier ecris(int16_t *donnees, uint32_t n);
  entier ecris(const Veccf &x);
  entier ecris(const Vecf &x);

  entier ferme();
  float lis_position_secondes() const;
  std::string lis_infos() const;
private:
  struct Impl;
  sptr<Impl> impl;
};



/** @} */

}



