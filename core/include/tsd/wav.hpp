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
extern std::tuple<ArrayXf, float> wav_charge(const std::string &chemin);

/** @~french  @brief Chargement d'un fichier audio stéréo (format .wav)
 *  @returns Données, et fréquence d'échantillonnage
 *  @~english @brief Loading of a two channels (stereo) .wav audio file */
extern std::tuple<ArrayXcf, float> wav_charge_stereo(const std::string &chemin);

/** @~french @brief Enregistre un signal réel en tant que fichier .wav monovoie
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern int wav_enregistre(const std::string &chemin, float fech, const ArrayXf &x);

/** @brief Save a complex signal as a two channels  (stereo) .wav file
 *  @param chemin Output file path,
 *  @param fech Sampling frequency,
 *  @param x Signal to be saved.
 *  @returns Non-zero value in case of error.
 */
extern int wav_enregistre_stereo(const std::string &chemin, float fech, const ArrayXcf &x);

/** @cond undoc */
const auto wav_save = wav_enregistre;
/** @endcond */


/**  @~english @brief Save a real signal as a single channel .wav file
  *  @~english @param path Output file path,
  *  @~english @param fs Sampling frequency,
  *  @~english @param x Signal to be saved. */
inline int wave_save(const std::string &path, float fs, const ArrayXf &x)
{
  return wav_enregistre(path, fs, x);
}


const auto wav_load = wav_charge;

class WavLecteur
{
public:
  WavLecteur();

  int charge(const std::string &filename);

  /** @brief Duration in samples */
  uint32_t lis_nechantillons() const;

  /** @brief Duration in seconds */
  float lis_duree() const;

  int lis_config(WavConfig &res);
  int lis_donnees(int16_t *res, uint32_t n);
  ArrayXf lis_flottant(uint32_t n);

  //ArrayXXcf lis_cplx();

  int redemarre();
  int ferme();
  bool eof() const;
  uint32_t lis_position() const;
  float lis_position_secondes() const;

  int recherche(float temps);

  std::string lis_infos() const;



private:
  struct Impl;
  sptr<Impl> impl;
};



class WavEcrivain
{
public:
  WavEcrivain();
  int init(const std::string &chemin, const WavConfig &config);

  int ecris(int16_t *donnees, uint32_t n);
  int ecris(const ArrayXcf &x);
  int ecris(const ArrayXf &x);

  int ferme();
  float lis_position_secondes() const;
  std::string lis_infos() const;
private:
  struct Impl;
  sptr<Impl> impl;
};

/** @} */

}



