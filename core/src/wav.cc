#include "tsd/tsd.hpp"
#include "tsd/wav.hpp"
#include <cstdio>
extern "C"
{
# include <unistd.h>
}

namespace tsd::audio {


  typedef struct wave_header_struct
  {
    /* Constante "RIFF"  (0x52,0x49,0x46,0x46) */
    uint8_t file_type_bloc_id[4];

    /* Taille du fichier moins 8 octets */
    uint32_t file_size;

    /* Format = "WAVE"  (0x57,0x41,0x56,0x45) */
    uint8_t file_format_id[4];

    /** [Bloc d�crivant le format audio] */

    /* Identifiant "fmt" (0x66,0x6D, 0x74,0x20) */
    uint8_t format_bloc_id[4];

    /* Nombre d'octets du bloc - 8 */
    uint32_t bloc_size;

    /* Format du stockage dans le fichier (1: PCM, 3 : FLOAT32) */
    uint16_t audio_format;

    /* Nombre de canaux (de 1 � 6) */
    uint16_t nbr_canaux;

    /* Fréquence d'échantillonnage (en Hertz) [Valeurs standardis�es : 11025, 22050, 44100 et �ventuellement 48000 et 96000] */
    uint32_t frequency;

    /* Nombre d'octets à lire par seconde (i.e., Frequence * BytePerBloc). */
    uint32_t byte_per_sec;

    /* Nombre d'octets par bloc d'échantillonnage (i.e., tous canaux confondus : NbrCanaux * BitsPerSample/8). */
    uint16_t byte_per_bloc;

    /* Nombre de bits utilis�es pour le codage de chaque �chantillon (8, 16, 24) */
    uint16_t bits_per_sample;

    /** [Bloc des données] */

    /* Constante data  (0x64,0x61,0x74,0x61) */
    uint8_t data_bloc_id[4];

    /* Nombre d'octets des donn�es (i.e. "Data[]", i.e. taille_du_fichier - taille_de_l'ent�te  (qui fait 44 octets normalement). */
    uint32_t data_size;

    /*   DATAS[] : [Octets du Sample 1 du Canal 1] [Octets du Sample 1 du Canal 2]
     *             [Octets du Sample 2 du Canal 1] [Octets du Sample 2 du Canal 2]
     *             .... */
  } __attribute__((__packed__)) wave_entete_t;


std::tuple<ArrayXcf, float> wav_charge_stereo(const std::string &chemin)
{
  WavLecteur lecteur;
  if(lecteur.charge(chemin))
    abort();
  WavConfig cfg;
  lecteur.lis_config(cfg);

  auto n = lecteur.lis_nechantillons();
  ArrayXf x = lecteur.lis_flottant(2 * n);

  tsd_assert(x.rows() == (int) (2 *n));

  ArrayXcf y = Eigen::Map<ArrayXcf>((cfloat *) x.data(), n);

  return {y, (float) cfg.freq_ech};
}

int wav_enregistre(const std::string &chemin, float fech, const ArrayXf &x)
{
  WavEcrivain wa;
  wa.init(chemin, {1, fech, WavConfig::PCM_FLOAT});
  wa.ecris(x);
  return 0;
}

int wav_enregistre_stereo(const std::string &chemin, float fech, const ArrayXcf &x)
{
  WavEcrivain wa;
  wa.init(chemin, {2, fech, WavConfig::PCM_FLOAT});
  wa.ecris(x);
  return 0;
}

std::tuple<ArrayXf, float> wav_charge(const std::string &chemin)
{
  WavLecteur lecteur;
  if(lecteur.charge(chemin))
    abort();

  auto n = lecteur.lis_nechantillons();
  WavConfig cfg;
  lecteur.lis_config(cfg);


  auto buf = (int16_t *) malloc(n * sizeof(int16_t));
  if(buf == nullptr)
  {
    msg_erreur("Erreur d'allocation ({} éléments).", n);
    abort();
  }

  lecteur.lis_donnees(buf, n);

  ArrayXf x(n);
  for(auto i = 0u; i < n; i++)
    x(i) = buf[i];

  free(buf);

  return {x, (float) cfg.freq_ech};
}

struct WavLecteur::Impl
{
  std::string chemin;
  FILE *fd = nullptr;
  uint32_t pos = 0, nech = 0;

  wave_entete_t entete;

  ~Impl()
  {
    ferme();
  }

  int charge(const std::string &chemin_)
  {
    ferme();

    chemin = chemin_;
    msg("Ouverture fichier wav, chemin = [{}]...", chemin.c_str());

    //cutils::fichiers::remplacement_motif(chemin);
    pos = 0;

    /*if(!cutils::fichiers::fichier_existe(chemin))
    {
      erreur("Fichier non trouve : [%s].", chemin.c_str());
      return -1;
    }*/

    fd = fopen(chemin.c_str(), "rb");
    if(fd == nullptr)
    {
      msg_erreur("Echec ouverture du fichier wav.");
      return -1;
    }
    if(fread(&entete, 1, sizeof(entete), fd) != sizeof(entete))
    {
      msg_erreur("Read header failed.");
      ferme();
      return -1;
    }

    auto v32 = *((uint32_t *) entete.file_type_bloc_id);

    if(v32 != 0x46464952)
    {
      msg_erreur("En-tête incorrect : {}", v32);
      ferme();
      return -1;
    }
    retry:
    v32 = *((uint32_t *) entete.data_bloc_id);
    if(v32 != 0x61746164)
    {
      char c[5];
      c[4] = 0;
      memcpy(c, entete.data_bloc_id, 4);
      msg_avert("Bloc de données ignoré : {:x} [{}], dim = {}", v32, c, (int) entete.data_size);
      msg("FTELL = {}", ftell(fd));
      fseek(fd, entete.data_size, SEEK_CUR);
      if(feof(fd))
      {
        msg_erreur("Champs data non trouvé.");
        ferme();
        return -1;
      }
      if(fread(&(entete.data_bloc_id), 8, 1, fd) != 1)
      {
        echec("Erreur lecture fichier WAV.");
        return -1;
      }
      goto retry;
    }
    nech = (8 * entete.data_size) / (entete.bits_per_sample * entete.nbr_canaux);
    msg("Nb canaux = {}, freq. ech={}, bits/ech={}",
        (int) entete.nbr_canaux, (int) entete.frequency, (int) entete.bits_per_sample);
    msg("Format = {} ({})", (int) entete.audio_format, entete.audio_format == 1 ? "PCM 16 bits" : entete.audio_format == 3 ? "PCM float32" : "?");
    msg("Nb echantillons au total (par canal) = {}.", nech);
    msg("Durée {:.2f} secondes ({} octets).",
        ((float) nech) / entete.frequency,
        (int) entete.data_size);

    return 0;
  }

  int ferme()
  {
    if(fd != nullptr)
    {
      fclose(fd);
      fd = nullptr;
    }
    return 0;
  }

  ArrayXf lis_flottant(uint32_t n)
  {
    if(fd == nullptr)
    {
      msg_erreur("non ouvert.");
      return ArrayXf();
    }
    ArrayXf y(n);

    if(entete.audio_format == 1)
    {
      // Format 16 bits
      int16_t *t16 = (int16_t *) malloc(2*n);
      lis_donnees(t16, n);
      for(auto i = 0u; i < n; i++)
        y(i) = t16[i];
      free(t16);
      return y;
    }


    int nlus = fread(y.data(), 4, n, fd);
    if(nlus != (int) n)
    {
      msg_avert("fread: nlus = {}, ndemandés = {} (position initiale = {}).", nlus, n, pos);
    }

    pos += nlus;
    return y;
  }

  int lis_donnees(int16_t *res, uint32_t n)
  {
    if(fd == nullptr)
    {
      msg_erreur("non ouvert.");
      return -1;
    }

    //trace_verbeuse("n = %d, pos = %d / %d", n, pos, nech);

    int nlus = fread(res, 2, n, fd);
    if(nlus != (int) n)
    {
      msg_avert("nlus = {}, n = {} (pos = {}).", nlus, n, pos);
    }

    pos += nlus;

    //infos("pos += %d -> %d", nlus, pos);

    return nlus;
  }
};

WavLecteur::WavLecteur()
{
  impl = std::make_shared<Impl>();
}




int WavLecteur::ferme()
{
  return impl->ferme();
}

int WavLecteur::charge(const std::string &chemin_)
{
  return impl->charge(chemin_);
}

float WavLecteur::lis_duree() const
{
  if(impl->fd == nullptr)
    return 0;
  return ((float) impl->nech) / impl->entete.frequency;
}

int WavLecteur::lis_config(WavConfig &res)
{
  res.freq_ech = impl->entete.frequency;
  res.ncanaux  = impl->entete.nbr_canaux;
  return 0;
}

std::string WavLecteur::lis_infos() const
{
  if(impl->fd == nullptr)
    return "(aucun fichier wav ouvert)";
  return format("{}, {} cana{}, {} Hz, {:.1f} sec.",
      impl->chemin.c_str(),
      (int) impl->entete.nbr_canaux,
      impl->entete.nbr_canaux > 1 ? "ux" : "l",
          (int) impl->entete.frequency,
          ((float) impl->nech) / impl->entete.frequency);
}

uint32_t WavLecteur::lis_nechantillons() const
{
  if(impl->fd == nullptr)
    return 0;
  return impl->nech;
}

bool WavLecteur::eof() const
{
  return (impl->fd == nullptr) || (feof(impl->fd));
}

uint32_t WavLecteur::lis_position() const
{
  return impl->pos;
}

float WavLecteur::lis_position_secondes() const
{
  return ((float) impl->pos) / impl->entete.frequency;
}

/*ArrayXXcf WavLecteur::lis_cplx()
{
  return impl->lis_cplx();
}*/

ArrayXf WavLecteur::lis_flottant(uint32_t n)
{
  return impl->lis_flottant(n);
}

int WavLecteur::lis_donnees(int16_t *res, uint32_t n)
{
  return impl->lis_donnees(res, n);
}

int WavLecteur::redemarre()
{
  if(impl->fd == nullptr)
  {
    msg_erreur("non ouvert.");
    return -1;
  }
  msg("Remise a zero lecteur wav.");
  fseek(impl->fd, sizeof(impl->entete), SEEK_SET);
  impl->pos = 0;
  return 0;
}

int WavLecteur::recherche(float temps)
{
  if(impl->fd == nullptr)
  {
    msg_erreur("non ouvert.");
    return -1;
  }

  impl->pos = (int) (temps * impl->entete.frequency);

  msg("Recherche, pos = {} secondes ({} echan / {} total) ({} %)", temps, impl->pos, impl->nech, ((float) impl->pos) / impl->nech);

  fseek(impl->fd, impl->pos * 2 + sizeof(impl->entete), SEEK_SET);

  return 0;
}

struct WavEcrivain::Impl
{
  std::string chemin;
  WavConfig cfg;
  FILE *fd = nullptr;
  uint32_t pos = 0;

  ~Impl()
  {
    ferme();
  }

  int init(const std::string &chemin, const WavConfig &config)
  {
    wave_entete_t entete;

    ferme();

    this->chemin = chemin;

    msg("Ouverture fichier wav en ecriture [{}], ncanaux = {}, fech = {:.1f} kHz.",
          chemin.c_str(), config.ncanaux, config.freq_ech * 1e-3f);

    /*auto ext = cutils::fichiers::extrait_extension(chemin);
    if(ext != "wav")
    {
      avertissement("Ceci n'est pas un fichier WAV (extension = %s).", ext.c_str());
      return ERREUR_WAV_EXTENSION;
    }*/

    fd = fopen(chemin.c_str(), "wb");
    if(fd == nullptr)
    {
      msg_erreur("Echec ouverture fichier wav en ecriture ({}).", chemin);
      return ERREUR_WAV_FICHIER_INACCESSIBLE;
    }
    // génération d'un entete vide, à remplir à la fermeture du fichier
    memset(&entete, 0, sizeof(wave_entete_t));
    fwrite(&entete, sizeof(wave_entete_t), 1, fd);

    /* copie de la config passï¿½e */
    this->cfg = config;

    pos = 0;
    return 0;
  }

  int ferme()
  {
    if(fd != nullptr)
    {
      msg("Fermeture fichier wav, nech = {}...", pos);

      bool flottant = cfg.format == WavConfig::PCM_FLOAT;
      int octets_par_echantillon = flottant ? 4 : 2;

      wave_entete_t entete;
      // Ecriture de l'en-tête
      entete.audio_format = flottant ? 3 : 1;
      entete.bloc_size = 16;
      entete.file_type_bloc_id[0] = 'R';
      entete.file_type_bloc_id[1] = 'I';
      entete.file_type_bloc_id[2] = 'F';
      entete.file_type_bloc_id[3] = 'F';
      entete.file_format_id[0] = 'W';
      entete.file_format_id[1] = 'A';
      entete.file_format_id[2] = 'V';
      entete.file_format_id[3] = 'E';
      entete.format_bloc_id[0] = 'f';
      entete.format_bloc_id[1] = 'm';
      entete.format_bloc_id[2] = 't';
      entete.format_bloc_id[3] = 32;
      entete.nbr_canaux       = cfg.ncanaux;
      entete.frequency        = cfg.freq_ech;
      entete.byte_per_sec     = cfg.freq_ech * octets_par_echantillon;
      entete.byte_per_bloc    = octets_par_echantillon;
      entete.bits_per_sample  = octets_par_echantillon * 8;
      entete.data_bloc_id[0] = 'd';
      entete.data_bloc_id[1] = 'a';
      entete.data_bloc_id[2] = 't';
      entete.data_bloc_id[3] = 'a';
      entete.data_size = pos * octets_par_echantillon * cfg.ncanaux;
      fseek(fd, 0, SEEK_SET);
      fwrite(&entete, sizeof(wave_entete_t), 1, fd);
      fclose(fd);
      fd = nullptr;
    }
    return 0;
  }

  int ecris(const ArrayXcf &x)
  {
    //fwrite(x.data(), 8, x.rows(), fd);
    //pos += x.rows();// / cfg.ncanaux;

    ArrayXf x2 = Eigen::Map<ArrayXf>((float *) x.data(), x.rows() * 2);
    return ecris(x2);
  }

  int ecris(const ArrayXf &x)
  {
    //infos("Ecriture %d ech float...", x.rows());
    if(fd != nullptr)
    {
      if(cfg.format == WavConfig::PCM_16)
      {
        auto p = (int16_t *) malloc(2 * x.rows());
        for(auto i = 0; i < x.rows(); i++)
          p[i] = (int16_t) x(i);
        ecris(p, x.rows());
        free(p);
      }
      else
      {
        fwrite(x.data(), 4, x.rows(), fd);
        pos += x.rows() / cfg.ncanaux;
      }
    }
    return 0;
  }

  int ecris(int16_t *donnees, uint32_t n)
  {
    if(fd != nullptr)
    {
      if(cfg.format == WavConfig::PCM_16)
      {
        fwrite(donnees, 2, n, fd);
        pos += n / cfg.ncanaux;
      }
      else
      {
        ArrayXf x(n);
        for(auto i = 0u; i < n; i++)
          x(i) = donnees[i];
        ecris(x);
      }
    }
    return 0;
  }

};


WavEcrivain::WavEcrivain()
{
  impl = std::make_shared<Impl>();
}

std::string WavEcrivain::lis_infos() const
{
  float duree = ((float) impl->pos) / impl->cfg.freq_ech;
  return impl->chemin + format(" ({} cana{}, {:.1f} kHz, {:.1f} secondes)",
      impl->cfg.ncanaux,
      impl->cfg.ncanaux > 1 ? "ux" : "l",
          (float) (impl->cfg.freq_ech * 1e-3), duree);
}

int WavEcrivain::init(const std::string &chemin, const WavConfig &config)
{
  return impl->init(chemin, config);
}

float WavEcrivain::lis_position_secondes() const
{
  return ((float) impl->pos) / impl->cfg.freq_ech;
}

int WavEcrivain::ecris(int16_t *donnees, uint32_t n)
{
  return impl->ecris(donnees, n);
}

int WavEcrivain::ecris(const ArrayXf &x)
{
  return impl->ecris(x);
}

int WavEcrivain::ecris(const ArrayXcf &x)
{
  return impl->ecris(x);
}

int WavEcrivain::ferme()
{
  return impl->ferme();
}




}

