#include "tsd/tsd.hpp"
#include "tsd/wav.hpp"
#include <cstdio>
#include <cstring>
extern "C"
{
# include <unistd.h>
}

using namespace std;

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

    /** [Blocs des données] */

    /* Type de bloc ("data" = 0x64,0x61,0x74,0x61 pour les échantillons) */
    uint8_t data_bloc_id[4];

    /* Nombre d'octets des donn�es (i.e. "Data[]", i.e. taille_du_fichier - taille_de_l'ent�te  (qui fait 44 octets normalement). */
    uint32_t data_size;

    /*   DATAS[] : [Octets du Sample 1 du Canal 1] [Octets du Sample 1 du Canal 2]
     *             [Octets du Sample 2 du Canal 1] [Octets du Sample 2 du Canal 2]
     *             .... */
  } __attribute__((__packed__)) wave_entete_t;


tuple<Veccf, float> wav_charge_stereo(cstring chemin)
{
  WavLecteur lecteur;
  lecteur.charge(chemin);
  soit cfg = lecteur.lis_config();

  soit n = lecteur.lis_nechantillons();
  soit x = lecteur.lis_flottant(2 * n);

  assertion(x.rows() == 2 *n);

  Veccf y(n);
  memcpy(y.rawptr(), x.rawptr(), 2 * n * sizeof(float));

  retourne {y, (float) cfg.freq_ech};
}

void wav_enregistre(cstring chemin, float fech, const Vecf &x)
{
  WavEcrivain wa;
  wa.init(chemin, {1, fech, WavConfig::PCM_FLOAT});
  wa.ecris(x);
}

void wav_enregistre_stereo(cstring chemin, float fech, const Veccf &x)
{
  WavEcrivain wa;
  wa.init(chemin, {2, fech, WavConfig::PCM_FLOAT});
  wa.ecris(x);
}

tuple<Vecf, float> wav_charge(cstring chemin)
{
  WavLecteur lecteur;
  lecteur.charge(chemin);

  soit n   = lecteur.lis_nechantillons();
  soit cfg = lecteur.lis_config();


  soit buf = (int16_t *) malloc(n * sizeof(int16_t));
  si(!buf)
    échec("Erreur d'allocation ({} éléments).", n);

  soit x = lecteur.lis_flottant(n);

  retourne {x, (float) cfg.freq_ech};
}

struct WavLecteur::Impl
{
  string chemin;
  FILE *fd = nullptr;
  entier pos = 0, nech = 0;
  wave_entete_t entete;

  ~Impl()
  {
    ferme();
  }

  void charge(cstring chemin_)
  {
    ferme();

    chemin = chemin_;
    //msg("Ouverture fichier wav, chemin = [{}]...", chemin);

    pos = 0;

    /*si(!cutils::fichiers::fichier_existe(chemin))
    {
      erreur("Fichier non trouve : [%s].", chemin);
      retourne -1;
    }*/

    fd = fopen(chemin.c_str(), "rb");
    si(!fd)
      échec("Echec ouverture du fichier wav.");
    si(fread(&entete, 1, sizeof(entete), fd) != sizeof(entete))
    {
      ferme();
      échec("Read header failed.");
    }

    soit v32 = *((uint32_t *) entete.file_type_bloc_id);

    si(v32 != 0x46464952)
    {
      ferme();
      échec("En-tête incorrect : {}", v32);
    }
    retry:
    v32 = *((uint32_t *) entete.data_bloc_id);
    si(v32 != 0x61746164) // "data" en ASCII
    {
      // Position avant : début des données

      char c[5];
      c[4] = 0;
      memcpy(c, entete.data_bloc_id, 4);
      //msg_avert("Bloc de données ignoré : {:x} [{}], dim = {}", v32, c, (entier) entete.data_size);
      //msg("FTELL = {0} = 0x{0:x}", ftell(fd));

      entier nskips = entete.data_size;

      // Alignemnt 16 bits
      si(nskips & 1)
        nskips++;

      fseek(fd, nskips, SEEK_CUR);

      // Position après : début entête bloc suivant

      //msg("FTELL début bloc suivant = {0} = 0x{0:x}", ftell(fd));

      si(feof(fd))
      {
        ferme();
        échec("Champs data non trouvé.");
      }
      si(fread(&(entete.data_bloc_id), 8, 1, fd) != 1)
        échec("Erreur lecture fichier WAV.");

      // Position après : début données

      goto retry;
    }
    nech = (8 * entete.data_size) / (entete.bits_per_sample * entete.nbr_canaux);
    /*msg("Nb canaux = {}, freq. ech={}, bits/ech={}",
        (entier) entete.nbr_canaux, (entier) entete.frequency, (entier) entete.bits_per_sample);
    msg("Format = {} ({})", (entier) entete.audio_format, entete.audio_format == 1 ? "PCM 16 bits" : entete.audio_format == 3 ? "PCM float32" : "?");
    msg("Nb echantillons au total (par canal) = {}.", nech);
    msg("Durée {:.2f} secondes ({} octets).",
        ((float) nech) / entete.frequency,
        (entier) entete.data_size);*/
  }

  void ferme()
  {
    si(fd)
    {
      fclose(fd);
      fd = nullptr;
    }
  }

  Vecf lis_flottant(entier n)
  {
    assertion_msg(fd, "fichier wav non ouvert.");

    Vecf y(n);

    si(entete.audio_format == 1)
    {
      // Format 16 bits
      soit t16 = (int16_t *) malloc(2*n);
      lis_donnees(t16, n);
      pour(auto i = 0; i < n; i++)
        y(i) = t16[i];
      free(t16);
      retourne y;
    }

    entier nlus = fread(y.data(), 4, n, fd);
    si(nlus != n)
      msg_avert("fread: nlus = {}, ndemandés = {} (position initiale = {}).", nlus, n, pos);

    pos += nlus;
    retourne y;
  }

  entier lis_donnees(int16_t *res, entier n)
  {
    assertion_msg(fd, "fichier wav non ouvert.");

    //msg_verb("n = {}, pos = {} / {}", n, pos, nech);

    entier nlus = fread(res, 2, n, fd);
    si(nlus != n)
      msg_avert("nlus = {}, n = {} (pos = {}).", nlus, n, pos);

    pos += nlus;
    retourne nlus;
  }
};

WavLecteur::WavLecteur()
{
  impl = make_shared<Impl>();
}


void WavLecteur::ferme()
{
  impl->ferme();
}

void WavLecteur::charge(cstring chemin_)
{
  impl->charge(chemin_);
}

float WavLecteur::lis_duree() const
{
  si(!impl->fd)
    retourne 0;
  retourne ((float) impl->nech) / impl->entete.frequency;
}

WavConfig WavLecteur::lis_config() const
{
  retourne {impl->entete.nbr_canaux, (double) impl->entete.frequency};
}

string WavLecteur::lis_infos() const
{
  si(!impl->fd)
    retourne "(aucun fichier wav ouvert)";
  retourne format("{}, {} cana{}, {} Hz, {:.1f} sec.",
      impl->chemin,
      (entier) impl->entete.nbr_canaux,
      impl->entete.nbr_canaux > 1 ? "ux" : "l",
          (entier) impl->entete.frequency,
          ((float) impl->nech) / impl->entete.frequency);
}

entier WavLecteur::lis_nechantillons() const
{
  si(!impl->fd)
    retourne 0;
  retourne impl->nech;
}

bouléen WavLecteur::eof() const
{
  retourne (!impl->fd) || (feof(impl->fd));
}

entier WavLecteur::lis_position() const
{
  retourne impl->pos;
}

float WavLecteur::lis_position_secondes() const
{
  retourne ((float) impl->pos) / impl->entete.frequency;
}

/*ArrayXXcf WavLecteur::lis_cplx()
{
  retourne impl->lis_cplx();
}*/

Vecf WavLecteur::lis_flottant(entier n)
{
  retourne impl->lis_flottant(n);
}

entier WavLecteur::lis_donnees(int16_t *res, entier n)
{
  retourne impl->lis_donnees(res, n);
}

void WavLecteur::redemarre()
{
  si(!impl->fd)
    échec("wav-lecteur/redemarre: non ouvert.");
  msg("Remise a zero lecteur wav.");
  fseek(impl->fd, sizeof(impl->entete), SEEK_SET);
  impl->pos = 0;
}

void WavLecteur::recherche(float temps)
{
  si(!impl->fd)
    échec("non ouvert.");

  impl->pos = (entier) (temps * impl->entete.frequency);

  msg("Recherche, pos = {} secondes ({} echan / {} total) ({} %)", temps, impl->pos, impl->nech, ((float) impl->pos) / impl->nech);

  fseek(impl->fd, impl->pos * 2 + sizeof(impl->entete), SEEK_SET);
}

struct WavEcrivain::Impl
{
  string chemin;
  WavConfig cfg;
  FILE *fd = nullptr;
  uint32_t pos = 0;

  ~Impl()
  {
    ferme();
  }

  void init(cstring chemin, const WavConfig &config)
  {
    wave_entete_t entete;

    ferme();

    this->chemin = chemin;

    msg("Ouverture fichier wav en ecriture [{}], ncanaux = {}, fech = {:.1f} kHz.",
          chemin, config.ncanaux, config.freq_ech * 1e-3f);

    /*auto ext = cutils::fichiers::extrait_extension(chemin);
    si(ext != "wav")
    {
      avertissement("Ceci n'est pas un fichier WAV (extension = %s).", ext);
      retourne ERREUR_WAV_EXTENSION;
    }*/

    fd = fopen(chemin.c_str(), "wb");
    si(!fd)
      échec("Echec ouverture fichier wav en ecriture ({}). errno={} ({}).", chemin, errno, strerror(errno));
    // génération d'un entete vide, à remplir à la fermeture du fichier
    memset(&entete, 0, sizeof(wave_entete_t));
    fwrite(&entete, sizeof(wave_entete_t), 1, fd);

    /* copie de la config passï¿½e */
    this->cfg = config;

    pos = 0;
  }

  void ferme()
  {
    si(fd)
    {
      msg("Fermeture fichier wav, nech = {}...", pos);

      bouléen flottant = cfg.format == WavConfig::PCM_FLOAT;
      entier octets_par_echantillon = flottant ? 4 : 2;

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
  }

  void ecris(const Veccf &x)
  {
    //ArrayXf x2 = Eigen::Map<ArrayXf>((float *) x.data(), x.rows() * 2);
    // TODO: no copy
    entier n = x.rows();
    Vecf x2(2*n);
    memcpy(x2.rawptr(), x.rawptr(), 2 * n * sizeof(float));
    ecris(x2);
  }

  void ecris(const Vecf &x)
  {
    //infos("Ecriture %d ech float...", x.rows());
    si(fd)
    {
      si(cfg.format == WavConfig::PCM_16)
      {
        soit p = (int16_t *) malloc(2 * x.rows());
        pour(auto i = 0; i < x.rows(); i++)
          p[i] = (int16_t) x(i);
        ecris(p, x.rows());
        free(p);
      }
      sinon
      {
        fwrite(x.data(), 4, x.rows(), fd);
        pos += x.rows() / cfg.ncanaux;
      }
    }
  }

  void ecris(int16_t *donnees, entier n)
  {
    si(fd)
    {
      si(cfg.format == WavConfig::PCM_16)
      {
        fwrite(donnees, 2, n, fd);
        pos += n / cfg.ncanaux;
      }
      sinon
      {
        Vecf x(n);
        pour(auto i = 0; i < n; i++)
          x(i) = donnees[i];
        ecris(x);
      }
    }
  }

};


WavEcrivain::WavEcrivain()
{
  impl = make_shared<Impl>();
}

string WavEcrivain::lis_infos() const
{
  float duree = ((float) impl->pos) / impl->cfg.freq_ech;
  retourne impl->chemin + format(" ({} cana{}, {:.1f} kHz, {:.1f} secondes)",
      impl->cfg.ncanaux,
      impl->cfg.ncanaux > 1 ? "ux" : "l",
          (float) (impl->cfg.freq_ech * 1e-3), duree);
}

void WavEcrivain::init(cstring chemin, const WavConfig &config)
{
  impl->init(chemin, config);
}

float WavEcrivain::lis_position_secondes() const
{
  retourne ((float) impl->pos) / impl->cfg.freq_ech;
}

void WavEcrivain::ecris(int16_t *donnees, entier n)
{
  impl->ecris(donnees, n);
}

void WavEcrivain::ecris(const Vecf &x)
{
  impl->ecris(x);
}

void WavEcrivain::ecris(const Veccf &x)
{
  impl->ecris(x);
}

void WavEcrivain::ferme()
{
  impl->ferme();
}




}

