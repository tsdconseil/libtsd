#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include <vector>
#include <cstdarg>
#include <cassert>
#include <iostream>
#include <algorithm>

namespace tsd::vue {

  static bool est_deci(char c)
  {
    return ((c >= '0') && (c <= '9'));
  }

  // D'après libcutil
  std::vector<int> parse_liste_entiers(const std::string &str)
  {
    std::vector<int> res;
    auto s = str.c_str();
    res.clear();
    int current = 0;
    if(!est_deci(s[0]))
      return std::vector<int>();
    while(strlen(s) > 0)
    {
      if(est_deci(s[0]))
      {
        current = current * 10 + (s[0] - '0');
      }
      else
      {
        auto c = s[0];
        if(((c == '.') || (c == ',') || (c == '/') || ((c == ':'))) && (strlen(s) > 1))
        {
          res.push_back(current);
          current = 0;
        }
        else
        {
          return std::vector<int>();
        }
      }
      s++;
    }
    res.push_back(current);
    return res;
  }

  Couleur::Couleur(float r, float g, float b, float alpha)
  {
    this->r = std::clamp(r, 0.0f, 255.0f);
    this->g = std::clamp(g, 0.0f, 255.0f);
    this->b = std::clamp(b, 0.0f, 255.0f);
    this->alpha = std::clamp(alpha, 0.0f, 255.0f);
  }

  /*Couleur::Couleur(unsigned char r, unsigned char g, unsigned char b, unsigned char alpha)
  {
    this->r = r;
    this->g = g;
    this->b = b;
    this->alpha = alpha;
  }*/

std::ostream& operator<<(std::ostream& ss, const Couleur &t)
{
  ss << fmt::format("(r={},v={},b={},a={})", t.r, t.g, t.b, t.alpha);
  return ss;
}

std::string Couleur::vers_chaine() const
{
  return fmt::format("{}.{}.{}", r, g, b);
}

Couleur::Couleur(const std::string &s)
{
  alpha = 255;
  if(!s.empty())
  {
    auto lst = parse_liste_entiers(s);
    if(lst.size() == 3)
    {
      r = lst[0];
      g = lst[1];
      b = lst[2];
    }
  }
}

Couleur::Couleur(int32_t rgba)
{
  alpha = (rgba >> 24) & 0xff;
  r     = (rgba >> 16) & 0xff;
  g     = (rgba >> 8) & 0xff;
  b     = rgba & 0xff;
}

int32_t Couleur::vers_rgba() const
{
  return (alpha << 24) | (r << 16) | (g << 8) | b;
}

unsigned char Couleur::lumi() const
{
  return 0.299 * r + 0.587 * g + 0.114 * b;
}

const Couleur
  Couleur::Blanc  {255,255,255,255},
  Couleur::Noir   {0,0,0,255},
  Couleur::Rouge  {255,0,0,255},
  Couleur::Vert   {0,255,0,255},
  Couleur::Bleu   {0,0,255,255},
  Couleur::Violet {128,0,128,255},
  Couleur::Orange {255,165,0,255},
  Couleur::Jaune  {255,255,0,255},
  Couleur::Cyan   {0,255,255,255},
  Couleur::Marron {139,69,19,255},
  Couleur::Gris   {128,128,128,255};


  Couleur Couleur::eclaircir(float ratio) const
  {
    Couleur y;
    y.r = 255 - (255 - r) * ratio;
    y.g = 255 - (255 - g) * ratio;
    y.b = 255 - (255 - b) * ratio;
    return y;
  }

  Couleur Couleur::assombrir(float ratio) const
  {
    Couleur y;
    y.r = r * ratio;
    y.g = g * ratio;
    y.b = b * ratio;
    return y;
  }






  Image affiche_texte(
                   const std::vector<std::string> &lignes,
                   const TexteConfiguration &config,
                   TexteProps *props = nullptr)
{
    Image O(1, 1);
    O.remplir(config.couleur_fond);

  // On applique d'abord les paramètres demandés, et on regarde si tout rentre
  // Sinon, on diminue le scale, et on recommence, jusqu'à ce que tout rentre.

  float scale = config.scale;
  int nl = lignes.size();

  std::vector<int> ypos(nl), xdim(nl), ydim(nl);

  if(props != nullptr)
  {
    props->scale_out = 0;
    props->xdim.clear();
    props->ypos.clear();
  }

  Dim dim_max = config.dim_max;

  if(dim_max.l == -1)
  {
    dim_max.l = 1000;
    dim_max.h = 800;
  }

  if((dim_max.l < 1) || (dim_max.h < 1))
  {
    //msg_avert("Pas possible de placer du texte dans un rectangle de dimension ({},{}) pixels", dim_max.l, dim_max.h);
    return O;
  }


  //msg("affiche texte : start.");
  Dim dim_totale;
  int nitr = 0;
  for(;;)
  {
    nitr++;
    dim_totale = {0,0};
    for(auto i = 0u; i < lignes.size(); i++)
    {
      ypos[i] = dim_totale.h;
      Dim dim = O.texte_dim(lignes[i], scale);

      dim.h += 5; // Ajoute un peu d'espace entre les lignes

      dim_totale.l = std::max(dim.l, dim_totale.l);

      ydim[i] = dim.h;
      xdim[i] = dim.l;
      dim_totale.h += dim.h;
    }
    if(config.dim_max.l == -1)
      break;
    if((dim_totale.l <= dim_max.l) && (dim_totale.h <= dim_max.h))
      break;
    if(scale < 1e-3)
    {
      msg_avert("Image::affiche_texte({} lignes, [{}]). Echec placement @scale = {}. Dim max = {}.", lignes.size(), lignes[0], scale, dim_max);
      return O;
    }
    scale *= 0.9;
  }
  //msg("affiche texte : {} itr.", nitr);
  //infos("dim_totale = %d, %d", dim_totale.width, dim_totale.height);
  //dim_out = dim_totale;

  if(props != nullptr)
  {
    props->scale_out = scale;
    props->xdim = xdim;
    props->ydim = ydim;
    props->ypos = ypos;
  }

  // A supprimer
  //dim_totale.l += 10;
  //dim_totale.h += 10;

  /*dim_totale.width  += 5;
  dim_totale.height += 2;*/

  O.resize(dim_totale);
  O.remplir(config.couleur_fond);
  O.def_couleur_dessin(config.couleur);
  for(auto i = 0u; i < lignes.size(); i++)
  {
    //msg("ypos[{}] = {}, dessin {}", i, ypos[i], lignes[i]);
    auto x = 0;//config.org.x;
    if(config.alignement == TexteConfiguration::ALIGN_CENTRE)
    {
      x = dim_totale.l/2 - xdim[i]/2;
    }

    Dim d = O.texte_dim(lignes[i], scale);
    if((d.l + x > O.sx()) || (d.h + ypos[i] > O.sy()))
    {
      msg_erreur("Dépassement texte.");
      msg("dim texte : {}, dim O : {}, x={}, ypos[{}] = {}",
          d, O.get_dim(), x, i, ypos[i]);
      break;
    }
    //infos("dim texte : %d x %d", d.l, d.h);

    O.puts({x, ypos[i]}, lignes[i], scale);
  }
  //msg("affiche texte : end.");
  return O;
}


void texte_ajoute(Image O, const TexteConfiguration &config, const std::string &s, ...)
{
  va_list ap;
  va_start(ap, s);
  char buf[5000];
  vsnprintf(buf, 5000, s.c_str(), ap);
  va_end(ap);

  TexteConfiguration tc = config;
  if(tc.dim_max.l == -1)
    tc.dim_max.l = O.sx() - tc.org.x;
  tc.dim_max.l = std::min(tc.dim_max.l, O.sx() - tc.org.x);
  if(tc.dim_max.h == -1)
    tc.dim_max.h = O.sy() - tc.org.y;
  tc.dim_max.h = std::min(tc.dim_max.h, O.sy() - tc.org.y);

  auto T = texte_creation_image(buf, tc, nullptr);

  //infos("O.cols = %d, tc.org = %d, T.cols = %d", O.cols, tc.org.x, T.cols);
  tsd_assert(T.sx() + tc.org.x <= O.sx());
  tsd_assert(T.sy() + tc.org.y <= O.sy());

  O.puti(Point{tc.org.x, tc.org.y}, T);
  //rdi_cible = tc.transparence * rdi_cible + (1.0 - tc.transparence) * T;
}


Image texte_creation_image(const std::string &s,
    const TexteConfiguration &config,
    TexteProps *props)
{

  //infos("Text creation image[%s] : scale = %f", s.c_str(), config.scale);

  // Séparation des différentes lignes du texte
  std::vector<std::string> lignes;
  unsigned int i = 0;

  //infos("affiche texte [%s]", s.c_str());

  for(;;)
  {
    auto pos_rl = s.find('\n', i);
    if(pos_rl == std::string::npos)
    {
      auto l = s.substr(i);
      lignes.push_back(l);
      break;
    }
    auto l = s.substr(i, pos_rl - i);
    lignes.push_back(l);
    i = pos_rl + 1;
  }
  //for(auto i = 0u; i < lignes.size(); i++)
    //infos("Ligne[%d] : [%s]", i, lignes[i].c_str());

  return affiche_texte(lignes, config, props);
}
}

