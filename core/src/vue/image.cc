#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include <vector>
#include <cstdarg>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;


namespace tsd::vue {

static bouléen est_deci(char c)
{
  retourne ((c >= '0') && (c <= '9'));
}

// D'après libcutil
vector<entier> parse_liste_entiers(cstring str)
{
  soit n = (entier) str.size();
  vector<entier> res;
  entier current = 0;
  si(!est_deci(str[0]))
    retourne res;
  pour(auto i = 0; i < n; i++)
  {
    si(est_deci(str[i]))
      current = current * 10 + (str[i] - '0');
    sinon
    {
      soit c = str[i];
      si(((c == '.') || (c == ',') || (c == '/') || ((c == ':'))) && (i + 1 < n))
      {
        res.push_back(current);
        current = 0;
      }
      sinon
        retourne vector<entier>();
    }
  }
  res.push_back(current);
  retourne res;
}

Couleur Couleur::mélange(const Couleur &a, const Couleur &b, float α)
{
  retourne
  {
    α * a.r + (1 - α) * b.r,
    α * a.g + (1 - α) * b.g,
    α * a.b + (1 - α) * b.b
  };
}

Couleur Couleur::rand()
{
  soit r = randi(256, 3);
  retourne Couleur{(float) r(0), (float) r(1), (float) r(2)};
}

Couleur::Couleur(float r, float g, float b, float alpha)
{
  this->r     = clamp(r, 0.0f, 255.0f);
  this->g     = clamp(g, 0.0f, 255.0f);
  this->b     = clamp(b, 0.0f, 255.0f);
  this->alpha = clamp(alpha, 0.0f, 255.0f);
}

ostream& operator<<(ostream& ss, const Couleur &t)
{
  ss << sformat("(r={},v={},b={},a={})", t.r, t.g, t.b, t.alpha);
  retourne ss;
}

string Couleur::vers_chaine() const
{
  retourne sformat("{}.{}.{}", r, g, b);
}

Couleur::Couleur(cstring s)
{
  alpha = 255;
  si(!s.empty())
  {
    soit lst = parse_liste_entiers(s);
    si(lst.size() == 3)
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
  retourne (alpha << 24) | (r << 16) | (g << 8) | b;
}

unsigned char Couleur::lumi() const
{
  retourne 0.299 * r + 0.587 * g + 0.114 * b;
}

const Couleur
  Couleur::Blanc        {255,255,255},
  Couleur::Noir         {0,0,0},
  Couleur::Rouge        {255,0,0},
  Couleur::Vert         {0,255,0},
  Couleur::Bleu         {0,0,255},
  Couleur::Violet       {128,0,128},
  Couleur::Orange       {255,165,0},
  Couleur::Jaune        {255,255,0},
  Couleur::Cyan         {0,255,255},
  Couleur::Marron       {139,69,19},
  Couleur::Gris         {128,128,128},
  Couleur::BleuSombre   {0,0,180},
  Couleur::VertSombre   {0,100,0},
  Couleur::RougeSombre  {210,0,0},
  Couleur::CyanSombre   {0,192,192},
  Couleur::VioletSombre {192,0,192},
  Couleur::JauneSombre  {192,192,0},
  Couleur::MarronSombre {128,70,0},
  Couleur::OrangeSombre {210,118,0};


  Couleur Couleur::eclaircir(float α) const
  {
    retourne {
        255 - (255 - r) * α,
        255 - (255 - g) * α,
        255 - (255 - b) * α};
  }

  Couleur Couleur::assombrir(float α) const
  {
    retourne
    {
      r * α,
      g * α,
      b * α
    };
  }






Image affiche_texte(
    const vector<string> &lignes,
    const TexteConfiguration &config,
    TexteProps *props = nullptr)
{
  Image O(1, 1);
  O.remplir(config.couleur_fond);

  // On applique d'abord les paramètres demandés, et on regarde si tout rentre
  // Sinon, on diminue le scale, et on recommence, jusqu'à ce que tout rentre.
  soit échelle = config.échelle;
  entier nl = lignes.size();

  vector<entier> ypos(nl), xdim(nl), ydim(nl);

  si(props)
  {
    props->échelle_appliquée = 0;
    props->xdim.clear();
    props->ypos.clear();
  }

  Dim dim_max = config.dim_max;

  si(dim_max.l == -1)
    dim_max = {1000, 800};

  si((dim_max.l < 1) || (dim_max.h < 1))
  {
    //msg_avert("Pas possible de placer du texte dans un rectangle de dimension ({},{}) pixels", dim_max.l, dim_max.h);
    retourne O;
  }

  Dim dim_totale;
  soit nitr = 0;
  pour(;;)
  {
    nitr++;
    dim_totale = {0,0};
    pour(auto i = 0u; i < lignes.size(); i++)
    {
      ypos[i] = dim_totale.h;

      soit dim = O.texte_dim(lignes[i], échelle);

      si(i + 1 < lignes.size())
        dim.h += 5; // Ajoute un peu d'espace entre les lignes

      dim_totale.l = max(dim.l, dim_totale.l);

      ydim[i] = dim.h;
      xdim[i] = dim.l;
      dim_totale.h += dim.h;
    }
    si(config.dim_max.l == -1)
      break;
    si((dim_totale.l <= dim_max.l) && (dim_totale.h <= dim_max.h))
      break;
    si(échelle < 1e-3)
    {
      //msg_avert("Image::affiche_texte({} lignes, [{}]). Echec placement @scale = {}. Dim max = {}.", lignes.size(), lignes[0], échelle, dim_max);
      retourne O;
    }
    // Essai en plus petit
    échelle *= 0.9;
  }

  //msg("affiche_texte(): {} itr, échelle réelle appliquée: {}, dim out={}", nitr, échelle, dim_totale);

  si(props)
  {
    props->échelle_appliquée = échelle;
    props->xdim              = xdim;
    props->ydim              = ydim;
    props->ypos              = ypos;
  }

  O.resize(dim_totale);
  O.remplir(config.couleur_fond);
  O.def_couleur_dessin(config.couleur);
  pour(auto i = 0u; i < lignes.size(); i++)
  {
    soit x = 0;
    si(config.alignement == TexteConfiguration::ALIGN_CENTRE)
      x = dim_totale.l/2 - xdim[i]/2;

    Dim d = O.texte_dim(lignes[i], échelle);
    si((d.l + x > O.sx()) || (d.h + ypos[i] > O.sy()))
    {
      msg("dim texte : {}, dim O : {}, x={}, ypos[{}] = {}", d, O.get_dim(), x, i, ypos[i]);
      msg_erreur("Dépassement texte.");
      break;
    }
    //msg("dim texte : {}", d);

    O.puts({x, ypos[i]}, lignes[i], échelle);
  }
  //msg("affiche texte : end.");
  retourne O;
}


Image texte_creation_image(cstring s,
                           const TexteConfiguration &config,
                           TexteProps *props)
{
  // Séparation des différentes lignes du texte
  vector<string> lignes;
  soit i = 0;

  pour(;;)
  {
    soit pos_rl = s.find('\n', i);
    si(pos_rl == string::npos)
    {
      soit l = s.substr(i);
      lignes.push_back(l);
      break;
    }
    soit l = s.substr(i, pos_rl - i);
    lignes.push_back(l);
    i = pos_rl + 1;
  }

  retourne affiche_texte(lignes, config, props);
}
}

