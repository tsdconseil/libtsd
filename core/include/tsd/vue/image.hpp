#pragma once

/** @cond undoc */

#include "tsd/tsd.hpp"
#include <vector>
#include <cstdint>

/** Interface abstraite vers fonction de dessin */

namespace tsd::vue {

struct Couleur
{
  unsigned char r = 0, g = 0, b = 0, alpha = 255;

  Couleur eclaircir(float ratio) const;
  Couleur assombrir(float ratio) const;
  Couleur(int32_t rgba);
  int32_t vers_rgba() const;
  unsigned char lumi() const;
  Couleur(){}
  Couleur(const std::string &s);
  //Couleur(unsigned char r, unsigned char g, unsigned char b, unsigned char alpha = 255);
  Couleur(float r, float g, float b, float alpha = 255);

  std::string vers_chaine() const;

  static const Couleur Blanc, Noir, Rouge, Vert, Bleu, Violet, Jaune, Orange, Cyan, Marron, Gris;

  std::strong_ordering operator <=>(const Couleur &c) const = default;
};

extern std::ostream& operator<<(std::ostream& ss, const Couleur &t);

template<typename T>
struct Point_
{
  T x = 0, y = 0;

  Point_(){}

  Point_(T x, T y){this->x = x; this->y = y;}

  template<typename T2>
    Point_(const Point_<T2> &r)
    {
      x = r.x;
      y = r.y;
    }

  std::strong_ordering operator<=>(const Point_<T>&) const = default;
};

template<typename T>
  Point_<T> operator +(const Point_<T> &p0, const Point_<T> &p1)
{
  return Point_<T>{p0.x+p1.x, p0.y+p1.y};
}

template<typename T>
  Point_<T> operator -(const Point_<T> &p0, const Point_<T> &p1)
{
  return Point_<T>{p0.x-p1.x, p0.y-p1.y};
}

using Point = Point_<int>;
using Pointf = Point_<float>;

template<typename T>
struct Dim_
{
  T l = 0, h = 0;
  std::partial_ordering operator<=>(const Dim_<T>&) const = default;
};

using Dim = Dim_<int>;
using Dimf = Dim_<float>;



template<typename T>
struct Rect_
{
  T x = 0, y = 0, l = 0, h = 0;

  Rect_(){}
  Rect_(T x, T y, T l, T h)
  {
    this->x = x;
    this->y = y;
    this->l = l;
    this->h = h;
  }
  Rect_(const Point_<T> &p0, const Point_<T> &p1)
  {
    *this = englobant(p0, p1);
  }

  Rect_(const Point_<T> &p, const Dim_<T> &d)
  {
    x = p.x;
    y = p.y;
    l = d.l;
    h = d.h;
  }
  template<typename T2>
    Rect_(const Rect_<T2> &r)
    {
      x = r.x;
      y = r.y;
      l = r.l;
      h = r.h;
    }

  Point_<T> tl() const
  {
    return {x, y};
  }

  Point_<T> centre() const
  {
    return {x + l/2, y + h/2};
  }

  Point_<T> br() const
  {
    if constexpr (std::is_integral<T>::value)
      return {x+l-1, y+h-1};
    else
      return {x+l, y+h};
  }

  static Rect_<T> englobant(const Point_<T> &p0, const Point_<T> &p1)
  {
    if constexpr (std::is_integral<T>::value)
      return Rect_<T>{std::min(p0.x,p1.x), std::min(p0.y,p1.y), std::abs(p1.x-p0.x)+1, std::abs(p1.y-p0.y)+1};
    else
      return Rect_<T>{std::min(p0.x,p1.x), std::min(p0.y,p1.y), std::abs(p1.x-p0.x), std::abs(p1.y-p0.y)};
  }

  Dim_<T> dim() const
  {
    return {l, h};
  }

  bool contient(const Point_<T> &p) const
  {
    return (p.x >= x) && (p.y >= y) && (p.x <= x + l) && (p.y <= y + h);
  }

  std::strong_ordering operator<=>(const Rect_<T>&) const = default;
};


using Rect = Rect_<int>;
using Rectf = Rect_<float>;


template<typename T>
  std::ostream& operator<<(std::ostream& s, const Rect_<T>& t)
{
  s << fmt::format("rect({}x{}, {}x{})", t.x, t.y, t.l, t.h);
  return s;
}

template<typename T>
  std::ostream& operator<<(std::ostream& s, const Point_<T>& t)
{
  s << fmt::format("{}x{}", t.x, t.y);
  return s;
}

template<typename T>
  std::ostream& operator<<(std::ostream& s, const Dim_<T>& t)
{
  s << fmt::format("{}x{}", t.l, t.h);
  return s;
}

enum class Orientation
{
  HORIZONTALE,
  VERTICALE
};


struct CMap
{
  virtual void calc(float t, float &r, float &v, float &b) = 0;
  Couleur couleur(float t);
  std::string nom;
};

extern sptr<CMap> cmap_parse(const std::string &nom);

// Classe g??n??rique pour le dessin
struct Image
{
  Image(int sx = 0, int sy = 0, void *data = nullptr);
  Image(const Dim &dim, void *data = nullptr) : Image(dim.l, dim.h, data){}

  const void *data() const;
  void *data();

  Image clone() const;

  void remplir(const Couleur &c);

  void resize(const Dim &d);
  void resize(int sx, int sy);

  Image redim(int l, int h) const;
  Image redim(const Dim &dim) const {return redim(dim.l, dim.h);}

  void def_couleur_dessin(const Couleur &c);
  void def_couleur_remplissage(const Couleur &c);
  void def_epaisseur(int ep);

  enum StyleLigne
  {
    PLEINE, POINTILLEE, TIRETS
  };

  void ligne_aa(float x0, float y0, float x1, float y1, StyleLigne style = PLEINE);
  void ligne(const Point &p0, const Point &p1, StyleLigne style = PLEINE);


  void fleche(const Point &p0, const Point &p1, StyleLigne style = PLEINE, float lg = 5);

  void point(const Point &p0);
  void point(const Point &p0, const Couleur &c);
  void point(const Point &p0, const Couleur &c, float alpha);
  void rectangle(const Point &p0, const Point &p1);
  void rectangle(const Rect &r);
  void ellipse(const Point &p0, const Point &p1);
  void cercle(const Point &p0, int r);

  void rectangle_plein(const Point &p0, const Point &p1);
  void ellipse_pleine(const Point &p0, const Point &p1);
  void cercle_plein(const Point &p0, int r);

  Dim texte_dim(const std::string &s, float dim);

  enum Alignement
  {
    DEBUT,
    MILIEU,
    FIN
  };

  void puts(const Point &p0, const std::string &s, float dim, Alignement align_horizontal = Alignement::DEBUT, Alignement align_vertical = Alignement::DEBUT);

  void puti(const Point &pos, Image src, Rect rdi_source);
  void puti(const Rect &rdi, Image src);
  void puti(const Point &pos, Image src);
  void puti_avec_gamma(const Point &pos, Image src);
  void puti(const Point &pos, Image src, float ??);
  // Pose uniquement le gamma, avec la couleur en cours
  void put_gamma(const Point &pos, Image src);
  Image rotation_90() const;

  // Blend, in-place
  void blend(Image src);
  Image blend_nv(Image src);

  void enregister(const std::string &chemin);
  void charger(const std::string &chemin);
  Image sous_image(int x0, int y0, int l, int h);

  Dim get_dim() const;
  inline int sx() const{return get_dim().l;}
  inline int sy() const{return get_dim().h;}

  inline bool empty() const{return sx() * sy() == 0;}

private:
  struct Impl;
  sptr<Impl> impl;
};

struct Font
{
  virtual ~Font(){}
  /** Sur l'image r??sultante, seul le gamma est significatif
   *  (toutes les valeurs RVB sont ?? z??ro).
   */
  virtual Image rendre(const std::string &s, float scale = 1.0f) = 0;
};



extern sptr<Font> fonte_ft_creation();


struct TexteConfiguration
{
  Couleur couleur      = Couleur{0,0,0};
  Couleur couleur_fond = Couleur{255,255,255,0};
  float scale             = 0.5;
  int thickness           = 1;
  int fontface            = 1;
  Dim dim_max             = {-1,-1};
  Point org               = {0,0};
  float transparence      = 0.3;
  enum
  {
    ALIGN_GAUCHE = 0,
    ALIGN_CENTRE,
    ALIGN_DROIT
  } alignement = ALIGN_GAUCHE;
};

struct TexteProps
{
  float scale_out;
  std::vector<int> xdim, ypos, ydim;
};

// Affiche du texte, ??ventuellement sur plusieurs lignes
// L'image de sortie est cr????e.
extern Image texte_creation_image(const std::string &s,
                                        const TexteConfiguration &config,
                                        TexteProps *props = nullptr);


// Affiche du texte dans un cadre semi-transparent au dessus d'une image existante.
extern void texte_ajoute(Image O, const TexteConfiguration &config,
    const std::string &s, ...);


}


/** @endcond */
