//#pragma once
#ifndef IMAGE_HPP
#define IMAGE_HPP

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

  static Couleur mélange(const Couleur &a, const Couleur &b, float alpha);
  static Couleur rand();

  std::string vers_chaine() const;

  static const Couleur Blanc, Noir, Rouge, Vert, Bleu, Violet, Jaune, Orange, Cyan, Marron, Gris;

  static const Couleur BleuSombre, VertSombre, RougeSombre, CyanSombre, VioletSombre, JauneSombre, MarronSombre, OrangeSombre;

  std::strong_ordering operator <=>(const Couleur &c) const = default;
};

extern std::ostream& operator<<(std::ostream& ss, const Couleur &t);

template<typename T>
struct Point_
{
  T x = 0, y = 0;

  Point_(){}

  template<typename T2, typename T3>
  Point_(T2 x, T3 y){this->x = (T) x; this->y = (T) y;}

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

using Point = Point_<entier>;
using Pointf = Point_<float>;

template<typename T>
struct Dim_
{
  T l = 0, h = 0;

  Dim_(){}

  template<typename T1, typename T2>
  Dim_(T1 l, T2 h)
  {
    this->l = (T) l;
    this->h = (T) h;
  }

  std::partial_ordering operator<=>(const Dim_<T>&) const = default;
};

using Dim = Dim_<entier>;
using Dimf = Dim_<float>;



template<typename T>
struct Rect_
{
  T x = 0, y = 0, l = 0, h = 0;

  Rect_(){}
  template<typename T1, typename T2, typename T3, typename T4>
  Rect_(T1 x, T2 y, T3 l, T4 h)
  {
    this->x = (T) x;
    this->y = (T) y;
    this->l = (T) l;
    this->h = (T) h;
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
      return Rect_<T>{min(p0.x,p1.x), min(p0.y,p1.y), abs(p1.x-p0.x)+1, abs(p1.y-p0.y)+1};
    else
      return Rect_<T>{min(p0.x,p1.x), min(p0.y,p1.y), abs(p1.x-p0.x), abs(p1.y-p0.y)};
  }

  Dim_<T> dim() const
  {
    return {l, h};
  }

  bouléen contient(const Point_<T> &p) const
  {
    return (p.x >= min(x, x + l)) && (p.y >= min(y, y + h)) && (p.x <= max(x, x + l)) && (p.y <= max(y, y + h));
  }

  std::strong_ordering operator<=>(const Rect_<T>&) const = default;
};


using Rect = Rect_<entier>;
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

// Classe générique pour le dessin
struct Image
{
  Image(entier sx = 0, entier sy = 0, void *data = nullptr);
  Image(const Dim &dim, void *data = nullptr) : Image(dim.l, dim.h, data){}

  const void *data() const;
  void *data();

  Image clone() const;

  void remplir(const Couleur &c);

  void resize(const Dim &d);
  void resize(entier sx, entier sy);

  Image redim(entier l, entier h) const;
  Image redim(const Dim &dim) const {return redim(dim.l, dim.h);}

  void def_couleur_dessin(const Couleur &c);
  void def_couleur_remplissage(const Couleur &c);
  void def_epaisseur(entier ep);

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
  void cercle(const Point &p0, entier r);

  void rectangle_plein(const Point &p0, const Point &p1);
  void ellipse_pleine(const Point &p0, const Point &p1);
  void cercle_plein(const Point &p0, entier r);

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
  void puti(const Point &pos, Image src, float γ);
  // Pose uniquement le gamma, avec la couleur en cours
  void put_gamma(const Point &pos, Image src);
  Image rotation_90() const;

  // Blend, in-place
  void blend(Image src);
  Image blend_nv(Image src);

  void enregister(const std::string &chemin) const;
  void charger(const std::string &chemin);
  Image sous_image(entier x0, entier y0, entier l, entier h);

  Dim get_dim() const;
  inline entier sx() const{return get_dim().l;}
  inline entier sy() const{return get_dim().h;}

  inline bouléen empty() const{return sx() * sy() == 0;}

private:
  struct Impl;
  sptr<Impl> impl;
};

struct Font
{
  virtual ~Font(){}
  /** Sur l'image résultante, seul le gamma est significatif
   *  (toutes les valeurs RVB sont à zéro).
   */
  virtual Image rendre(const std::string &s, float scale = 1.0f) = 0;
};



extern sptr<Font> fonte_ft_creation();


struct TexteConfiguration
{
  Couleur couleur      = Couleur{0,0,0};
  Couleur couleur_fond = Couleur{255,255,255,0};
  float scale             = 0.5;
  entier thickness           = 1;
  entier fontface            = 1;
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
  std::vector<entier> xdim, ypos, ydim;
};

// Affiche du texte, éventuellement sur plusieurs lignes
// L'image de sortie est créée.
extern Image texte_creation_image(const std::string &s,
                                        const TexteConfiguration &config,
                                        TexteProps *props = nullptr);


// Affiche du texte dans un cadre semi-transparent au dessus d'une image existante.
extern void texte_ajoute(Image O, const TexteConfiguration &config,
    const std::string &s, ...);


}


ostream_formater(tsd::vue::Point)
ostream_formater(tsd::vue::Pointf)
ostream_formater(tsd::vue::Rect)
ostream_formater(tsd::vue::Rectf)
ostream_formater(tsd::vue::Dim)
ostream_formater(tsd::vue::Dimf)
ostream_formater(tsd::vue::Couleur)

#if 0
template <> struct fmt::formatter<tsd::vue::Point> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Point& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "{}x{}", t.x, t.y);
  }
};

template <> struct fmt::formatter<tsd::vue::Pointf> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Pointf& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "{}x{}", t.x, t.y);
  }
};


template <> struct fmt::formatter<tsd::vue::Rect> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Rect& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "rect({}x{}, {}x{})", t.x, t.y, t.l, t.h);
  }
};

template <> struct fmt::formatter<tsd::vue::Rectf> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Rectf& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "rect({}x{}, {}x{})", t.x, t.y, t.l, t.h);
  }
};


template <> struct fmt::formatter<tsd::vue::Dim> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Dim& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "{}x{}", t.l, t.h);
  }
};

template <> struct fmt::formatter<tsd::vue::Dimf> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Dimf& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "{}x{}", t.l, t.h);
  }
};

template <> struct fmt::formatter<tsd::vue::Couleur> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::vue::Couleur& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "(r={},v={},b={},a={})", t.r, t.g, t.b, t.alpha);
  }
};
#endif

/** @endcond */
#endif
