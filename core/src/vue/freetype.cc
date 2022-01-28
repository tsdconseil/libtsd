#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"

#include <cstdio>
#include <string>
#include <cmath>
#include <cassert>
#include <thread>
#include <mutex>
#include <map>
#include <filesystem>

extern "C"
{
# include <ft2build.h>
# include FT_FREETYPE_H
}

#define VERBOSE(AA)


namespace tsd::vue {

struct FreeTypeFont: Font
{
  bool init_ok = false;
  FT_Library library;
  FT_Face face;
  float echelle = 1;

  std::mutex mutex;

  struct CarSpec
  {
    uint32_t charcode;
    int echelle;
    std::strong_ordering operator<=>(const CarSpec&) const = default;
  };

  struct CarRendu
  {
    Image img;
    int inc_x = 0, inc_y = 0;
    int bitmap_left = 0, bitmap_top = 0;
  };

  std::map<CarSpec, CarRendu> cache;

  ~FreeTypeFont()
  {
    FT_Done_Face(face);
    FT_Done_FreeType(library);
  }

  FreeTypeFont()
  {
    init();
  }

  int init()
  {
    const std::lock_guard<std::mutex> lock(mutex);
    init_ok = false;
    if(FT_Init_FreeType(&library))
    {
      msg_erreur("FT_Init_FreeType");
      return -1;
    }

    /* create face object */
    //std::string fn = "./build/debug/data/fonts/goudy_bookletter_1911.otf";
    //std::string fn = "./build/debug/data/fonts/OpenSans-Light.ttf";
    // std::string fn = "./build/debug/data/fonts/OpenSans-Regular.ttf";
    /*
#   ifdef WIN
    // TODO
    std::string fn = "C:/msys64/usr/share/tsd/OpenSans-Regular.ttf";
#   else
    std::string fn = "/usr/share/tsd/OpenSans-Regular.ttf";
#   endif
  */

    std::string nom = "OpenSans-Regular.ttf";

    std::vector<std::string> cds =
      {
          "./data/fonts",             // Installeur win32
          "/usr/share/tsd",           // Installeur linux
          "./build/debug/data/fonts",  // Mode dév
          "./build/debug-linux/data/fonts"  // Mode dév
      };

    std::string fn;
    for(auto &s: cds)
    {
      if(std::filesystem::exists(s + "/" + nom))
      {
        fn = s + "/" + nom;
        break;
      }
    }

    if(fn.empty())
    {
      msg_erreur("Fichier de fonte non trouvé ({})", nom);
      return -1;
    }

    if(FT_New_Face(library, fn.c_str(), 0, &face))
    {
      msg_erreur("FT_New_Face (fichier non trouvé : [{}]).", fn);
      return -1;
    }

    /* use 50pt at 100dpi */
    /* set character size
     *  width, heigth */
    int res = FT_Set_Char_Size(face, 50 * 64, 0, 50, 0);
    if(res)
    {
      msg_erreur("FT_Set_Char_Size : {:x}", res);
      return -1;
    }
    echelle = 1;
    init_ok = true;
    return 0;
  }

  CarRendu rendre_car_cache(const CarSpec &sp)
  {
    if(cache.find(sp) == cache.end())
      cache[sp] = rendre_car(sp);
    return cache[sp];
  }

  CarRendu rendre_car(const CarSpec &sp)
  {
    CarRendu res;
    int gl_index = FT_Get_Char_Index(face, sp.charcode);
    VERBOSE(msg("load glyph...");)
    // C'est ça qu'il faudrait "cacher"
    auto error = FT_Load_Glyph(face, gl_index, FT_LOAD_RENDER);
    VERBOSE(msg("ok.");)
    /* ignore errors */
    if(error)
      return res;

    auto slot = face->glyph;
    auto &bmp = slot->bitmap;

    int sx = bmp.width, sy = bmp.rows;

    //if(/*(s[i] != ' ')*/ /*&&*/ (sx * sy == 0))
    //{
      //msg_avert("Caractère 0x{:02x} : sx = {}, sy = {}.", s[i], sx, sy);
    //}

    if(sx * sy > 0)
    {
      VERBOSE(msg("dessin gl...");)
      res.img.resize(sx, sy);

      //img.remplir(Couleur{255,255,255,0});
      auto ptr  = (int32_t *) res.img.data();
      for(auto y = 0; y < sy; y++)
      {
        for(auto x = 0; x < sx; x++)
        {
          int32_t val = (bmp.buffer[y * sx + x] & 0xff);
          *ptr++ = 0 | (0 << 8) | (0 << 16) | (val << 24);
        }
      }

      res.bitmap_left = slot->bitmap_left;
      res.bitmap_top = slot->bitmap_top;
      //Point pos2{pos.x + (int) slot->bitmap_left, slot->bitmap_top};
      //dim.l   = std::max(pos2.x + sx, dim.l);
      //top_max = std::max(top_max, (int) slot->bitmap_top);
      //imgs.push_back(img);
      //poss.push_back(pos2);
      VERBOSE(msg("ok.");)
    }

    res.inc_x = slot->advance.x / 64;
    res.inc_y = slot->advance.y / 64;

    return res;
  }

  Image rendre(const std::string &s, float echelle)
  {
    const std::lock_guard<std::mutex> lock(mutex);

    if(!init_ok || (s.size() == 0))
      return Image();

    if(echelle != this->echelle)
    {
      this->echelle = echelle;
      //msg("nv echelle : set char size...");
      int res = FT_Set_Char_Size(face, 50 * 64, 0, echelle * 50, 0);
      //msg("ok");
      if(res)
      {
        msg_erreur("FT_Set_Char_Size : 0x{:x}", res);
        return Image();
      }
    }

    int n = s.size();
    std::vector<Image> imgs;
    std::vector<Point> poss;
    Point pos{0,0};
    Dim dim{0,0};
    int top_max = 0;

    //msg("Rendu freetype : [{}]...", s);

    for(auto i = 0; i < n; i++)
    {
      uint32_t charcode;

      if((s[i] & (1 << 7)) == 0)
      {
        charcode = s[i];
      }
      else if((s[i] & (1 << 5)) == 0)
      {
        charcode = (s[i] & 0b11111);
        if(i+1 >= n)
        {
          msg_erreur("PB UTF8");
          break;
        }
        charcode = (charcode << 6) | (s[i+1] & 0b111111);
        i++;
      }
      else
      {
        msg_erreur("TODO : UTF8 3 / 4 octets.");
        break;
      }

      CarSpec sp;
      sp.charcode = charcode;
      sp.echelle  = echelle * 50;
      auto rendu = rendre_car_cache(sp);

      if(!rendu.img.empty())
      {
        Point pos2{pos.x + rendu.bitmap_left, rendu.bitmap_top};
        dim.l   = std::max(pos2.x + rendu.img.sx(), dim.l);
        top_max = std::max(top_max, rendu.bitmap_top);
        imgs.push_back(rendu.img);
        poss.push_back(pos2);
      }

      pos.x += rendu.inc_x;
      pos.y += rendu.inc_y;

#     if 0
      int gl_index = FT_Get_Char_Index(face, charcode);
      VERBOSE(msg("load glyph...");)
      // C'est ça qu'il faudrait "cacher"
      auto error = FT_Load_Glyph(face, gl_index, FT_LOAD_RENDER);
      VERBOSE(msg("ok.");)
      /* ignore errors */
      if(error)
        continue;

      auto slot = face->glyph;
      auto &bmp = slot->bitmap;

      int sx = bmp.width, sy = bmp.rows;

      if((s[i] != ' ') && (sx * sy == 0))
      {
        msg_avert("Caractère 0x{:02x} : sx = {}, sy = {}.", s[i], sx, sy);
      }

      if(sx * sy > 0)
      {
        VERBOSE(msg("dessin gl...");)
        Image img(sx, sy);

        //img.remplir(Couleur{255,255,255,0});
        auto ptr  = (int32_t *) img.data();
        for(auto y = 0; y < sy; y++)
        {
          for(auto x = 0; x < sx; x++)
          {
            int32_t val = (bmp.buffer[y * sx + x] & 0xff);
            *ptr++ = 0 | (0 << 8) | (0 << 16) | (val << 24);
          }
        }

        Point pos2{pos.x + (int) slot->bitmap_left, slot->bitmap_top};
        dim.l   = std::max(pos2.x + sx, dim.l);
        top_max = std::max(top_max, (int) slot->bitmap_top);
        imgs.push_back(img);
        poss.push_back(pos2);
        VERBOSE(msg("ok.");)
      }

      pos.x += slot->advance.x / 64;
      pos.y += slot->advance.y / 64;
#     endif
    }

    for(auto i = 0u; i < imgs.size(); i++)
    {
      poss[i].y = top_max - poss[i].y;
      dim.h = std::max(dim.h, poss[i].y + imgs[i].sy());
    }

    VERBOSE(msg("assemblage...");)
    Image img(dim.l, dim.h);
    img.remplir(Couleur{255,255,255,0});
    if(img.empty())
    {
      msg_avert("Rendu freetype : dim totale = {} x {} (echelle = {}), texte = '{}', dim texte = {} caractères.", dim.l, dim.h, echelle, s, s.size());
      msg("imgs ({}) :", imgs.size());
      for(auto &i: imgs)
        msg("  {}", i.get_dim());
      return img;
    }

    for(auto i = 0u; i < imgs.size(); i++)
    {
      tsd_assert(poss[i].x + imgs[i].sx() <= dim.l);
      tsd_assert(poss[i].y + imgs[i].sy() <= dim.h);
      img.puti(poss[i], imgs[i]);
    }
    VERBOSE(msg("ok.");)
    return img;
  }

};

sptr<Font> fonte_ft_creation()
{
  return std::make_shared<FreeTypeFont>();
}
}





