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

using namespace std;

namespace tsd::vue {

struct FreeTypeFont: Font
{
  bouléen init_ok = non;
  FT_Library library;
  FT_Face face;
  float echelle = 1;

  mutex mut;

  struct CarSpec
  {
    uint32_t charcode;
    entier echelle;
    strong_ordering operator<=>(const CarSpec&) const = default;
  };

  struct CarRendu
  {
    Image img;
    entier inc_x = 0, inc_y = 0;
    entier bitmap_left = 0, bitmap_top = 0;
  };

  map<CarSpec, CarRendu> cache;

  ~FreeTypeFont()
  {
    FT_Done_Face(face);
    FT_Done_FreeType(library);
  }

  FreeTypeFont()
  {
    init();
  }

  void init()
  {
    const lock_guard<mutex> lock(mut);
    init_ok = non;
    si(FT_Init_FreeType(&library))
      échec("FT_Init_FreeType");

    string nom = "OpenSans-Regular.ttf";

    vector<string> cds =
      {
          "./data/fonts",             // Installeur win32
          "c:/msys64/usr/share/tsd",  // Installeur win32
          "/usr/share/tsd",           // Installeur linux
          "./build/debug/data/fonts",  // Mode dév
          "./build/debug-linux/data/fonts"  // Mode dév
      };

    string fn;
    pour(auto &s: cds)
    {
      si(filesystem::exists(s + "/" + nom))
      {
        fn = s + "/" + nom;
        break;
      }
    }

    si(fn.empty())
      échec("Fichier de fonte non trouvé ({})", nom);

    si(FT_New_Face(library, fn.c_str(), 0, &face))
      échec("FT_New_Face (fichier non trouvé : [{}]).", fn);

    FT_Select_Charmap(face , ft_encoding_unicode);

    /* use 50pt at 100dpi */
    /* set character size
     *  width, heigth */
    entier res = FT_Set_Char_Size(face, 50 * 64, 0, 50, 0);
    si(res)
      échec("FT_Set_Char_Size : {:x}", res);
    echelle = 1;
    init_ok = oui;
  }

  CarRendu rendre_car_cache(const CarSpec &sp)
  {
    si(cache.find(sp) == cache.end())
      cache[sp] = rendre_car(sp);
    retourne cache[sp];
  }

  CarRendu rendre_car(const CarSpec &sp)
  {
    CarRendu res;
    entier gl_index = FT_Get_Char_Index(face, sp.charcode);
    VERBOSE(msg("load glyph...");)
    // C'est ça qu'il faudrait "cacher"
    soit error = FT_Load_Glyph(face, gl_index, FT_LOAD_RENDER);
    VERBOSE(msg("ok.");)
    /* ignore errors */
    si(error)
      retourne res;

    soit slot = face->glyph;
    soit &bmp = slot->bitmap;

    entier sx = bmp.width, sy = bmp.rows;

    //si(/*(s[i] != ' ')*/ /*&&*/ (sx * sy == 0))
    //{
      //msg_avert("Caractère 0x{:02x} : sx = {}, sy = {}.", s[i], sx, sy);
    //}

    si(sx * sy > 0)
    {
      VERBOSE(msg("dessin gl...");)
      res.img.resize(sx, sy);

      //img.remplir(Couleur{255,255,255,0});
      soit ptr  = (int32_t *) res.img.data();
      pour(auto y = 0; y < sy; y++)
      {
        pour(auto x = 0; x < sx; x++)
        {
          int32_t val = (bmp.buffer[y * sx + x] & 0xff);
          *ptr++ = 0 | (0 << 8) | (0 << 16) | (val << 24);
        }
      }

      res.bitmap_left = slot->bitmap_left;
      res.bitmap_top = slot->bitmap_top;
      VERBOSE(msg("ok.");)
    }

    res.inc_x = slot->advance.x / 64;
    res.inc_y = slot->advance.y / 64;

    retourne res;
  }

  Image rendre(cstring s, float echelle)
  {
    const lock_guard<mutex> lock(mut);

    si(!init_ok || (s.size() == 0))
      retourne Image();

    si(echelle != this->echelle)
    {
      this->echelle = echelle;
      //msg("nv echelle : set char size...");
      entier res = FT_Set_Char_Size(face, 50 * 64, 0, echelle * 50, 0);
      //msg("ok");
      si(res)
      {
        msg_erreur("FT_Set_Char_Size : 0x{:x}", res);
        retourne Image();
      }
    }

    entier n = s.size();
    vector<Image> imgs;
    vector<Point> poss;
    Point pos{0,0};
    Dim dim{0,0};
    soit top_max = 0;

    //msg("Rendu freetype : [{}]...", s);

    pour(auto i = 0; i < n; i++)
    {
      uint32_t charcode;

      // UTF8 sur un octet
      si((s[i] & (1 << 7)) == 0)
      {
        charcode = s[i];
      }
      // UTF8 sur deux octets
      sinon si((s[i] & (1 << 5)) == 0)
      {
        charcode = (s[i] & 0b11111);
        si(i+1 >= n)
        {
          msg_erreur("PB UTF8");
          break;
        }
        charcode = (charcode << 6) | (s[i+1] & 0b111111);
        i++;
      }
      // UTF8 sur trois octets
      sinon si((s[i] & (1 << 4)) == 0)
      {
        si(i+2 >= n)
        {
          msg_erreur("PB UTF8");
          break;
        }
        charcode = (s[i] & 0b1111);
        charcode = (charcode << 6) | (s[i+1] & 0b111111);
        charcode = (charcode << 6) | (s[i+2] & 0b111111);
        i += 2;
      }
      // UTF8 sur quatre octets
      sinon //si((s[i] & (1 << 2)) == 0)
      {
        si(i+3 >= n)
        {
          msg_erreur("PB UTF8");
          break;
        }
        charcode = (s[i] & 0b11);
        charcode = (charcode << 6) | (s[i+1] & 0b111111);
        charcode = (charcode << 6) | (s[i+2] & 0b111111);
        charcode = (charcode << 6) | (s[i+3] & 0b111111);
        i += 3;
      }
#     if 0
      sinon
      {
        si(i+3 >= n)
        {
          msg_erreur("PB UTF8");
          break;
        }
        msg_erreur("TODO : UTF8 4 octets / cas spécial : s = {:b}.{:b}.{:b}.{:b}.",
            (unsigned char) s[i], (unsigned char) s[i+1], (unsigned char) s[i+2], (unsigned char) s[i+3]);
        break;
        /*
        charcode = (s[i] & 0b11);
        charcode = (charcode << 6) | (s[i+1] & 0b111111);
        charcode = (charcode << 6) | (s[i+2] & 0b111111);
        charcode = (charcode << 6) | (s[i+3] & 0b111111);*/
      }
      /*sinon
      {
        msg_erreur("TODO : UTF8 4 octets.");
        break;
      }*/
#     endif

      CarSpec sp;
      sp.charcode = charcode;
      sp.echelle  = echelle * 50;
      soit rendu = rendre_car_cache(sp);

      si(!rendu.img.empty())
      {
        Point pos2{pos.x + rendu.bitmap_left, rendu.bitmap_top};
        dim.l   = max(pos2.x + rendu.img.sx(), dim.l);
        top_max = max(top_max, rendu.bitmap_top);
        imgs.push_back(rendu.img);
        poss.push_back(pos2);
      }

      pos.x += rendu.inc_x;
      pos.y += rendu.inc_y;
    }

    pour(auto i = 0u; i < imgs.size(); i++)
    {
      poss[i].y = top_max - poss[i].y;
      dim.h = max(dim.h, poss[i].y + imgs[i].sy());
    }

    VERBOSE(msg("assemblage...");)
    Image img(dim.l, dim.h);
    img.remplir(Couleur{255,255,255,0});
    si(img.empty())
    {
      msg_avert("Rendu freetype : dim totale = {} x {} (echelle = {}), texte = '{}', dim texte = {} caractères.", dim.l, dim.h, echelle, s, s.size());
      msg("imgs ({}) :", imgs.size());
      pour(auto &i: imgs)
        msg("  {}", i.get_dim());
      retourne img;
    }

    pour(auto i = 0u; i < imgs.size(); i++)
    {
      assertion(poss[i].x + imgs[i].sx() <= dim.l);
      assertion(poss[i].y + imgs[i].sy() <= dim.h);
      img.puti(poss[i], imgs[i]);
    }
    VERBOSE(msg("ok.");)
    retourne img;
  }

};

sptr<Font> fonte_ft_creation()
{
  retourne make_shared<FreeTypeFont>();
}
}





