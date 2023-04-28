#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <cassert>
#include <utility>

#if LIBTSD_USE_PNG
# include <png.h>
#endif

#ifdef WIN
# include "windows.h"
#else
# include <sys/stat.h>
#endif


#define DBG_AXES(AA)

using namespace std;

namespace tsd::vue {

struct bgra
{
  unsigned char bgra[4];
};


// [chemin, fichier]
static tuple<string, string>  analyse_chemin(cstring chemin_complet)
{
  si(chemin_complet == ".")
    retourne {".", ""};

  soit n = (entier) chemin_complet.size();
  // Index dernier '/'
  soit j = -1;


  pour(auto i = 0; i < n; i++)
    si((chemin_complet[i] == '/') || (chemin_complet[i] == '\\'))
      j = i;

  // Pas de '/'
  si(j == -1)
    retourne {"", chemin_complet};

  // '/' exclu
  retourne {chemin_complet.substr(0, j), chemin_complet.substr(j+1)};
}

static bouléen fichier_existe(const string &nom)
{
  soit f = fopen(nom.c_str(), "r");
  si(f == nullptr)
    retourne non;
  fclose(f);
  retourne oui;
}

static void creation_dossier(const string &chemin)
{
# ifdef WIN
  CreateDirectory(chemin.c_str(), nullptr);
# else
  entier res = mkdir(chemin.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  si(!res)
    msg_avert("Echec création dossier {}.", chemin);
# endif
}

static void creation_dossier_si_inexistant(const string &chemin)
{
  si(!fichier_existe(chemin))
    creation_dossier(chemin);
}



struct Image::Impl
{
  entier sx = 0, sy = 0;
  bgra *bitmap = nullptr;
  bgra cd = {0}, cf = {0};
  entier ep = 1;
  static sptr<Font> ftfonte;


  inline float fpart(float x)
  {
    si(x > 0)
      retourne x - (entier) x;
    sinon
      retourne x - ((entier) (x)+1);
  }

  Dim dim() const
  {
    retourne {sx, sy};
  }

  void charger(const string &chemin)
  {

#   if LIBTSD_USE_PNG
    //msg("lecture png [{}]...", chemin);
    entier width, height;
    png_byte color_type;
    png_byte bit_depth;
    png_bytep *row_pointers = NULL;

    FILE *fp = fopen(chemin.c_str(), "rb");

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    si(!png)
      echec("png_create_read_struct");

    png_infop info = png_create_info_struct(png);
    si(!info)
      echec("png_create_info_struct");

    si(setjmp(png_jmpbuf(png)))
      echec("png_jmpbuf(png)");

    png_init_io(png, fp);

    png_read_info(png, info);

    width      = png_get_image_width(png, info);
    height     = png_get_image_height(png, info);
    color_type = png_get_color_type(png, info);
    bit_depth  = png_get_bit_depth(png, info);

    // Read any color_type into 8bit depth, RGBA format.
    // See http://www.libpng.org/pub/png/libpng-manual.txt

    si(bit_depth == 16)
      png_set_strip_16(png);

    si(color_type == PNG_COLOR_TYPE_PALETTE)
      png_set_palette_to_rgb(png);

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    si(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
      png_set_expand_gray_1_2_4_to_8(png);

    si(png_get_valid(png, info, PNG_INFO_tRNS))
      png_set_tRNS_to_alpha(png);

    // These color_type don't have an alpha channel then fill it with 0xff.
    si(color_type == PNG_COLOR_TYPE_RGB ||
       color_type == PNG_COLOR_TYPE_GRAY ||
       color_type == PNG_COLOR_TYPE_PALETTE)
      png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    si(color_type == PNG_COLOR_TYPE_GRAY ||
       color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
      png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    //soit nrb = png_get_rowbytes(png,info);

    //msg("mallocs (h = {}, w = {}, nrb = {})...", height, width, nrb);

    row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    pour(entier y = 0; y < height; y++) {
      row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
    }

    msg("Read image...");
    png_read_image(png, row_pointers);

    fclose(fp);

    msg("sx : {}, 4*sx: {}", sx, 4 * sx);

    //msg("Destroy read struct...");
    png_destroy_read_struct(&png, &info, NULL);

    msg("Resize({}, {})..", width, height);
    resize(width, height);

    msg("lecture octets..");
    pour(auto y = 0; y < height; y++)
    {
      soit iptr = (unsigned char *) row_pointers[y];
      soit optr = ((unsigned char *) ((char *) this->bitmap)) + y * sx * 4;
      pour(auto x = 0; x < width; x++)
      {
        optr[2] = *iptr++;
        optr[1] = *iptr++;
        optr[0] = *iptr++;
        optr[3] = 255;
        optr += 4;
        iptr++;
      }

      //memcpy(((char *) this->bitmap) + y * sx * 4, row_pointers[y], 4 * sx);
      free(row_pointers[y]);
    }
    free(row_pointers);
    msg("ok.");

#   endif
  }


  void enregister(const string &chemin) const
  {
#   if LIBTSD_USE_PNG

    string ch = chemin;

    si((chemin.size() < 4) || (chemin[chemin.size()-4] != '.'))
      ch += ".png";

    soit [dossier, fichier] = analyse_chemin(ch);
    creation_dossier_si_inexistant(dossier);

    //si(ch.ends_with(".png"))
    si((ch.size() >= 4) && (ch.substr(ch.size() - 4, 4) == ".png"))
    {
      entier width = sx, height = sy;
      png_byte color_type = PNG_COLOR_TYPE_RGBA;
      png_byte bit_depth = 8;

      png_structp png_ptr;
      png_infop info_ptr;
      //entier number_of_passes;
      png_bytep * row_pointers = (png_bytep *) malloc(sy * sizeof(png_bytep *));
      soit ptr = (unsigned char *) data();

      soit ptri = (unsigned char *) malloc(sx * sy * 4);
      soit tmp = ptri;
      pour(auto i = 0; i < sx * sy; i++)
      {
        *tmp++ = ptr[2];
        *tmp++ = ptr[1];
        *tmp++ = ptr[0];
        *tmp++ = 255;

        ptr += 4;
      }

      pour(auto y = 0; y < sy; y++)
        row_pointers[y] = ptri + y * sx * 4;


      /* create file */
      FILE *fp = fopen(ch.c_str(), "wb");
      si (!fp)
      {
        msg_erreur("[write_png_file] File {} could not be opened pour writing", ch.c_str());
        retourne;
      }


      /* initialize stuff */
      png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

      si (!png_ptr)
      {
        msg_erreur("[write_png_file] png_create_write_struct failed");
        retourne;
      }

      info_ptr = png_create_info_struct(png_ptr);
      si (!info_ptr)
      {
        msg_erreur("[write_png_file] png_create_info_struct failed");
        retourne;
      }

      si (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during init_io");
        retourne;
      }

      png_init_io(png_ptr, fp);


      /* write header */
      si(setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during writing header");
        retourne;
      }

      png_set_IHDR(png_ptr, info_ptr, width, height,
                   bit_depth, color_type, PNG_INTERLACE_NONE,
                   PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

      png_write_info(png_ptr, info_ptr);


      /* write bytes */
      si (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during writing bytes");
        retourne;
      }

      png_write_image(png_ptr, row_pointers);


      /* end write */
      si (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during end of write");
        retourne;
      }

      png_write_end(png_ptr, NULL);

      /* cleanup heap allocation */
      //pour(y=0; y<height; y++)
        // free(row_pointers[y]);
      free(row_pointers);
      free(ptri);

      fclose(fp);
      retourne;
    }
#   endif

    msg_erreur("Enregistrement fichier image ({}) : pas de codec configuré pour cette extension.", chemin);
  }

  // On suppose que le signal d'entrée a été interpolé
  // pour des coordonnées entières de x
  void courbe_xy(entier xmin, const Vecf &Y)
  {
    entier xmax = Y.rows() + xmin - 1;
    pour(auto x = xmin; x < xmax; x++)
    {
      // Plutôt vertical, il faut utiliser les pixels adjacents horizontaux
      si(abs(Y(x+1)-Y(x)) > 1)
      {

      }
      sinon
      {

      }



      /*pour(auto y = ymin; y <= ymax; y++)
      {

      }*/
    }
  }


  void pta(entier x, entier y, float a)
  {
    si((x < 0) || (y < 0) || (x >= sx) || (y >= sy))
      retourne;

    soit ptr = bitmap + x + y * sx;
    soit ip = (bgra *) &cd;
    float a0 = sqrt(a) * cd.bgra[3] / 255.0f;
    float b0 = 1 - a0;

    pour(auto k = 0u; k < 4; k++)
      ptr->bgra[k] = ptr->bgra[k] * b0 + ip->bgra[k] * a0;
  }

  bouléen est_entier(float x) const
  {
    retourne abs(x - floor(x)) < 1e-6;
  }

  // https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
  // (p1,p2) <-> (p3,p4)
  static Pointf intersecte_lignes(
      const Pointf &p1, const Pointf &p2,
      const Pointf &p3, const Pointf &p4)
  {
    soit d = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x),
         a = p1.x*p2.y-p1.y*p2.x,
         b = p3.x*p4.y - p3.y*p4.x;
    retourne
    {
      (a * (p3.x - p4.x) - (p1.x - p2.x) * b) / d,
      (a * (p3.y - p4.y) - (p1.y - p2.y) * b) / d
    };
  }

  // Première intersection en partant de p0
  tuple<Pointf, int> intersecte_ligne_rdi(const Pointf &p0, const Pointf &p1) const
  {
    Pointf r{-1,-1};
    int côté = -1;

    // Intersection avec ligne du haut
    soit pnord  = intersecte_lignes(p0, p1, {0,0}, {sx-1,0}),
         psud   = intersecte_lignes(p0, p1, {0,sy-1}, {sx-1,sy-1}),
         pest   = intersecte_lignes(p0, p1, {sx-1,0}, {sx-1,sy-1}),
         pouest = intersecte_lignes(p0, p1, {0,0}, {0,sy-1});

    soit t = [&](const Pointf &p)
    {
      si((p.x >= 0) && (p.x <= sx-1) && (p.y >= 0) && (p.y <= sy-1))
        r = p;
    };

    si((p0.y < 0) && (p0.x >= 0) && (p0.x <= sx-1))
      t(pnord);
    sinon si((p0.y >= sy) && (p0.x >= 0) && (p0.x <= sx-1))
      t(psud);
    sinon si((p0.x < 0) && (p0.y >= 0) && (p0.y <= sy-1))
      t(pouest);
    sinon si((p0.x >= sx-1) && (p0.y >= 0) && (p0.y <= sy-1))
      t(pest);
    sinon si((p0.x <= 0) && (p0.y <= 0))
    {
      t(pouest);
      t(pnord);
    }
    sinon si((p0.x >= sx-1) && (p0.y <= 0))
    {
      t(pest);
      t(pnord);
    }
    sinon si((p0.x <= 0) && (p0.y >= sy-1))
    {
      t(pouest);
      t(psud);
    }
    sinon si((p0.x >= sx-1) && (p0.y >= sy-1))
    {
      t(pest);
      t(psud);
    }

    retourne {r, côté};
  }

  // retourne vrai si en dehors
  bouléen restreint_segment_rdi(float &x0, float &y0, float &x1, float &y1) const
  {
    // Trouve les points d'intersection avec la zone d'intérêt,
    // et ne trace que dans cet intervalle

    Pointf p0{x0,y0}, p1{x1,y1};

    soit i0 = est_dans_rdi(p0);
    soit i1 = est_dans_rdi(p1);

    si(i0 && i1)
      retourne non;

    int c0 = -1, c1 = -1;

    si(!i0)
    {
      std::tie(p0,c0) = intersecte_ligne_rdi(p0, p1);
    }
    si(!i1)
    {
      std::tie(p1,c1) = intersecte_ligne_rdi(p1, p0);
    }

    si((p0.x == -1) || (p1.x == -1))
      retourne oui;

    si(!i0 && !i1)
    {
      si(c0 == c1)
        retourne oui;
    }

    x0 = p0.x;
    y0 = p0.y;
    x1 = p1.x;
    y1 = p1.y;

    retourne non;
  }


  void ligne_aa(float x0, float y0, float x1, float y1, StyleLigne style)
  {
    si(restreint_segment_rdi(x0, y0, x1, y1))
     retourne;

    si(((x0 == x1) || (y0 == y1))
        && est_entier(x0) && est_entier(y0)
        && est_entier(y1) && est_entier(x1))
      ligne({(entier) x0, (entier) y0}, {(entier) x1, (entier) y1}, style);
    sinon
    {
      pour(auto k = -ep/2; k <= (ep-1)/2; k++)
      {
        si(abs(x1 - x0) > abs(y1 - y0))
          ligne_aa_int(x0, y0 + k, x1, y1 + k, style);
        sinon
          ligne_aa_int(x0 + k, y0, x1 + k, y1, style);
      }
    }
  }

  void ligna_aa_nv(float x0, float y0, float x1, float y1)
  {

  }

  void ligne_aa_int(float x0, float y0, float x1, float y1, StyleLigne style = PLEINE)
  {
    //infos("Ligne aa : (%dx%d) -> (%dx%d)", x0, y0, x1, y1);



    bouléen plutot_verticale = abs(y1 - y0) > abs(x1 - x0);

    si(plutot_verticale)
    {
      swap(x0 , y0);
      swap(x1 , y1);
    }
    si(x0 > x1)
    {
      swap(x0,x1);
      swap(y0,y1);
    }

    //compute the slope
    float dx = x1 - x0;
    float dy = y1 - y0;
    float pas = dx == 0 ? 1 : dy / dx;
    float y = y0;

    //Couleur c(cd);

    entier xpxl1, xpxl2;

    {
      // handle first endpoint
      entier xend = (entier) round(x0);
      float yend = y0 + pas * (xend - x0);
      float xgap = 1 - fpart(x0 + 0.5);
      xpxl1 = xend; // this will be used in the main loop
      entier ypxl1 = (entier) (yend);
      si(plutot_verticale)
      {
        pta(ypxl1,   xpxl1, (1 - fpart(yend)) * xgap);
        pta(ypxl1+1, xpxl1, fpart(yend) * xgap);
      }
      sinon
      {
        pta(xpxl1, ypxl1  , (1 - fpart(yend)) * xgap);
        pta(xpxl1, ypxl1+1,  fpart(yend) * xgap);
      }
      y = yend + pas;
    }



    //intery := yend + gradient // first y-intersection pour the main loop

    {
      // handle second endpoint
      entier xend = round(x1);
      float yend = y1 + pas * (xend - x1);
      float xgap = fpart(x1 + 0.5);
      xpxl2 = xend; //this will be used in the main loop
      entier ypxl2 = (entier) (yend);
      si(plutot_verticale)
      {
        pta(ypxl2  , xpxl2, (1-fpart(yend)) * xgap);
        pta(ypxl2+1, xpxl2,  fpart(yend) * xgap);
      }
      sinon
      {
        pta(xpxl2, ypxl2,  (1-fpart(yend)) * xgap);
        pta(xpxl2, ypxl2+1, fpart(yend) * xgap);
      }
    }


    entier inc = 1;
    si(style == POINTILLEE)
    {
      inc = 3;
    }

    // main loop
    si(plutot_verticale)
    {
      pour(auto x = xpxl1+1; x <= xpxl2-1; x += inc)
      {
        float r = fpart(y);
        pta((entier) y, x, 1-r);
        pta((entier) y+1, x, r);
        y += pas * inc;
      }
    }
    sinon
    {
      pour(auto x = xpxl1+1; x <= xpxl2-1; x += inc)
      {
        float r = fpart(y);
        pta(x, (entier) y, 1-r);
        pta(x, (entier) y+1, r);
        y += pas * inc;
      }
    }

  }

  ~Impl()
  {
    si(bitmap != nullptr)
      free(bitmap);
    bitmap = nullptr;
    sx = sy = 0;
  }

  Impl(entier sx, entier sy, void *data = nullptr)
  {
    this->sx = sx;
    this->sy = sy;
    si(sx * sy > 0)
    {
      bitmap = (bgra *) malloc(4*((size_t)sx)*((size_t)sy));

      si(bitmap == nullptr)
        echec("Echec allocation image : {} x {} ({:.1f} Mo)", sx, sy, 4.0f*sx*sy*1e-6);

      si(data != nullptr)
        memcpy(bitmap, data, sx * sy * 4);
      //sinon
        //memset(bitmap, 255, sx * sy * 4);
    }

    si(!ftfonte)
    {
      ftfonte = fonte_ft_creation();
    }
  }

  Image clone() const
  {
    Image res(sx, sy);
    memcpy(res.impl->bitmap, bitmap, sx * sy * 4);
    retourne res;
  }

  const void *data() const
  {
    retourne bitmap;
  }
  void *data()
  {
    retourne bitmap;
  }

  Image rotation_90() const
  {
    Image res(sy, sx);

    soit iptr = (const int32_t *) bitmap;
    soit optr = (int32_t *) res.impl->bitmap;

    optr += sx * sy - 1;
    pour(auto x = 0; x < sx; x++)
      pour(auto y = 0; y < sy; y++)
        *optr-- = iptr[(sx - 1 - x) + y * sx];

    retourne res;
  }




  inline bgra &pixel(entier x, entier y)
  {
    retourne *(bitmap + y * sx + x);
  }

  inline const bgra &pixel(entier x, entier y) const
  {
    retourne *(bitmap + y * sx + x);
  }


  // Pose uniquement le gamma, avec la couleur en cours
  void put_gamma(const Point &pos, Image src)
  {
    soit isx = src.impl->sx, isy = src.impl->sy;

    si((isx + pos.x > sx) || (isy + pos.y > sy) || (pos.x < 0) || (pos.y < 0))
    {
      //msg_avert("ImageBmp:: put gamma : dépassement.");
      retourne;
    }

    pour(auto y = 0; y < isy; y++)
    {
      pour(auto x = 0; x < isx; x++)
      {
        const soit &iv = src.impl->pixel(x, y);
        soit &ov = pixel(x + pos.x, y + pos.y);
        soit a    = iv.bgra[3] / 255.0f,
             a_or = ov.bgra[3] / 255.0f;

        pour(auto k = 0; k < 3; k++)
          ov.bgra[k]      = (1 - a) * ov.bgra[k] + a * cd.bgra[k];

        ov.bgra[3]  = max(iv.bgra[3], ov.bgra[3]);

        si((iv.bgra[3] > 0) && (a_or == 0))
          ov.bgra[3]  = 255 - (255 - iv.bgra[3]) / 2;
      }
    }
  }

  void puti_avec_gamma(const Point &pos, Image src)
  {
    soit isx = src.impl->sx, isy = src.impl->sy;

    si(pos.x < 0)
    {
      puti_avec_gamma({0, pos.y}, src.sous_image(-pos.x, 0, src.sx() + pos.x, src.sy()));
      retourne;
    }
    si(pos.y < 0)
    {
      puti_avec_gamma({pos.x, 0}, src.sous_image(0, -pos.y, src.sx(), src.sy() + pos.y));
      retourne;
    }

    si((isx + pos.x > sx) || (isy + pos.y > sy) || (pos.x < 0) || (pos.y < 0))
    {
      msg_avert("ImageBmp:: put image avec gamma : dépassement : dim cible = {}x{}, pos={}, dim src={},{}.",
          sx, sy, pos, isx, isy);
      retourne;
    }

    si((pos.x >= 2e9) || (pos.y >= 2e9))
    {
      msg_erreur("pos = {}", pos);
    }

    pour(auto y = 0; y < isy; y++)
    {
      pour(auto x = 0; x < isx; x++)
      {
        // Enlever cette assertion !

        si(!(((x + pos.x >= 0) && (y + pos.y >= 0) && (x + pos.x < sx) && (y + pos.y < sy))))
        {
          echec("x={},y={},isx={},isy={},sx={},sy={},p.x={},p.y={},isy + pos.y={}",x,y,isx,isy,sx,sy,pos.x,pos.y,isy + pos.y);
        }

        tsd_assert((x + pos.x >= 0) && (y + pos.y >= 0) && (x + pos.x < sx) && (y + pos.y < sy));

        soit iv  = src.impl->pixel(x, y);
        soit a   = iv.bgra[3] / 255.0f;
        soit &ov = pixel(x + pos.x, y + pos.y);

        pour(auto k = 0; k < 3; k++)
          ov.bgra[k] = (1 - a) * ov.bgra[k] + a * iv.bgra[k];
        ov.bgra[3]  = iv.bgra[3];
      }
    }
  }

  inline bgra *ptr(entier x, entier y)
  {
    retourne bitmap + x + y * sx;
  }

  Image sous_image(entier x0, entier y0, entier l, entier h)
  {
    si((x0 < 0) || (y0 < 0) || (x0 + l > sx) || (y0 + h > sy) || (l < 0) || (h < 0))
    {
      msg_avert("sous_image : hors borne.");
      retourne Image();
    }

    Image res(l, h);
    pour(auto y = 0; y < h; y++)
      memcpy(res.impl->bitmap + y * l, this->bitmap + x0 + (y0 + y) * sx, l * 4);
    retourne res;
  }

  void puti(const Point &pos, Image src, Rect rdi_source)
  {
    pour(auto y = 0; y < rdi_source.h; y++)
      memcpy(ptr(pos.x, pos.y + y), src.impl->ptr(rdi_source.x, rdi_source.y + y), rdi_source.l * 4);
  }

  void puti(const Point &pos, Image src)
  {
    si((pos.x < 0) || (pos.y < 0) || (pos.x + src.impl->sx > sx) || (pos.y + src.impl->sy > sy))
    {
      DBG_AXES(msg("puti : pos = {}", pos);)
      Point pos2 = pos;

      soit x0 = max(-pos.x, 0), y0 = max(-pos.y, 0);

      si(pos2.x < 0)
        pos2.x = 0;
      si(pos2.y < 0)
        pos2.y = 0;

      soit l = src.impl->sx - x0, h = src.impl->sy - y0;

      si(l + pos2.x > sx)
        l = sx - pos2.x;
      si(h + pos2.y > sy)
        h = sy - pos2.y;

      puti(pos2, src.sous_image(x0, y0, l, h));
      retourne;
    }

    soit dimx = src.impl->sx;
    si(dimx + pos.x > sx)
      dimx = sx - pos.x;

    pour(auto y = pos.y; (y < pos.y + src.impl->sy) && (y < sy); y++)
      si(dimx > 0)
        memcpy(ptr(pos.x, y), src.impl->ptr(0, y - pos.y), dimx * 4);
  }

  void puti(const Rect &rdi, Image src)
  {
    si((rdi.l == src.impl->sx) && (rdi.h == src.impl->sy))
    {
      puti(Point{rdi.x, rdi.y}, src);
      retourne;
    }

    si((rdi.l <= 0) || (rdi.h <= 0))
    {
      msg_avert("RDI = {} x {} !!!", rdi.l, rdi.h);
      retourne;
    }

    // Image vide
    si(src.sx() * src.sy() == 0)
      retourne;

    msg_avert("Attention: redim {} -> {}x{}", src.get_dim(), rdi.l, rdi.h);
    puti(rdi, src.redim(rdi.l, rdi.h));
  }

  void puti(const Point &pos, Image src, float γ)
  {

  }



  Image blend_nv(Image src)
  {
    soit s1 = this;
    soit s2 = src.impl;
    tsd_assert(s2);
    si((s1->sx != s2->sx) || (s1->sy != s2->sy))
    {
      msg_erreur("ImageBmp::blend() : image source de dim différente : ({}) != ({})", s1->dim(), src.get_dim());
      retourne Image();
    }

    Image res(s1->sx, s2->sy);

    soit optr  = (bgra *) res.impl->bitmap;
    soit iptr1 = (const bgra *) s1->bitmap,
         iptr2 = (const bgra *) s2->bitmap;

    pour(auto i = 0; i < sy*sx; i++)
    {
      soit a = (*iptr2).bgra[3] / 255.0f;
      pour(auto k = 0; k < 3; k++)
      {
        soit &ov  = (*optr).bgra[k];
        soit &iv1 = (*iptr1).bgra[k];
        soit &iv2 = (*iptr2).bgra[k];
        ov = iv1 * (1 - a) + iv2 * a;
      }
      optr++;
      iptr1++;
      iptr2++;
    }
    retourne res;
  }


  void blend(Image src)
  {
    si((src.impl->sx != sx) || (src.impl->sy != sy))
    {
      msg_erreur("ImageBmp::blend() : image source de dim différente : ({}x{}) != ({}x{})", sx,sy,src.impl->sx,src.impl->sy);
      retourne;
    }

    soit optr = (bgra *) bitmap;
    soit iptr = (const bgra *) src.impl->bitmap;

    pour(auto x = 0; x < sx; x++)
    {
      pour(auto y = 0; y < sy; y++)
      {
        float a = (*iptr).bgra[3] / 255.0f;
        pour(auto k = 0; k < 3; k++)
        {
          soit &ov = (*optr).bgra[k];
          soit &iv = (*iptr).bgra[k];
          ov = ov * (1 - a) + iv * a;
        }
        optr++;
        iptr++;
      }
    }
  }

  void resize(entier sx, entier sy)
  {
    this->sx = sx;
    this->sy = sy;
    si(bitmap != nullptr)
      free(bitmap);
    bitmap = nullptr;
    si(sx * sy)
    {
      bitmap = (bgra *) malloc(((size_t) sx)*((size_t) sy)*((size_t) 4));
      si(bitmap == nullptr)
        echec("Image resize({}x{} = {:.1f} Mo) : malloc error", sx, sy, 4.0f*sx*(sy*1e-6f));
    }
  }

  void def_couleur_dessin(const Couleur &c)
  {
    cd = {c.b, c.g, c.r, c.alpha};
  }
  void def_couleur_remplissage(const Couleur &c)
  {
    cf = {c.b, c.g, c.r, c.alpha};
  }
  void def_epaisseur(entier ep)
  {
    this->ep = ep;
  }


  void fleche(const Point &p0, const Point &p1, StyleLigne style, float lg = 5)
  {
    ligne(p0, p1, style);

    float α = atan2(p1.y-p0.y,p1.x-p0.x);

    float θ1 = 0.85*π  + α;
    float θ2 = -0.85*π + α;
    Point p;
    p.x = p1.x + lg * cos(θ1);
    p.y = p1.y + lg * sin(θ1);
    ligne(p1, p, style);

    p.x = p1.x + lg * cos(θ2);
    p.y = p1.y + lg * sin(θ2);
    ligne(p1, p, style);
  }

  void ligne(const Point &p0, const Point &p1, StyleLigne style)
  {
    si(p0.x == p1.x)
    {
      soit inc = 1;
      si(style == StyleLigne::POINTILLEE)
        inc *= 3;
      soit ymin = max(0, min(p0.y, p1.y)),
           ymax = min(sy - 1, max(p0.y, p1.y));

      pour(auto y = ymin; y <= ymax; y += inc)
        pour(auto k = -ep/2; k <= (ep-1)/2; k++)
          point({p0.x+k, y});
    }
    sinon si(p0.y == p1.y)
    {
      soit inc = 1;
      si(style == StyleLigne::POINTILLEE)
        inc *= 3;
      soit xmin = max(0, min(p0.x, p1.x)),
           xmax = min(sx - 1, max(p0.x, p1.x));
      pour(auto x = xmin; x <= xmax; x += inc)
        pour(auto k = -ep/2; k <= (ep-1)/2; k++)
          point({x, p0.y+k});
    }
    sinon
      pour(auto k = -ep/2; k <= (ep-1)/2; k++)
        ligne_aa_int(p0.x, p0.y+k, p1.x, p1.y+k, style);
  }

  inline bool est_dans_rdi(const Pointf &p)  const
  {
    retourne ((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy));
  }

  inline void point(const Point &p)
  {
    si((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
      pixel(p.x, p.y) =  cd;
  }
  inline void point(const Point &p, const Couleur &c)
  {
    si((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
      pixel(p.x, p.y) = {c.b, c.g, c.r, c.alpha};
  }

  bgra couleur_vers_bgra(const Couleur &c)
  {
    bgra res;
    res.bgra[0] = c.b;
    res.bgra[1] = c.g;
    res.bgra[2] = c.r;
    res.bgra[3] = c.alpha;
    retourne res;
  }

  inline bgra pixel_e(entier x, entier y) const
  {
    x = max(0,min(x, sx-1));
    y = max(0,min(y, sy-1));
    retourne *(bitmap + y * sx + x);
  }



  Image redim(entier l, entier h) const
  {
    si((l < sx/2) && (h < sy/2))
    {
      Image nv(sx/2, sy/2);

      // Decimation rapport = 1/2
      pour(entier y = 0; y + 1 < sy; y += 2)
      {
        pour(entier x = 0; x + 1 < sx; x += 2)
        {
          const soit &p00 = pixel(x, y),
                     &p01 = pixel(x+1,y),
                     &p10 = pixel(x,y+1),
                     &p11 = pixel(x+1,y+1);

          soit &po = nv.impl->pixel(x/2, y/2);
          pour(auto k = 0; k < 4; k++)
            po.bgra[k] = (p00.bgra[k]+p01.bgra[k]+p10.bgra[k]+p11.bgra[k])/4;
        }
      }

      retourne nv.redim(l, h);
    }


    // Interpolation bilinéaire
    Image res(l, h);
    pour(auto y = 0; y < h; y++)
    {
      pour(auto x = 0; x < l; x++)
      {
        soit xf = (x * 1.0f) / (l-1),
             yf = (y * 1.0f) / (h-1);

        soit xe = xf * sx, ye = yf * sy;

        soit xei = (entier) floor(xe), yei = (entier) floor(ye);
        soit rx = xe - xei, ry = ye - yei;

        soit p00 = pixel_e(xei, yei),
             p01 = pixel_e(xei+1, yei),
             p10 = pixel_e(xei, yei+1),
             p11 = pixel_e(xei+1, yei+1);

        pour(auto k = 0; k < 4; k++)
        {
          soit i1 = p00.bgra[k] * (1-rx) + p01.bgra[k] * rx,
               i2 = p10.bgra[k] * (1-rx) + p11.bgra[k] * rx;
          res.impl->pixel(x, y).bgra[k] = i1 * (1-ry) + i2 * ry;
        }
      }
    }

    retourne res;
  }

  void point(const Point &p, const Couleur &c, float α)
  {
    si((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
    {
      soit &ov = pixel(p.x, p.y);
      soit nv = couleur_vers_bgra(c);

      pour(auto k = 0; k < 4; k++)
        ov.bgra[k] = α * nv.bgra[k] + (1-α) * ov.bgra[k];
    }
  }

  void point(const Point &p, const bgra &c, float α)
  {
    si((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
    {
      soit &ov = pixel(p.x, p.y);
      pour(auto k = 0; k < 4; k++)
        ov.bgra[k] = α * c.bgra[k] + (1-α) * ov.bgra[k];
    }
  }

  void rectangle_plein_alpha(const Point &p0, const Point &p1, const bgra &c)
  {
    float α = c.bgra[3] / 255.0f;
    //soit nv = couleur_vers_bgra(c);
    soit xmin = max(0, min(p0.x, p1.x)),
         ymin = max(0, min(p0.y, p1.y)),
         xmax = min(sx-1,max(p0.x, p1.x)),
         ymax = min(sy-1,max(p0.y, p1.y));

    si(xmax >= xmin)
    {
      pour(auto y = ymin; y <= ymax; y++)
      {
        pour(auto x = xmin; x <= xmax; x++)
        {
          soit &ov = pixel(x, y);
          pour(auto k = 0; k < 4; k++)
            ov.bgra[k] = α * c.bgra[k] + (1-α) * ov.bgra[k];
        }
      }
    }
  }

  void rectangle(const Rect &r)
  {
    rectangle(r.tl(), r.br());
  }

  void rectangle(const Point &p0, const Point &p1)
  {
    ligne(p0, {p1.x, p0.y}, StyleLigne::PLEINE);
    ligne(p0, {p0.x, p1.y}, StyleLigne::PLEINE);
    ligne(p1, {p1.x, p0.y}, StyleLigne::PLEINE);
    ligne(p1, {p0.x, p1.y}, StyleLigne::PLEINE);
  }

  void cercle8(const Point &ctr, entier x_, entier y_, const bgra &c, float α)
  {
    // Dessine les 8 octants
    pour(auto k = 0; k < 4; k++)
    {
      entier x = x_, y = y_;
      si(k & 1)
        x = -x;
      si((k / 2) & 1)
        y = -y;

      pour(auto l = 0; l < 2; l++)
      {
        point({ctr.x + x, ctr.y + y}, c, α);
        swap(x, y);
      }
    }
  }

  void ellipse4(const Point &ctr, entier x_, entier y_, const bgra &c, float α)
  {
    // Dessine les 4 quadrants
    pour(auto k = 0; k < 4; k++)
    {
      entier x = x_, y = y_;
      si(k & 1)
        x = -x;
      si((k / 2) & 1)
        y = -y;
      point({ctr.x + x, ctr.y + y}, c, α);
    }
  }

  struct Courbe2
  {
    float x, y;
    // xy = f(t)
    virtual tuple<float,float> xy(float t) = 0;
    // xy = f(t)
    virtual tuple<float,float> der(float t) = 0;
  };

  struct Ellipse2: Courbe2
  {
    entier xc, yc, rx, ry;
    Ellipse2(entier xc, entier yc, entier rx, entier ry)
    {
      this->xc = xc;
      this->yc = yc;
      this->rx = rx;
      this->ry = ry;
      x = rx;
      y = 0;
    }
    tuple<float,float> xy(float t)
    {
      retourne {xc + rx * cos(2*π*t), yc + ry * sin(2*π*t)};
    }
    tuple<float,float> der(float t)
    {
      retourne {-2*π*rx * sin(2*π*t), 2*π*ry * cos(2*π*t)};
    }
  };

  void dessine_courbe_aa(Courbe2 *courbe, bgra cd)
  {
    soit nbitr = 0;
    soit t = 0.0f;
    tantque(t < 1)
    {
      si(nbitr++ > 100000)
      {
        msg_avert("Dessine courbe AA : trop d'itérations.");
        retourne;
      }
      soit [x, y]   = courbe->xy(t);
      soit [dx, dy] = courbe->der(t);
      //msg("t = {}, pos={}x{}, der={}x{}", t, x, y, dx, dy);

      // OK
      si(abs(dx) > abs(dy))
      {
        // Utilise deux pixels verticaux
        // cherche x entier
        entier xi   = floor(x);
        float dt = (xi - x) / dx;
        float yf = y + dt * dy;
        entier yi   = floor(yf);
        float α = yf - yi;

        point({xi, yi}, cd, 1-α);
        point({xi, yi+1}, cd, α);
        t += abs(0.5f / dx);
      }
      // NOK
      sinon
      {
        // Utilise deux pixels horizontaux
        // cherche x entier
        entier yi   = floor(y);
        float dt = (yi - y) / dy; // dt pour aller de y à floor(y)
        float xf = x + dt * dx;   // xf = point correspondant à floor(y)
        entier xi   = floor(xf);     // xi = point précédent
        float α = xf - xi;

        point({xi, yi}, cd, 1-α);
        point({xi+1, yi}, cd, α);
        t += abs(0.5f / dy);
      }
    }
  }

  void ellipse(const Point &p_0, const Point &p_1, bgra cd, bouléen avec_alpha)
  {
    Point p0{(p_0.x+p_1.x)/2, (p_0.y+p_1.y)/2};
    soit rx = (p_1.x - p_0.x) / 2.0f,
         ry = (p_1.y - p_0.y) / 2.0f;

    Ellipse2 el2(p0.x, p0.y, rx, ry);

    dessine_courbe_aa(&el2, cd);
  }

  void ellipse(const Point &p0, const Point &p1)
  {
    ellipse(p0, p1, cd, oui);
  }

  // D'après An Efficient Antialiasing Technique, Xiaolin Wu, 1991
  void cercle(const Point &p0, entier r, bgra cd, bouléen avec_alpha)
  {
    soit i = r, T = 0;

    cercle8(p0, i, 0, cd, 1.0f);

    pour(auto j = 1; j <= i; j++)
    {
      float rac = sqrt(r*r - j*j);
      entier d1 = ceil(rac);
      si(d1 < T)
        i--;

      float α = d1 - rac;

      cercle8(p0, i,   j, cd, avec_alpha ? 1.0f-α : 1.0f);
      cercle8(p0, i-1, j, cd, avec_alpha ? α : 1.0f);
      T = d1;
    }
  }

  void cercle(const Point &p0, entier r)
  {
    cercle(p0, r, cd, oui);
  }



  void rectangle_plein(const Point &p0, const Point &p1)
  {
    if(cf.bgra[3] != 255)
    {
      rectangle_plein_alpha(p0, p1, cf);
      return;
    }


    soit xmin = max(0, min(p0.x, p1.x)),
         ymin = max(0, min(p0.y, p1.y)),
         xmax = min(sx-1,max(p0.x, p1.x)),
         ymax = min(sy-1,max(p0.y, p1.y));

    si(xmax >= xmin)
      pour(auto y = ymin; y <= ymax; y++)
        fill(&(bitmap[y*sx+xmin]), &(bitmap[y*sx+xmax+1]), cf);
  }
  void ellipse_pleine(const Point &p0, const Point &p1)
  {
    msg_avert("TODO: image / ellipse pleine().");
  }
  void cercle_plein(const Point &p0, entier r)
  {
    pour(auto r2 = r; r2 >= 0; r2--)
      cercle(p0, r2, cf, r2 == r);
  }
  Dim texte_dim(const string &s, float dim)
  {
    retourne ftfonte->rendre(s, dim).get_dim();
  }

  entier calc_align(entier p0, entier dim, Alignement align)
  {
    si(align == Alignement::DEBUT)
      retourne p0;
    sinon si(align == Alignement::FIN)
      retourne p0 - dim;
    // Milieu
    retourne p0 - dim / 2;
  }

  void puts(const Point &p0, const string &s, float dim, Alignement align_horizontal, Alignement align_vertical)
  {
    si(s.empty())
      retourne;

    soit i = ftfonte->rendre(s, dim);
    si(i.empty())
    {
      msg_avert("Image::puts() : ftfonte->rendre : échec, s = [{}]", s);
      retourne;
    }

    Point p;
    p.x = calc_align(p0.x, i.sx(), align_horizontal);
    p.y = calc_align(p0.y, i.sy(), align_vertical);

    put_gamma(p, i);
  }
};

sptr<Font> Image::Impl::ftfonte;

Image::Image(entier sx, entier sy, void *data)
{
  impl = make_shared<Image::Impl>(sx, sy, data);
}

void Image::remplir(const Couleur &c)
{
  si(!empty())
  {
    def_couleur_remplissage(c);
    fill(impl->bitmap, impl->bitmap + sx()*sy(), impl->cf);
  }
}

void Image::resize(const Dim &d)
{
  resize(d.l, d.h);
}

const void *Image::data() const{retourne impl->data();}
void *Image::data(){retourne impl->data();}
Image Image::clone() const{retourne impl->clone();}
void Image::resize(entier sx, entier sy){impl->resize(sx,sy);}
Image Image::redim(entier l, entier h) const {retourne impl->redim(l,h);}
void Image::def_couleur_dessin(const Couleur &c){impl->def_couleur_dessin(c);}
void Image::def_couleur_remplissage(const Couleur &c){impl->def_couleur_remplissage(c);}
void Image::def_epaisseur(entier ep){impl->def_epaisseur(ep);}
void Image::ligne_aa(float x0, float y0, float x1, float y1, StyleLigne style){impl->ligne_aa(x0, y0, x1, y1, style);}
void Image::ligne(const Point &p0, const Point &p1, StyleLigne style){impl->ligne(p0, p1, style);}
void Image::fleche(const Point &p0, const Point &p1, StyleLigne style, float lg){impl->fleche(p0, p1, style, lg);}
//void Image::ligne(const Pointf &p0, const Pointf &p1, StyleLigne style){impl->ligne(p0, p1, style);}
void Image::point(const Point &p0){impl->point(p0);}
void Image::point(const Point &p0, const Couleur &c, float α){impl->point(p0, c, α);}
void Image::point(const Point &p0, const Couleur &c){impl->point(p0, c);}
void Image::rectangle(const Point &p0, const Point &p1){impl->rectangle(p0, p1);}
void Image::rectangle(const Rect &r){impl->rectangle(r);}
void Image::ellipse(const Point &p0, const Point &p1){impl->ellipse(p0, p1);}
void Image::cercle(const Point &p0, entier r){impl->cercle(p0, r);}
void Image::rectangle_plein(const Point &p0, const Point &p1){impl->rectangle_plein(p0, p1);}
void Image::ellipse_pleine(const Point &p0, const Point &p1){impl->ellipse_pleine(p0, p1);}
void Image::cercle_plein(const Point &p0, entier r){impl->cercle_plein(p0, r);}
Dim Image::texte_dim(const string &s, float dim){retourne impl->texte_dim(s, dim);}
void Image::puts(const Point &p0, const string &s, float dim, Alignement align_horizontal, Alignement align_vertical){impl->puts(p0, s, dim, align_horizontal, align_vertical);}
void Image::puti(const Point &pos, Image src, Rect rdi_source){impl->puti(pos, src, rdi_source);}
void Image::puti(const Rect &rdi, Image src){impl->puti(rdi, src);}
void Image::puti(const Point &pos, Image src){impl->puti(pos, src);}
void Image::puti_avec_gamma(const Point &pos, Image src){impl->puti_avec_gamma(pos, src);}
void Image::puti(const Point &pos, Image src, float γ){impl->puti(pos, src, γ);}
Image Image::rotation_90() const{retourne impl->rotation_90();}
void Image::blend(Image src){impl->blend(src);}
Image Image::blend_nv(Image src){retourne impl->blend_nv(src);}
void Image::enregister(const string &chemin) const {impl->enregister(chemin);}
void Image::charger(const string &chemin){impl->charger(chemin);}
Image Image::sous_image(entier x0, entier y0, entier l, entier h){retourne impl->sous_image(x0, y0, l, h);}
Dim Image::get_dim() const {retourne impl->dim();}
void Image::put_gamma(const Point &pos, Image src){impl->put_gamma(pos, src);}
}






