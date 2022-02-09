#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include <cassert>
#include <utility>

#if LIBTSD_USE_PNG
# include <png.h>
#endif


#define DBG_AXES(AA)


namespace tsd::vue {

  struct bgra
  {
    unsigned char bgra[4];
  };


struct Image::Impl
{
  int sx = 0, sy = 0;
  bgra *bitmap = nullptr;
  bgra cd = {0}, cf = {0};
  int ep = 1;
  static sptr<Font> ftfonte;


  inline float fpart(float x)
  {
      if(x > 0)
        return x - (int) x;
      else
        return x - ((int) (x)+1);

  }

  Dim dim() const
  {
    return {sx, sy};
  }

  void charger(const std::string &chemin)
  {

#   if LIBTSD_USE_PNG
    //msg("lecture png [{}]...", chemin);
    int width, height;
    png_byte color_type;
    png_byte bit_depth;
    png_bytep *row_pointers = NULL;

    FILE *fp = fopen(chemin.c_str(), "rb");

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png)
      echec("png_create_read_struct");

    png_infop info = png_create_info_struct(png);
    if(!info)
      echec("png_create_info_struct");

    if(setjmp(png_jmpbuf(png)))
      echec("png_jmpbuf(png)");

    png_init_io(png, fp);

    png_read_info(png, info);

    width      = png_get_image_width(png, info);
    height     = png_get_image_height(png, info);
    color_type = png_get_color_type(png, info);
    bit_depth  = png_get_bit_depth(png, info);

    // Read any color_type into 8bit depth, RGBA format.
    // See http://www.libpng.org/pub/png/libpng-manual.txt

    if(bit_depth == 16)
      png_set_strip_16(png);

    if(color_type == PNG_COLOR_TYPE_PALETTE)
      png_set_palette_to_rgb(png);

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
      png_set_expand_gray_1_2_4_to_8(png);

    if(png_get_valid(png, info, PNG_INFO_tRNS))
      png_set_tRNS_to_alpha(png);

    // These color_type don't have an alpha channel then fill it with 0xff.
    if(color_type == PNG_COLOR_TYPE_RGB ||
       color_type == PNG_COLOR_TYPE_GRAY ||
       color_type == PNG_COLOR_TYPE_PALETTE)
      png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    if(color_type == PNG_COLOR_TYPE_GRAY ||
       color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
      png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    //auto nrb = png_get_rowbytes(png,info);

    //msg("mallocs (h = {}, w = {}, nrb = {})...", height, width, nrb);

    row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for(int y = 0; y < height; y++) {
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
    for(auto y = 0; y < height; y++)
    {
      auto iptr = (unsigned char *) row_pointers[y];
      auto optr = ((unsigned char *) ((char *) this->bitmap)) + y * sx * 4;
      for(auto x = 0; x < width; x++)
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

  void enregister(const std::string &chemin)
  {

#   if LIBTSD_USE_PNG

    std::string ch = chemin;

    if((chemin.size() < 4) || (chemin[chemin.size()-4] != '.'))
      ch += ".png";

    //if(ch.ends_with(".png"))
    if((ch.size() >= 4) && (ch.substr(ch.size() - 4, 4) == ".png"))
    {
      int width = sx, height = sy;
      png_byte color_type = PNG_COLOR_TYPE_RGBA;
      png_byte bit_depth = 8;

      png_structp png_ptr;
      png_infop info_ptr;
      //int number_of_passes;
      png_bytep * row_pointers = (png_bytep *) malloc(sy * sizeof(png_bytep *));
      auto ptr = (unsigned char *) data();

      auto ptri = (unsigned char *) malloc(sx * sy * 4);
      auto tmp = ptri;
      for(auto i = 0; i < sx * sy; i++)
      {
        *tmp++ = ptr[2];
        *tmp++ = ptr[1];
        *tmp++ = ptr[0];
        *tmp++ = 255;

        ptr += 4;
      }

      for(auto y = 0; y < sy; y++)
        row_pointers[y] = ptri + y * sx * 4;


      /* create file */
      FILE *fp = fopen(ch.c_str(), "wb");
      if (!fp)
      {
        msg_erreur("[write_png_file] File {} could not be opened for writing", ch.c_str());
        return;
      }


      /* initialize stuff */
      png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

      if (!png_ptr)
      {
        msg_erreur("[write_png_file] png_create_write_struct failed");
        return;
      }

      info_ptr = png_create_info_struct(png_ptr);
      if (!info_ptr)
      {
        msg_erreur("[write_png_file] png_create_info_struct failed");
        return;
      }

      if (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during init_io");
        return;
      }

      png_init_io(png_ptr, fp);


      /* write header */
      if(setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during writing header");
        return;
      }

      png_set_IHDR(png_ptr, info_ptr, width, height,
                   bit_depth, color_type, PNG_INTERLACE_NONE,
                   PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

      png_write_info(png_ptr, info_ptr);


      /* write bytes */
      if (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during writing bytes");
        return;
      }

      png_write_image(png_ptr, row_pointers);


      /* end write */
      if (setjmp(png_jmpbuf(png_ptr)))
      {
        msg_erreur("[write_png_file] Error during end of write");
        return;
      }

      png_write_end(png_ptr, NULL);

      /* cleanup heap allocation */
      //for(y=0; y<height; y++)
        // free(row_pointers[y]);
      free(row_pointers);
      free(ptri);

      fclose(fp);
      return;
    }
#   endif

    msg_erreur("Enregistrement fichier image ({}) : pas de codec configuré pour cette extension.", chemin);
  }

  // On suppose que le signal d'entrée a été interpolé
  // pour des coordonnées entières de x
  void courbe_xy(int xmin, const ArrayXf &Y)
  {
    int xmax = Y.rows() + xmin - 1;
    for(auto x = xmin; x < xmax; x++)
    {
      // Plutôt vertical, il faut utiliser les pixels adjacents horizontaux
      if(std::abs(Y(x+1)-Y(x)) > 1)
      {

      }
      else
      {

      }



      /*for(auto y = ymin; y <= ymax; y++)
      {

      }*/
    }
  }


  void pta(int x, int y, float a)
  {
    if((x < 0) || (y < 0) || (x >= sx) || (y >= sy))
      return;

    auto ptr = bitmap + x + y * sx;
    auto ip = (bgra *) &cd;
    float a0 = std::sqrt(a);
    float b0 = 1 - a0;

    for(auto k = 0u; k < 4; k++)
      ptr->bgra[k] = ptr->bgra[k] * b0 + ip->bgra[k] * a0;
  }

  void ligne_aa(float x0, float y0, float x1, float y1, StyleLigne style)
  {
    //if( ((x0 == x1) || (y0 == y1)) (abs(x0 - floor(x0) < 1e-6)
    //    && (abs(x1 - floor(x1) < 1e-6)

    if(((x0 == x1) || (y0 == y1))
        && (abs(x0 - floor(x0) < 1e-6))
        && (abs(y0 - floor(y0) < 1e-6))
        && (abs(y1 - floor(y1) < 1e-6))
        && (abs(x1 - floor(x1) < 1e-6)))
      ligne({(int) x0, (int) y0}, {(int) x1, (int) y1}, style);
    else
    {
      for(auto k = -ep/2; k <= (ep-1)/2; k++)
        ligne_aa_int(x0, y0 + k, x1, y1 + k, style);
    }
  }

  void ligna_aa_nv(float x0, float y0, float x1, float y1)
  {

  }

  void ligne_aa_int(float x0, float y0, float x1, float y1, StyleLigne style = PLEINE)
  {
    //infos("Ligne aa : (%dx%d) -> (%dx%d)", x0, y0, x1, y1);

    {
      // TODO : trouver les points d'intersection avec la zone d'intérêt,
      // et ne tracer que dans cet intervalle

    }

    bool plutot_verticale = std::abs(y1 - y0) > std::abs(x1 - x0);

    if(plutot_verticale)
    {
      std::swap(x0 , y0);
      std::swap(x1 , y1);
    }
    if(x0 > x1)
    {
      std::swap(x0,x1);
      std::swap(y0,y1);
    }

    //compute the slope
    float dx = x1 - x0;
    float dy = y1 - y0;
    float pas = dx == 0 ? 1 : dy / dx;
    float y = y0;

    //Couleur c(cd);

    int xpxl1, xpxl2;

    {
      // handle first endpoint
      int xend = (int) round(x0);
      float yend = y0 + pas * (xend - x0);
      float xgap = 1 - fpart(x0 + 0.5);
      xpxl1 = xend; // this will be used in the main loop
      int ypxl1 = (int) (yend);
      if(plutot_verticale)
      {
        pta(ypxl1,   xpxl1, (1 - fpart(yend)) * xgap);
        pta(ypxl1+1, xpxl1, fpart(yend) * xgap);
      }
      else
      {
        pta(xpxl1, ypxl1  , (1 - fpart(yend)) * xgap);
        pta(xpxl1, ypxl1+1,  fpart(yend) * xgap);
      }
      y = yend + pas;
    }



    //intery := yend + gradient // first y-intersection for the main loop

    {
      // handle second endpoint
      int xend = round(x1);
      float yend = y1 + pas * (xend - x1);
      float xgap = fpart(x1 + 0.5);
      xpxl2 = xend; //this will be used in the main loop
      int ypxl2 = (int) (yend);
      if(plutot_verticale)
      {
        pta(ypxl2  , xpxl2, (1-fpart(yend)) * xgap);
        pta(ypxl2+1, xpxl2,  fpart(yend) * xgap);
      }
      else
      {
        pta(xpxl2, ypxl2,  (1-fpart(yend)) * xgap);
        pta(xpxl2, ypxl2+1, fpart(yend) * xgap);
      }
    }


    int inc = 1;
    if(style == POINTILLEE)
    {
      inc = 3;
    }

    // main loop
    if(plutot_verticale)
    {
      for(auto x = xpxl1+1; x <= xpxl2-1; x += inc)
      {
        float r = fpart(y);
        pta((int) y, x, 1-r);
        pta((int) y+1, x, r);
        y += pas * inc;
      }
    }
    else
    {
      for(auto x = xpxl1+1; x <= xpxl2-1; x += inc)
      {
        float r = fpart(y);
        pta(x, (int) y, 1-r);
        pta(x, (int) y+1, r);
        y += pas * inc;
      }
    }

  }

  ~Impl()
  {
    if(bitmap != nullptr)
      free(bitmap);
    bitmap = nullptr;
    sx = sy = 0;
  }

  Impl(int sx, int sy, void *data = nullptr)
  {
    this->sx = sx;
    this->sy = sy;
    if(sx * sy > 0)
    {
      bitmap = (bgra *) malloc(4*((size_t)sx)*((size_t)sy));

      if(bitmap == nullptr)
        echec("Echec allocation image : {} x {} ({:.1f} Mo)", sx, sy, 4.0f*sx*sy*1e-6);

      if(data != nullptr)
        memcpy(bitmap, data, sx * sy * 4);
      //else
        //memset(bitmap, 255, sx * sy * 4);
    }

    if(!ftfonte)
    {
      ftfonte = fonte_ft_creation();
    }
  }

  Image clone() const
  {
    Image res(sx, sy);
    memcpy(res.impl->bitmap, bitmap, sx * sy * 4);
    return res;
  }

  const void *data() const
  {
    return bitmap;
  }
  void *data()
  {
    return bitmap;
  }

  Image rotation_90() const
  {
    Image res(sy, sx);

    auto iptr = (const int32_t *) bitmap;
    auto optr = (int32_t *) res.impl->bitmap;

    optr += sx * sy - 1;
    for(auto x = 0; x < sx; x++)
      for(auto y = 0; y < sy; y++)
        *optr-- = iptr[(sx - 1 - x) + y * sx];

    return res;
  }




  inline bgra &pixel(int x, int y)
  {
    return *(bitmap + y * sx + x);
  }

  inline const bgra &pixel(int x, int y) const
  {
    return *(bitmap + y * sx + x);
  }


  // Pose uniquement le gamma, avec la couleur en cours
  void put_gamma(const Point &pos, Image src)
  {
    int isx = src.impl->sx, isy = src.impl->sy;

    if((isx + pos.x > sx) || (isy + pos.y > sy) || (pos.x < 0) || (pos.y < 0))
    {
      //msg_avert("ImageBmp:: put gamma : dépassement.");
      return;
    }

    for(auto y = 0; y < isy; y++)
    {
      for(auto x = 0; x < isx; x++)
      {
        auto iv = src.impl->pixel(x, y);
        float a = iv.bgra[3] / 255.0f;
        auto &ov = pixel(x + pos.x, y + pos.y);

        for(auto k = 0; k < 3; k++)
          ov.bgra[k]      = (1 - a) * ov.bgra[k] + a * cd.bgra[k];
        //ov.bgra[3]  = 255;//iv.bgra[3];
        //ov.bgra[3]  = 255 - (255 - iv.bgra[3]) / 4;
        ov.bgra[3]  = iv.bgra[3];
        if(iv.bgra[3] > 0)
          ov.bgra[3]  = 255 - (255 - iv.bgra[3]) / 2;
      }
    }
  }

  void puti_avec_gamma(const Point &pos, Image src)
  {

    /*if((src->sx != sx) || (src->sy != sy))
    {
      erreur("ImageBmp::blend() : image source de dim différente : (%dx%d) != (%dx%d)", sx,sy,src->sx,src->sy);
      return;
    }*/

    int isx = src.impl->sx, isy = src.impl->sy;


    //if((pos.x < 0) || (pos.y < 0))
    if(pos.x < 0)
    {
      Image src2 = src.sous_image(-pos.x, 0, src.sx() + pos.x, src.sy());
      puti_avec_gamma({0, pos.y}, src2);
      return;
    }
    if(pos.y < 0)
    {
      Image src2 = src.sous_image(0, -pos.y, src.sx(), src.sy() + pos.y);
      puti_avec_gamma({pos.x, 0}, src2);
      return;
    }

    /*if(isx + pos.x > sx)
    {
      int d =
    }*/


    if((isx + pos.x > sx) || (isy + pos.y > sy) || (pos.x < 0) || (pos.y < 0))
    {
      msg_avert("ImageBmp:: put image avec gamma : dépassement : dim cible = {}x{}, pos={}, dim src={},{}.",
          sx, sy, pos, isx, isy);
      return;
    }

    if((pos.x >= 2e9) || (pos.y >= 2e9))
    {
      msg_erreur("pos = {}", pos);
    }

    for(auto y = 0; y < isy; y++)
    {
      for(auto x = 0; x < isx; x++)
      {
        // Enlever cette assertion !

        if(!(((x + pos.x >= 0) && (y + pos.y >= 0) && (x + pos.x < sx) && (y + pos.y < sy))))
        {
          echec("x={},y={},isx={},isy={},sx={},sy={},p.x={},p.y={},isy + pos.y={}",x,y,isx,isy,sx,sy,pos.x,pos.y,isy + pos.y);
        }

        assert((x + pos.x >= 0) && (y + pos.y >= 0) && (x + pos.x < sx) && (y + pos.y < sy));


        auto iv = src.impl->pixel(x, y);
        float a = iv.bgra[3] / 255.0f;
        auto &ov = pixel(x + pos.x, y + pos.y);



        for(auto k = 0; k < 3; k++)
          ov.bgra[k] = (1 - a) * ov.bgra[k] + a * iv.bgra[k];
        ov.bgra[3]  = iv.bgra[3];
      }
    }
  }

  inline bgra *ptr(int x, int y)
  {
    return bitmap + x + y * sx;
  }

  Image sous_image(int x0, int y0, int l, int h)
  {
    if((x0 < 0) || (y0 < 0) || (x0 + l > sx) || (y0 + h > sy) || (l < 0) || (h < 0))
    {
      msg_avert("sous_image : hors borne.");
      return Image();
    }


    Image res(l, h);
    for(auto y = 0; y < h; y++)
      memcpy(res.impl->bitmap + y * l, this->bitmap + x0 + (y0 + y) * sx, l * 4);
    return res;
  }

  void puti(const Point &pos, Image src, Rect rdi_source)
  {
    for(auto y = 0; y < rdi_source.h; y++)
    {
      memcpy(ptr(pos.x, pos.y + y), src.impl->ptr(rdi_source.x, rdi_source.y + y), rdi_source.l * 4);
    }
  }

  void puti(const Point &pos, Image src)
  {
    if((pos.x < 0) || (pos.y < 0) || (pos.x + src.impl->sx > sx) || (pos.y + src.impl->sy > sy))
    {
      DBG_AXES(msg("puti : pos = {}", pos);)
      Point pos2 = pos;

      int x0 = std::max(-pos.x, 0);
      int y0 = std::max(-pos.y, 0);

      if(pos2.x < 0)
        pos2.x = 0;
      if(pos2.y < 0)
        pos2.y = 0;

      int l = src.impl->sx - x0;
      int h = src.impl->sy - y0;

      if(l + pos2.x > sx)
        l = sx - pos2.x;
      if(h + pos2.y > sy)
        h = sy - pos2.y;

      auto extrait = src.sous_image(x0, y0, l, h);
      puti(pos2, extrait);
      return;
    }
    //msg("puti : src dim = {}x{}, dst dim = {}x{}, pos={}x{} (place = {}x{})", src->sx, src->sy, sx, sy, pos.x, pos.y, sx-pos.x, sy-pos.y);

    auto dimx = src.impl->sx;
    if(dimx + pos.x > sx)
      dimx = sx - pos.x;

    for(auto y = pos.y; (y < pos.y + src.impl->sy) && (y < sy); y++)
    {
      if(dimx > 0)
        memcpy(ptr(pos.x, y), src.impl->ptr(0, y - pos.y), dimx * 4);
    }
  }

  void puti(const Rect &rdi, Image src)
  {
    if((rdi.l == src.impl->sx) && (rdi.h == src.impl->sy))
    {
      puti(Point{rdi.x, rdi.y}, src);
      return;
    }

    if((rdi.l <= 0) || (rdi.h <= 0))
    {
      msg_avert("RDI = {} x {} !!!", rdi.l, rdi.h);
      return;
    }

    // Image vide
    if(src.sx() * src.sy() == 0)
      return;


    //auto srci = std::dynamic_pointer_cast<ImageBmp>(src);
    //tsd_assert(srci);


    msg_avert("Attention: redim {} -> {}x{}", src.get_dim(), rdi.l, rdi.h);
    puti(rdi, src.redim(rdi.l, rdi.h));

    //msg_erreur("TODO: puti (rdi = {}x{}, source = {}x{}).", rdi.l, rdi.h, src->sx, src->sy);
  }

  void puti(const Point &pos, Image src, float γ)
  {

  }



  Image blend_nv(Image src)
  {
    auto s1 = this;
    auto s2 = src.impl;
    tsd_assert(s2);
    if((s1->sx != s2->sx) || (s1->sy != s2->sy))
    {
      msg_erreur("ImageBmp::blend() : image source de dim différente : ({}) != ({})", s1->dim(), src.get_dim());
      return Image();
    }

    Image res(s1->sx, s2->sy);
    //auto r2 = std::dynamic_pointer_cast<ImageBmp>(res);

    auto optr = (bgra *) res.impl->bitmap;
    auto iptr1 = (const bgra *) s1->bitmap;
    auto iptr2 = (const bgra *) s2->bitmap;

    for(auto i = 0; i < sy*sx; i++)
    {
      float a = (*iptr2).bgra[3] / 255.0f;
      for(auto k = 0; k < 3; k++)
      {
        auto &ov  = (*optr).bgra[k];
        auto &iv1 = (*iptr1).bgra[k];
        auto &iv2 = (*iptr2).bgra[k];
        ov = iv1 * (1 - a) + iv2 * a;
      }
      optr++;
      iptr1++;
      iptr2++;
    }
    return res;
  }


  void blend(Image src)
  {
    if((src.impl->sx != sx) || (src.impl->sy != sy))
    {
      msg_erreur("ImageBmp::blend() : image source de dim différente : ({}x{}) != ({}x{})", sx,sy,src.impl->sx,src.impl->sy);
      return;
    }

    auto optr = (bgra *) bitmap;
    auto iptr = (const bgra *) src.impl->bitmap;

    for(auto x = 0; x < sx; x++)
    {
      for(auto y = 0; y < sy; y++)
      {
        float a = (*iptr).bgra[3] / 255.0f;
        for(auto k = 0; k < 3; k++)
        {
          auto &ov = (*optr).bgra[k];
          auto &iv = (*iptr).bgra[k];
          ov = ov * (1 - a) + iv * a;
        }
        optr++;
        iptr++;
      }
    }
  }

  void resize(int sx, int sy)
  {
    this->sx = sx;
    this->sy = sy;
    if(bitmap != nullptr)
      free(bitmap);
    bitmap = nullptr;
    if(sx * sy)
    {
      bitmap = (bgra *) malloc(((size_t) sx)*((size_t) sy)*((size_t) 4));
      if(bitmap == nullptr)
      {
        echec("Image resize({}x{} = {:.1f} Mo) : malloc error", sx, sy, 4.0f*sx*(sy*1e-6f));
      }
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
  void def_epaisseur(int ep)
  {
    this->ep = ep;
  }


  void fleche(const Point &p0, const Point &p1, StyleLigne style, float lg = 5)
  {
    ligne(p0, p1, style);

    float alpha = atan2(p1.y-p0.y,p1.x-p0.x);

    float θ1 = 0.85*π  + alpha;
    float θ2 = -0.85*π + alpha;
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
    if(p0.x == p1.x)
    {
      int inc = 1;
      if(style == StyleLigne::POINTILLEE)
        inc *= 3;
      auto ymin = std::max(0, std::min(p0.y, p1.y));
      auto ymax = std::min(sy - 1, std::max(p0.y, p1.y));

      for(auto y = ymin; y <= ymax; y += inc)
        for(auto k = -ep/2; k <= (ep-1)/2; k++)
          point({p0.x+k, y});
    }
    else if(p0.y == p1.y)
    {
      int inc = 1;
      if(style == StyleLigne::POINTILLEE)
        inc *= 3;
      auto xmin = std::max(0, std::min(p0.x, p1.x));
      auto xmax = std::min(sx - 1, std::max(p0.x, p1.x));
      for(auto x = xmin; x <= xmax; x += inc)
        for(auto k = -ep/2; k <= (ep-1)/2; k++)
          point({x, p0.y+k});
    }
    else
    {
      for(auto k = -ep/2; k <= (ep-1)/2; k++)
        ligne_aa_int(p0.x, p0.y+k, p1.x, p1.y+k, style);
    }
  }
  inline void point(const Point &p)
  {
    if((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
      pixel(p.x, p.y) =  cd;
  }
  inline void point(const Point &p, const Couleur &c)
  {
    if((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
      pixel(p.x, p.y) = {c.b, c.g, c.r, c.alpha};
  }

  bgra couleur_vers_bgra(const Couleur &c)
  {
    bgra res;
    res.bgra[0] = c.b;
    res.bgra[1] = c.g;
    res.bgra[2] = c.r;
    res.bgra[3] = c.alpha;
    return res;
  }

  inline bgra pixel_e(int x, int y) const
  {
    x = std::min(x, sx-1);
    y = std::min(y, sy-1);
    return *(bitmap + y * sx + x);
  }



  Image redim(int l, int h) const
  {
    if((l < sx/2) && (h < sy/2))
    {
      Image nv(sx/2, sy/2);

      // Decimation rapport = 1/2
      for(int y = 0; y + 1 < sy; y += 2)
      {
        for(int x = 0; x + 1 < sx; x += 2)
        {
          const auto &p00 = pixel(x, y),
                     &p01 = pixel(x+1,y),
                     &p10 = pixel(x,y+1),
                     &p11 = pixel(x+1,y+1);

          auto &po = nv.impl->pixel(x/2, y/2);
          for(auto k = 0; k < 4; k++)
            po.bgra[k] = (p00.bgra[k]+p01.bgra[k]+p10.bgra[k]+p11.bgra[k])/4;
        }
      }

      return nv.redim(l, h);
    }


    // Interpolation bilinéaire
    Image res(l, h);
    for(auto y = 0; y < h; y++)
    {
      for(auto x = 0; x < l; x++)
      {
        float xf = (x * 1.0f) / (l-1);
        float yf = (y * 1.0f) / (h-1);

        float xe = xf * sx, ye = yf * sy;

        int xei = floor(xe), yei = floor(ye);
        float rx = xe - xei, ry = ye - yei;

        auto p00 = pixel_e(xei, yei);
        auto p01 = pixel_e(xei+1, yei);
        auto p10 = pixel_e(xei, yei+1);
        auto p11 = pixel_e(xei+1, yei+1);

        for(auto k = 0; k < 4; k++)
        {
          auto i1 = p00.bgra[k] * (1-rx) + p01.bgra[k] * rx;
          auto i2 = p10.bgra[k] * (1-rx) + p11.bgra[k] * rx;
          res.impl->pixel(x, y).bgra[k] = i1 * (1-ry) + i2 * ry;
        }
      }
    }

    return res;
  }

  void point(const Point &p, const Couleur &c, float alpha)
  {
    if((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
    {
      auto &ov = pixel(p.x, p.y);
      auto nv = couleur_vers_bgra(c);

      for(auto k = 0; k < 4; k++)
        ov.bgra[k] = alpha * nv.bgra[k] + (1-alpha) * ov.bgra[k];
    }
  }

  void point(const Point &p, const bgra &c, float alpha)
  {
    if((p.x >= 0) && (p.y >= 0) && (p.x < sx) && (p.y < sy))
    {
      auto &ov = pixel(p.x, p.y);
      for(auto k = 0; k < 4; k++)
        ov.bgra[k] = alpha * c.bgra[k] + (1-alpha) * ov.bgra[k];
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

  void cercle8(const Point &ctr, int x_, int y_, const bgra &c, float alpha)
  {
    // Dessinne les 8 octants
    for(auto k = 0; k < 4; k++)
    {
      int x = x_, y = y_;
      if(k & 1)
        x = -x;
      if((k / 2) & 1)
        y = -y;

      for(auto l = 0; l < 2; l++)
      {
        point({ctr.x + x, ctr.y + y}, c, alpha);
        std::swap(x, y);
      }
    }
  }

  void ellipse4(const Point &ctr, int x_, int y_, const bgra &c, float alpha)
  {
    // Dessinne les 4 quadrants
    for(auto k = 0; k < 4; k++)
    {
      int x = x_, y = y_;
      if(k & 1)
        x = -x;
      if((k / 2) & 1)
        y = -y;
      point({ctr.x + x, ctr.y + y}, c, alpha);
    }
  }

  struct Courbe2
  {
    float x, y;
    // xy = f(t)
    virtual std::tuple<float,float> xy(float t) = 0;
    // xy = f(t)
    virtual std::tuple<float,float> der(float t) = 0;
  };

  struct Ellipse2: Courbe2
  {
    int xc, yc, rx, ry;
    Ellipse2(int xc, int yc, int rx, int ry)
    {
      this->xc = xc;
      this->yc = yc;
      this->rx = rx;
      this->ry = ry;

      x = rx;
      y = 0;
    }
    std::tuple<float,float> xy(float t)
    {
      return {xc + rx * cos(2*π*t), yc + ry * sin(2*π*t)};
    }
    std::tuple<float,float> der(float t)
    {
      return {-2*π*rx * sin(2*π*t), 2*π*ry * cos(2*π*t)};
    }
  };

  void dessine_courbe_aa(Courbe2 *courbe, bgra cd)
  {
    int nbitr = 0;
    float t = 0;
    while(t < 1)
    {
      if(nbitr++ > 100000)
      {
        msg_avert("Dessine courbe AA : trop d'itérations.");
        return;
      }
      auto [x, y]   = courbe->xy(t);
      auto [dx, dy] = courbe->der(t);
      //msg("t = {}, pos={}x{}, der={}x{}", t, x, y, dx, dy);
      if(std::abs(dx) > std::abs(dy))
      {
        // Utilise deux pixels verticaux
        // cherche x entier
        int xi   = floor(x);
        float dt = (xi - x) / dx;
        float yf = y + dt * dy;
        int yi   = floor(yf);
        float alpha = yf - yi;

        point({xi, yi}, cd, 1-alpha);
        point({xi, yi+1}, cd, alpha);
        t += abs(0.5f / dx);
      }
      else
      {
        // Utilise deux pixels horizontaux
        // cherche x entier
        int yi   = floor(y);
        float dt = (yi - y) / dy; // dt pour aller de y à floor(y)
        float xf = x + dt * dx;   // xf = point correspondant à floor(y)
        int xi   = floor(xf);     // xi = point précédent
        float alpha = xf - xi;

        point({xi, yi}, cd, 1-alpha);
        point({xi+1, yi}, cd, alpha);
        t += abs(0.5f / dy);
      }
    }
  }

  void ellipse(const Point &p_0, const Point &p_1, bgra cd, bool avec_alpha)
  {
    Point p0{(p_0.x+p_1.x)/2, (p_0.y+p_1.y)/2};
    float rx = (p_1.x - p_0.x) / 2.0f;
    float ry = (p_1.y - p_0.y) / 2.0f;

    Ellipse2 el2(p0.x, p0.y, rx, ry);

    dessine_courbe_aa(&el2, cd);
  }

  void ellipse(const Point &p0, const Point &p1)
  {
    ellipse(p0, p1, cd, true);
  }

  // D'après An Efficient Antialiasing Technique, Xiaolin Wu, 1991
  void cercle(const Point &p0, int r, bgra cd, bool avec_alpha)
  {
    int i = r, T = 0;
    cercle8(p0, i, 0, cd, 1.0f);

    for(auto j = 1; j <= i; j++)
    {
      float rac = std::sqrt(r*r - j*j);
      int d1 = std::ceil(rac);
      if(d1 < T)
        i--;

      float alpha = d1 - rac;

      cercle8(p0, i,   j, cd, avec_alpha ? 1.0f-alpha : 1.0f);
      cercle8(p0, i-1, j, cd, avec_alpha ? alpha : 1.0f);
      T = d1;
    }
  }

  void cercle(const Point &p0, int r)
  {
    cercle(p0, r, cd, true);
  }

  void rectangle_plein(const Point &p0, const Point &p1)
  {
    int xmin = std::max(0, std::min(p0.x, p1.x));
    int ymin = std::max(0, std::min(p0.y, p1.y));
    int xmax = std::min(sx-1,std::max(p0.x, p1.x));
    int ymax = std::min(sy-1,std::max(p0.y, p1.y));

    if(xmax >= xmin)
      for(auto y = ymin; y <= ymax; y++)
        std::fill(&(bitmap[y*sx+xmin]), &(bitmap[y*sx+xmax+1]), cf);
  }
  void ellipse_pleine(const Point &p0, const Point &p1)
  {
    msg_avert("TODO: image / ellipse pleine().");
  }
  void cercle_plein(const Point &p0, int r)
  {
    for(auto r2 = r; r2 >= 0; r2--)
      cercle(p0, r2, cf, r2 == r);
  }
  Dim texte_dim(const std::string &s, float dim)
  {
    return ftfonte->rendre(s, dim).get_dim();
  }

  int calc_align(int p0, int dim, Alignement align)
  {
    if(align == Alignement::DEBUT)
      return p0;
    else if(align == Alignement::FIN)
      return p0 - dim;
    // Milieu
    return p0 - dim / 2;
  }

  void puts(const Point &p0, const std::string &s, float dim, Alignement align_horizontal, Alignement align_vertical)
  {
    if(s.empty())
      return;

    auto i = ftfonte->rendre(s, dim);
    if(i.empty())
    {
      msg_avert("Image::puts() : ftfonte->rendre : échec, s = [{}]", s);
      return;
    }

    Point p;
    p.x = calc_align(p0.x, i.sx(), align_horizontal);
    p.y = calc_align(p0.y, i.sy(), align_vertical);

    put_gamma(p, i);
  }
};

sptr<Font> Image::Impl::ftfonte;

Image::Image(int sx, int sy, void *data)
{
  impl = std::make_shared<Image::Impl>(sx, sy, data);
}

void Image::remplir(const Couleur &c)
{
  if(!empty())
  {
    def_couleur_remplissage(c);
    std::fill(impl->bitmap, impl->bitmap + sx()*sy(), impl->cf);
  }
}

void Image::resize(const Dim &d)
{
  resize(d.l, d.h);
}

const void *Image::data() const{return impl->data();}
void *Image::data(){return impl->data();}
Image Image::clone() const{return impl->clone();}
void Image::resize(int sx, int sy){impl->resize(sx,sy);}
Image Image::redim(int l, int h) const {return impl->redim(l,h);}
void Image::def_couleur_dessin(const Couleur &c){impl->def_couleur_dessin(c);}
void Image::def_couleur_remplissage(const Couleur &c){impl->def_couleur_remplissage(c);}
void Image::def_epaisseur(int ep){impl->def_epaisseur(ep);}
void Image::ligne_aa(float x0, float y0, float x1, float y1, StyleLigne style){impl->ligne_aa(x0, y0, x1, y1, style);}
void Image::ligne(const Point &p0, const Point &p1, StyleLigne style){impl->ligne(p0, p1, style);}
void Image::fleche(const Point &p0, const Point &p1, StyleLigne style, float lg){impl->fleche(p0, p1, style, lg);}
//void Image::ligne(const Pointf &p0, const Pointf &p1, StyleLigne style){impl->ligne(p0, p1, style);}
void Image::point(const Point &p0){impl->point(p0);}
void Image::point(const Point &p0, const Couleur &c, float alpha){impl->point(p0, c, alpha);}
void Image::point(const Point &p0, const Couleur &c){impl->point(p0, c);}
void Image::rectangle(const Point &p0, const Point &p1){impl->rectangle(p0, p1);}
void Image::rectangle(const Rect &r){impl->rectangle(r);}
void Image::ellipse(const Point &p0, const Point &p1){impl->ellipse(p0, p1);}
void Image::cercle(const Point &p0, int r){impl->cercle(p0, r);}
void Image::rectangle_plein(const Point &p0, const Point &p1){impl->rectangle_plein(p0, p1);}
void Image::ellipse_pleine(const Point &p0, const Point &p1){impl->ellipse_pleine(p0, p1);}
void Image::cercle_plein(const Point &p0, int r){impl->cercle_plein(p0, r);}
Dim Image::texte_dim(const std::string &s, float dim){return impl->texte_dim(s, dim);}
void Image::puts(const Point &p0, const std::string &s, float dim, Alignement align_horizontal, Alignement align_vertical){impl->puts(p0, s, dim, align_horizontal, align_vertical);}
void Image::puti(const Point &pos, Image src, Rect rdi_source){impl->puti(pos, src, rdi_source);}
void Image::puti(const Rect &rdi, Image src){impl->puti(rdi, src);}
void Image::puti(const Point &pos, Image src){impl->puti(pos, src);}
void Image::puti_avec_gamma(const Point &pos, Image src){impl->puti_avec_gamma(pos, src);}
void Image::puti(const Point &pos, Image src, float γ){impl->puti(pos, src, γ);}
Image Image::rotation_90() const{return impl->rotation_90();}
void Image::blend(Image src){impl->blend(src);}
Image Image::blend_nv(Image src){return impl->blend_nv(src);}
void Image::enregister(const std::string &chemin){impl->enregister(chemin);}
void Image::charger(const std::string &chemin){impl->charger(chemin);}
Image Image::sous_image(int x0, int y0, int l, int h){return impl->sous_image(x0, y0, l, h);}
Dim Image::get_dim() const {return impl->dim();}
void Image::put_gamma(const Point &pos, Image src){impl->put_gamma(pos, src);}
}






