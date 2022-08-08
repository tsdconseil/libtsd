#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/vue.hpp"


namespace tsd::vue
{
  struct CMapIlin: CMap
  {
    std::vector<std::vector<float>> pts;

    void calc(float t, float &r, float &v, float &b)
    {
      int i = 0;

      tsd_assert_msg(pts.size() > 1, "Carte de couleur linéaire : au moins deux éléments attendu.");

      if(t <= pts[0][0])
      {
        r = pts[0][1];
        v = pts[0][2];
        b = pts[0][3];
        return;
      }

      while(((i+1) < (int) pts.size()) && (t >= pts[i+1][0]))
        i++;

      if((i+1) >= (int) pts.size())
      {
        r = pts.back()[1];
        v = pts.back()[2];
        b = pts.back()[3];
        return;
      }

      // pts[i] <= t < pts[i+1]
      float a = (t - pts[i][0]) / (pts[i+1][0] - pts[i][0]);

      r = (1 - a) * pts[i][1] + a * pts[i+1][1];
      v = (1 - a) * pts[i][2] + a * pts[i+1][2];
      b = (1 - a) * pts[i][3] + a * pts[i+1][3];
    }
  };

  struct CMapParula: CMapIlin
  {
    CMapParula()
    {
      pts = {
          {0,     0.26710521,  0.03311059,  0.6188155},
          {0.125, 0.14990781,  0.28892902,  0.59136956},
          {0.25,  0.1668964 ,  0.41027528,  0.50959299},
          {0.375, 0.19159194,  0.516384  ,  0.47597744},
          {0.5,   0.16400433,  0.63146734,  0.39450263},
          {0.625, 0.3704552 ,  0.72530195,  0.18639333},
          {0.75,  0.72692898,  0.75698787,  0.1691954},
          {0.875, 0.99505988,  0.78542889,  0.32106514},
          {1.0,   0.98680727,  0.95697596,  0.12661626}
      };
    }
  };

//                              [0, 0, 0, 1, 1, 1, 0],
//                      'green':[0, 1, 1, 1, 0, 0, 0],
//                      'red':  [1, 1, 0, 0, 0, 1, 1]}
  struct CMapHSV: CMapIlin
  {
    CMapHSV()
    {
      pts = {
          {1.0f/6, 1, 0, 0},
          {2.0f/6, 1, 1, 0},
          {3.0f/6, 0, 1, 0},
          {4.0f/6, 0, 1, 1},
          {5.0f/6, 0, 0, 1},
          {1.0,    1, 0, 0}
      };
    }
  };

  struct CMapMono: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      if(t < 0)
        b = r = v = 1;
      else if (t < 1)
        b = r = v = 1 - t;
      else
        b = r = v = 0;
    }
  };

  struct CMapMonoInv: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      if(t < 0)
        b = r = v = 0;
      else if (t < 1)
        b = r = v = t;
      else
        b = r = v = 1;
    }
  };

  struct CMapMonoPer: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      if(t < 0)
        b = r = v = 1;
      else if (t <= 0.5)
        b = r = v = 2 * t;
      else if(t <= 1.0)
        b = r = v = 1 - 2 * (t - 0.5);
      else
        b = r = v = 0;
    }
  };



  struct CMapJet: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      r = v = b = 0;

      if(t < 0)
      {
        b = 0.5;
      }
      else if(t < 1.0/8)
      {
        b = 0.5 + t * 4;
      }
      else if(t < 3.0/8)
      {
        b = 1.0;
        v = (t - 1.0/8) * 4;
      }
      else if(t < 5.0/8)
      {
        v = 1.0;
        r = (t - 3.0/8) * 4;
        b = 1.0 - (t - 3.0/8) * 4;
      }
      else if(t < 7.0/8)
      {
        r = 1;
        v = 1.0 - (t - 5.0/8) * 4;
      }
      else if(t < 1)
      {
        r = 1.0 - (t - 7.0/8) * 4;
      }
      else
        r = 1.0f;
    }
  };


  Couleur CMap::couleur(float t)
  {
    float R = 0, V = 0, B = 0;
    calc(t, R, V, B);
    return Couleur{255*R, 255*V, 255*B};
  }

  sptr<CMap> cmap_parse(const std::string &nom)
  {
    if(nom == "mono")
      return std::make_shared<CMapMono>();
    else if(nom == "mono-inv")
      return std::make_shared<CMapMonoInv>();
    else if(nom == "mono-per")
      return std::make_shared<CMapMonoPer>();
    else if(nom == "hsv")
      return std::make_shared<CMapHSV>();
    else if(nom == "jet")
      return std::make_shared<CMapJet>();
    else if(nom == "parula")
      return std::make_shared<CMapParula>();
    msg_erreur("Cmap inconnue : {}", nom);
    return std::make_shared<CMapJet>();
  }

  void cmap_affiche(const std::string &nom, sptr<CMap> cmap)
  {
    Figure f(nom);

    auto n = 100u;
    ArrayXf x = linspace(0, 1, n), r(n), v(n), b(n);

    for(auto i = 0u; i < n; i++)
    {
      float R, V, B;
      cmap->calc(x(i), R, V, B);
      r(i) = R;
      v(i) = V;
      b(i) = B;
    }

    f.plot(x, r, "r-");
    f.plot(x, v, "g-");
    f.plot(x, b, "b-");

    f.afficher();
  }
}

