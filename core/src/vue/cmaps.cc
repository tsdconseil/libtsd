#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/vue.hpp"
#include <map>




namespace tsd::vue
{
  //template<typename T>
    //using std::make_shared<T>;

  struct CMapIlin: CMap
  {
    std::vector<std::vector<float>> pts;

    void calc(float t, float &r, float &v, float &b)
    {
      entier i = 0;

      tsd_assert_msg(pts.size() > 1, "Carte de couleur linéaire : au moins deux éléments attendu.");

      si(t <= pts[0][0])
      {
        r = pts[0][1];
        v = pts[0][2];
        b = pts[0][3];
        retourne;
      }

      tantque(((i+1) < (entier) pts.size()) && (t >= pts[i+1][0]))
        i++;

      si((i+1) >= (entier) pts.size())
      {
        r = pts.back()[1];
        v = pts.back()[2];
        b = pts.back()[3];
        retourne;
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


  struct CMapPM: CMapIlin
  {
    CMapPM()
    {
      pts = {

          {0.00, 0.0,   0.0,  0.5},
          {0.25, 0.0,   0.0,  1.0},
          {0.50, 1.0,   1.0,  1.0},
          {0.75, 1.0,   0.0,  0.0},
          {1.00, 0.5,   0.0,  0.0}


          /*{0.0, 0.0,   0.0,  1.0},
          {0.5, 1.0,   1.0,  1.0},
          {1.0, 1.0,   0.0,  0.0}*/
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
      si(t < 0)
        b = r = v = 1;
      sinon si (t < 1)
        b = r = v = 1 - t;
      sinon
        b = r = v = 0;
    }
  };

  struct CMapMonoInv: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      si(t < 0)
        b = r = v = 0;
      sinon si (t < 1)
        b = r = v = t;
      sinon
        b = r = v = 1;
    }
  };

  struct CMapMonoPer: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      si(t < 0)
        b = r = v = 1;
      sinon si (t <= 0.5)
        b = r = v = 2 * t;
      sinon si(t <= 1.0)
        b = r = v = 1 - 2 * (t - 0.5);
      sinon
        b = r = v = 0;
    }
  };



  struct CMapJet: CMap
  {
    void calc(float t, float &r, float &v, float &b)
    {
      r = v = b = 0;

      si(t < 0)
      {
        b = 0.5;
      }
      sinon si(t < 1.0/8)
      {
        b = 0.5 + t * 4;
      }
      sinon si(t < 3.0/8)
      {
        b = 1.0;
        v = (t - 1.0/8) * 4;
      }
      sinon si(t < 5.0/8)
      {
        v = 1.0;
        r = (t - 3.0/8) * 4;
        b = 1.0 - (t - 3.0/8) * 4;
      }
      sinon si(t < 7.0/8)
      {
        r = 1;
        v = 1.0 - (t - 5.0/8) * 4;
      }
      sinon si(t < 1)
      {
        r = 1.0 - (t - 7.0/8) * 4;
      }
      sinon
        r = 0.5;
    }
  };


  Couleur CMap::couleur(float t)
  {
    float R = 0, V = 0, B = 0;
    calc(t, R, V, B);
    retourne Couleur{255*R, 255*V, 255*B};
  }

  template<typename T>
  sptr<T> ms()
  {
    retourne std::make_shared<T>();
  }

  static std::map<std::string, sptr<CMap>> cmaps
    =
    {
        {"mono",      ms<CMapMono>()},
        {"mono-inv",  ms<CMapMonoInv>()},
        {"mono-per",  ms<CMapMonoPer>()},
        {"hsv",       ms<CMapHSV>()},
        {"jet",       ms<CMapJet>()},
        {"parula",    ms<CMapParula>()},
        {"pm",        ms<CMapPM>()}
    };

  sptr<CMap> cmap_parse(const std::string &nom)
  {
    si(cmaps.count(nom) == 0)
    {
      msg_erreur("Cmap inconnue : {}", nom);
      retourne cmaps["jet"];
    }
    retourne cmaps[nom];
  }

  void cmap_affiche(const std::string &nom, sptr<CMap> cmap)
  {
    Figure f(nom);

    soit n = 100u;
    soit x = linspace(0, 1, n);
    Vecf r(n), v(n), b(n);

    pour(auto i = 0u; i < n; i++)
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

