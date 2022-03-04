#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"

#define DBG(AA)


#include <deque>

using namespace std;


namespace tsd::vue
{

  // Définit les dimensions
  static const auto
    ECART_X_VALEURS   = 4,
    HAUTEUR_X_VALEURS = 20,
    HAUTEUR_X_LABEL   = 25,
    HAUTEUR_TITRE     = 25,
    HAUTEUR_Y_VALEUR  = 18;//13;

  static const auto
    FONTE_TITRE    = 1.0,
    FONTE_X_LABEL  = 0.6;

  void ConfigGrille::from_style(const string &style)
  {
    if(style == ".")
    {
      lg_tiret = 1;
      lg_trou  = 2;
    }
    else if(style == "-")
    {
      lg_tiret = 10;
      lg_trou  = 0;
    }
    else if(style == "--")
    {
      lg_tiret = 5;
      lg_trou  = 5;
    }
    else
    {
      lg_tiret = -1;
      lg_trou  = -1;
    }
  }

  string ConfigGrille::vers_style() const
  {
    if((lg_tiret == 1) && (lg_trou == 2))
      return ".";
    else if((lg_tiret == 10) && (lg_trou == 0))
      return "-";
    else if((lg_tiret == 5) && (lg_trou == 5))
      return "--";
    else
      return "";
  }

  bool ConfigAxes::est_mat_sur_clair() const
  {
    int score_grid  = grille_majeure.couleur.lumi();
    int score_back = couleur_arriere_plan.lumi(); // 255

    if(grille_majeure.couleur.vers_rgba() == 0) // Par défaut
      score_grid = 255 - score_back;

    return score_back > score_grid;
  }

  ConfigAxes::Legende::Position::Position(const string &id)
  {
    struct PosCart {string code; Type pos;};
    static const PosCart codes[] =
    {
        {"ne", NORD_EST},
        {"no", NORD_OUEST},
        {"so", SUD_OUEST},
        {"se", SUD_EST},
        {"ur", NORD_EST},
        {"ul", NORD_OUEST},
        {"ll", SUD_OUEST},
        {"lr", SUD_EST}
    };
    type = NORD_EST;
    for(auto &c: codes)
      if(c.code == id)
        type = c.pos;
  }

  ConfigAxes::Legende::Position::operator string() const
  {
    if(type == NORD_EST)
      return "ne";
    else if(type == SUD_EST)
      return "se";
    else if(type == NORD_OUEST)
      return "no";
    else if(type == SUD_OUEST)
      return "so";
    return "?";
  }


  void Axe::configure(const ConfigAxe &config, int dim)
  {
    DBG(msg("Configure axe : vmin={}, vmax={}", config.vmin, config.vmax);)

    if(config.vmin > config.vmax)
      echec("Configuration axe : vmin ({}) > vmax ({})", config.vmin, config.vmax);

    this->config = config;
    this->dim    = dim;
    decalage_zero   = (config.vmax + config.vmin) / 2;
    float diff = config.vmax - config.vmin;
    if(diff == 0)
      diff = 1;
    pixels_par_lsb = dim / diff;

    if(config.echelle_logarithmique)
    {
      float diff_db  = (log10(config.vmax) - log10(config.vmin));
      pixels_par_lsb = dim / diff_db;
      decalage_zero  = (log10(config.vmax) + log10(config.vmin)) / 2;
    }
  }


  bool Axe::tic_ok(float v)
  {
    int y = valeur_vers_pixel(v);
    if(y < 5)
      return false;
    if(y >= ((int) dim) - 5)
      return false;
    return true;
  }

  bool Axe::tic_ok_souple(float v)
  {
    int y = valeur_vers_pixel(v);
    if(y < 0)
      return false;
    if(y >= ((int) dim))
      return false;
    return true;
  }


  int Axe::valeur_vers_pixel(float v) const
  {
    if(config.echelle_logarithmique)
    {
      // vmin et vmax sont en linéaire, v aussi
      if(v == 0)
        return -1;
      v = log10(v);
    }
    int x = (int) ((v  - decalage_zero) * pixels_par_lsb);
    if(vertical)
      return dim/2-x;
    else
      return dim/2+x;
  }

  float Axe::valeur_vers_pixelf(float v) const
  {
    if(config.echelle_logarithmique)
    {
      // vmin et vmax sont en linéaire, v aussi
      v = log10(v);
    }
    float x = (float) ((v  - decalage_zero) * pixels_par_lsb);
    if(vertical)
      return ((float) dim)/2-x;
    else
      return ((float) dim)/2+x;
  }

  float Axe::pixel_vers_valeur(float p) const
  {
    /*if(vertical)
      return config.vmax - (p/dim) * (config.vmax - config.vmin);
    else
      return (p/dim) * (config.vmax - config.vmin);*/

    float tmp;
    if(vertical)
    {
      tmp = dim/2-p;
    }
    else
    {
      tmp = p-dim/2;
    }
    auto res = tmp / pixels_par_lsb + decalage_zero;

    if(config.echelle_logarithmique)
      return pow(10.0f, res);

    return res;
  }



// Ex: 100 MHz
// pnorm = 100 MHz
// Calcul une valeur de tic vérifiant tic >= tic_min,
// et "simple" en terme d'écriture décimale
static double calcule_tic(float tic_min)
{
  double res = 1;
  double pnorm = pow(10.0, (double) floor(log10(tic_min)));

  float v = tic_min / pnorm;
  // 0.15  -> 1.5  -> tic = 0.2
  // 23.15 -> 2.31 -> tic = 25



  if(v == 1)
  {
    // Dans ce cas, tic_min est déjà une puissance de 10
    res = 1;
  }
  else if(v <= 2)
    res = 2;
  else if(v <= 2.5)
    res = 2.5;
  else if(v <= 5.0)
    res = 5.0;
  else
    res = 10.0;

  res *= pnorm;

  //msg("tic min = {} -> tic = {} (pnorm = {})", tic_min, res, pnorm);

  return res;
}





struct Axes::Impl
{
  // TODO: regrouper les mutable
  mutable ConfigAxes config;

  // Position horizontale de l'axe vertical
  mutable int x_grad_vert = 0, y_grad_hor = 0;

  mutable double tic_majeur = 0, tic_mineur = 0;
  mutable unsigned int tic_majeur_pixels = 0, tic_mineur_pixels = 0;
  mutable Axe axe_horizontal, axe_vertical;
  /** Dimension de la vue, en pixels */
  mutable int sx = 1, sy = 1;
  /** Dimension allouée pour les courbes */
  mutable int sx_graphique = 1, sy_graphique = 1;
  mutable int marge_gauche = 0, marge_droite = 0;
  mutable int marge_basse = 0, marge_haute = 0;
  mutable   bool    mat_sur_clair = false;
  mutable Couleur couleur_avant_plan = {255,255,255};
  mutable Couleur fond_grille;

  //mutable string ech_x, ech_y;


  void post_rendu(Canva canva_) const
  {
    if(config.legende.afficher)
    {
      Canva canva = canva_.vue(Rect{0, 0, sx, sy});
      affiche_cartouche_legende(canva);
    }
  }

  Dim get_dim_graphique(const Dim &dim_totale) const
  {
    maj_dimensions(dim_totale);
    return {sx_graphique, sy_graphique};
  }







  // dim = allocation
  void maj_dimensions(const Dim &dim) const
  {
    auto &cah = config.axe_horizontal, &cav = config.axe_vertical;
    sx = dim.l;
    sy = dim.h;
    marge_haute = marge_basse  = marge_gauche = marge_droite = 0;

    if(!config.titre.empty())
    {
      if(0.05 * sy < HAUTEUR_TITRE)
      {
        marge_haute = HAUTEUR_TITRE - 0.05 * sy;
      }
    }

    if(cah.afficher)
    {
      marge_basse  = HAUTEUR_X_VALEURS + 2;

      if(!cah.label.empty())
        marge_basse += HAUTEUR_X_LABEL;
    }

    if(cav.afficher)
    {
      if(config.position_axe_vertical == ConfigAxes::AXE_VERTICAL_GAUCHE)
      {
        marge_gauche = 60;
        if(!cav.label.empty())
          marge_gauche += 25;
      }
      else
      {
        marge_droite = 60;
        if(!cav.label.empty())
          marge_droite += 25;
      }
    }

    int marge = marge_gauche + marge_droite;
    int margev = marge_basse + marge_haute;
    sx_graphique = sx > marge ? sx - marge : sx;
    sy_graphique = sy > margev ? sy - margev : sy;

    x_grad_vert = marge_gauche;

    if(config.position_axe_vertical == ConfigAxes::AXE_VERTICAL_DROITE)
    {
      x_grad_vert = sx - marge_droite;
    }

    y_grad_hor = sy - marge_basse;

    if(config.grille_majeure.pixels_par_tic_min == -1)
      config.grille_majeure.pixels_par_tic_min = 70;
    if(config.grille_mineure.pixels_par_tic_min == -1)
      config.grille_mineure.pixels_par_tic_min = 4;

    cah.id = "X";
    cav.id = "Y";
  }

  mutable Canva canva_pixels, canva_res;

  Canva rendre(Canva canva_, const Rectf &rdi) const
  {
    return rendre(canva_, rdi.x, rdi.x + rdi.l, rdi.y + rdi.h, rdi.y);// !!
  }

  Canva rendre(Canva canva_, float xmin, float xmax, float ymin, float ymax) const
  {
      // Ou alors, tricher dans le canva :
      // get_allocation renvoie la dim du sous-espace
      auto dim = canva_.get_allocation(); // Allocation en pixels

      //canva_.set_rdi({0, 0, (float) dim.l, (float) dim.h});

      maj_dimensions(dim);

      //DBG(msg("  rendu axes : allocation = {}x{}, rdi = {}", sx, sy, rdi));

      // Vue en pixels
      canva_pixels = canva_.vue(Rect{0, 0, sx, sy});

      // Est-ce qu'on a besoin de ça sachant que maintenant ces infos sont dans le canva ?

      auto &cah = config.axe_horizontal, &cav = config.axe_vertical;

      /*cah.vmin = rdi.x;
      cah.vmax = rdi.x + rdi.l;
      // Inversion
      cav.vmax = rdi.y;
      cav.vmin = rdi.y + rdi.h;*/

      cah.vmin = xmin;
      cah.vmax = xmax;
      cav.vmax = ymax;
      cav.vmin = ymin;




      DBG(msg("rendu axes : hmin={}, hmax={}, vmin={}, vmax={}, ecart pix={}, sx={}, sy={}",
          cah.vmin, cah.vmax,
          cav.vmin, cav.vmax,
          config.axe_horizontal.ecart_pixel,
          sx, sy);)

      // PB ICI SI ECHELLE LOG
      //canva_res = canva_pixels.clip(Rect{0, marge_haute, sx_graphique, sy_graphique}, rdi,
      //    config.axe_horizontal.echelle_logarithmique, config.axe_vertical.echelle_logarithmique);



      if(config.mode_isoview)
      {
        // Il faut se débrouiller pour que la résolution soit identique
        // en x et en y. Soit LSB / pixel identique.
        float dh = cah.vmax - cah.vmin, dv = cav.vmax - cav.vmin;
        float cth = (cah.vmax + cah.vmin) / 2, ctv = (cav.vmax + cav.vmin) / 2;

        // Résolution = nb lsb / pixel, en horizontal et en vertical
        float reso_x = dh / sx_graphique, reso_y = dv / sy_graphique;
        float ratio = (reso_x + 1e-200) / (reso_y + 1e-200);

        // Comment forcer les résolutions à être identique ?
        // On augmente l'intervalle pour la résolution la plus haute
        auto ar = config.aspect_ratio;

        // ar = ratio résolution y / x

        DBG(msg("RDI avant AR : extension_X/extension_Y = {} VS ratio d'aspet demandé = {}", ratio, ar);)
        if(ratio <  1.0 / ar)
        {
          //auto delta = (dh/2) *  ratio * ar;
          auto δ = sx_graphique * reso_y / (2 * ar);
          cah.vmin = cth - δ;
          cah.vmax = cth + δ;
          DBG(float nv_reso_x = 2 * δ / sx_graphique;)
          DBG(msg("cah : δ=+/- {}, nv ratio = {}", δ, nv_reso_x / reso_y);)
          // 2 (delta / sxg) / reso_y = 1/ar
          // delta = sxg * resoy / (2 * ar)
        }
        else
        {
          //auto delta = (dv/2) * ratio / ar;
          // nv_reso_y = 2 * delta / syg
          // resox / nv_reso_y  = 1/ar
          // nv_reso_y = ar/reso_x
          // delta = syg * ar / (2 * reso_x)
          auto δ = sy_graphique * reso_x / (2 * ar);
          cav.vmin = ctv - δ;
          cav.vmax = ctv + δ;
          DBG(auto nv_reso_y = 2 * δ / sy_graphique;)
          DBG(msg("cav : δ=+/- {}, ratio d'aspet final = {}", δ, reso_x / nv_reso_y);)
        }
      }


      xmin = cah.vmin;
      xmax = cah.vmax;
      ymax = cav.vmax;
      ymin = cav.vmin;

      Rectf rdi{xmin,ymax, xmax - xmin, ymin - ymax};
      canva_res = canva_pixels.clip_alt(
          Rect{0, marge_haute, sx_graphique, sy_graphique}, xmin, xmax, ymax, ymin,
          config.axe_horizontal.echelle_logarithmique, config.axe_vertical.echelle_logarithmique);


      DBG(msg("raxes: cfg axes..."));
      axe_horizontal.configure(cah, sx_graphique);
      axe_vertical.configure(cav, sy_graphique);

      DBG(msg("raxes: cpt tics..."));
      calcule_tics(axe_vertical);
      DBG(msg("Calcul tics horizontaux...");)
      calcule_tics(axe_horizontal);
      DBG(msg("ok.");)

      DBG(msg("Axes / couleurs :"));
      DBG(msg(" - ARP    = {}", config.couleur_arriere_plan));
      DBG(msg(" - grille (majeure) = {}", config.grille_majeure.couleur));

      mat_sur_clair = config.est_mat_sur_clair();
      couleur_avant_plan = mat_sur_clair ? Couleur::Noir : Couleur::Blanc;

      DBG(msg(" - mat sur clair = {}", mat_sur_clair));
      DBG(msg(" - => avant plan = {}", couleur_avant_plan));

      if(config.grille_majeure.afficher == -1)
        config.grille_majeure.afficher = 1;

      if(config.grille_mineure.afficher == -1)
        config.grille_mineure.afficher = 1;//0;

      if(config.grille_majeure.lg_tiret == -1)
      {
        config.grille_majeure.lg_tiret = 10;
        config.grille_majeure.lg_trou  = 0;
      }
      if(config.grille_mineure.lg_tiret == -1)
      {
        config.grille_mineure.lg_tiret = 1;
        config.grille_mineure.lg_trou  = 3;
      }

      if(config.grille_majeure.couleur.vers_rgba() == 0)
        config.grille_majeure.couleur = affaiblir_couleur(couleur_avant_plan, 0.3);

      if(config.grille_mineure.couleur.vers_rgba() == 0)
       config.grille_mineure.couleur = affaiblir_couleur(config.grille_majeure.couleur, 0.3);


      DBG(msg(" - grille majeure = {}", config.grille_majeure.couleur));
      DBG(msg(" - grille mineure = {}", config.grille_mineure.couleur));

      DBG(msg("raxes: ok, grille..."));

      //canva.set_coord_mode_pixels(true);

      affiche_grille(canva_pixels);

      DBG(msg("raxes: ok, titre..."));

      //float hauteur_max = 20;//sy_graphique * (0.25 / (1 + 0.25 + 0.1));


      canva_pixels.set_align(Align::CENTRE, Align::DEBUT);
      //canva.set_dim_fonte(config.titre_dim/10);
      canva_pixels.set_dim_fonte(FONTE_TITRE);
      canva_pixels.set_epaisseur(1);
      canva_pixels.set_couleur(couleur_avant_plan);
      canva_pixels.texte(0.5*sx_graphique, 0, config.titre, sx_graphique, HAUTEUR_TITRE);

      DBG(msg(" raxes: fin."));
      return canva_res;
    }

# if 0
  // Attention, dans quelle échelle il doit rendre ?
  Canva rendre(Canva canva_, const Rectf &rdi) const
  {
    // Ou alors, tricher dans le canva :
    // get_allocation renvoie la dim du sous-espace
    auto dim = canva_.get_allocation(); // Allocation en pixels


    maj_dimensions(dim);

    DBG(msg("  rendu axes : allocation = {}x{}, rdi = {}", sx, sy, rdi));

    // Vue en pixels
    canva_pixels = canva_.vue(Rect{0, 0, sx, sy});

    // Est-ce qu'on a besoin de ça sachant que maintenant ces infos sont dans le canva ?

    auto &cah = config.axe_horizontal, &cav = config.axe_vertical;

    cah.vmin = rdi.x;
    cah.vmax = rdi.x + rdi.l;
    // Inversion
    cav.vmax = rdi.y;
    cav.vmin = rdi.y + rdi.h;

    DBG(msg("rendu axes : hmin={}, hmax={}, vmin={}, vmax={}, ecart pix={}, sx={}, sy={}",
        cah.vmin, cah.vmax,
        cav.vmin, cav.vmax,
        config.axe_horizontal.ecart_pixel,
        sx, sy);)


    canva_res = canva_pixels.clip(Rect{0, marge_haute, sx_graphique, sy_graphique}, rdi,
        config.axe_horizontal.echelle_logarithmique, config.axe_vertical.echelle_logarithmique);

    if(config.mode_isoview)
    {
      // Il faut se débrouiller pour que la résolution soit identique
      // en x et en y. Soit LSB / pixel identique.
      float dh = cah.vmax - cah.vmin, dv = cav.vmax - cav.vmin;
      float cth = (cah.vmax + cah.vmin) / 2, ctv = (cav.vmax + cav.vmin) / 2;
      float reso_x = dh / sx_graphique, reso_y = dv / sy_graphique;
      float ratio = (reso_x + 1e-200) / (reso_y + 1e-200);

      // Comment forcer les résolutions à être identique ?
      // On augmente l'intervalle pour la résolution la plus haute
      auto ar = config.aspect_ratio;

      // ar = ratio résolution y / x

      DBG(msg("Ratio res = {}, ar = {}", ratio, ar);)
      if(ratio <  1.0 / ar)
      {
        //auto delta = (dh/2) *  ratio * ar;
        auto δ = sx_graphique * reso_y / (2 * ar);
        cah.vmin = cth - δ;
        cah.vmax = cth + δ;
        DBG(float nv_reso_x = 2 * δ / sx_graphique;)
        DBG(msg("cah : δ=+/- {}, nv ratio = {}", δ, nv_reso_x / reso_y);)
        // 2 (delta / sxg) / reso_y = 1/ar
        // delta = sxg * resoy / (2 * ar)
      }
      else
      {
        //auto delta = (dv/2) * ratio / ar;
        // nv_reso_y = 2 * delta / syg
        // resox / nv_reso_y  = 1/ar
        // nv_reso_y = ar/reso_x
        // delta = syg * ar / (2 * reso_x)
        auto δ = sy_graphique * ar / (2 * reso_x);
        cav.vmin = ctv - δ;
        cav.vmax = ctv + δ;
        DBG(auto nv_reso_y = 2 * δ / sy_graphique;)
        DBG(msg("cav : δ=+/- {}, nv ratio = {}", δ, reso_x / nv_reso_y);)
      }
    }


    DBG(msg("raxes: cfg axes..."));
    axe_horizontal.configure(cah, sx_graphique);
    axe_vertical.configure(cav, sy_graphique);

    DBG(msg("raxes: cpt tics..."));
    calcule_tics(axe_vertical);
    DBG(msg("Calcul tics horizontaux...");)
    calcule_tics(axe_horizontal);
    DBG(msg("ok.");)

    DBG(msg("Axes / couleurs :"));
    DBG(msg(" - ARP    = {}", config.couleur_arriere_plan));
    DBG(msg(" - grille (majeure) = {}", config.grille_majeure.couleur));

    mat_sur_clair = config.est_mat_sur_clair();
    couleur_avant_plan = mat_sur_clair ? Couleur::Noir : Couleur::Blanc;

    DBG(msg(" - mat sur clair = {}", mat_sur_clair));
    DBG(msg(" - => avant plan = {}", couleur_avant_plan));

    if(config.grille_majeure.afficher == -1)
      config.grille_majeure.afficher = 1;

    if(config.grille_mineure.afficher == -1)
      config.grille_mineure.afficher = 1;//0;

    if(config.grille_majeure.lg_tiret == -1)
    {
      config.grille_majeure.lg_tiret = 10;
      config.grille_majeure.lg_trou  = 0;
    }
    if(config.grille_mineure.lg_tiret == -1)
    {
      config.grille_mineure.lg_tiret = 1;
      config.grille_mineure.lg_trou  = 3;
    }

    if(config.grille_majeure.couleur.vers_rgba() == 0)
      config.grille_majeure.couleur = affaiblir_couleur(couleur_avant_plan, 0.3);

    if(config.grille_mineure.couleur.vers_rgba() == 0)
     config.grille_mineure.couleur = affaiblir_couleur(config.grille_majeure.couleur, 0.3);


    DBG(msg(" - grille majeure = {}", config.grille_majeure.couleur));
    DBG(msg(" - grille mineure = {}", config.grille_mineure.couleur));

    DBG(msg("raxes: ok, grille..."));

    //canva.set_coord_mode_pixels(true);

    affiche_grille(canva_pixels);

    DBG(msg("raxes: ok, titre..."));

    //float hauteur_max = 20;//sy_graphique * (0.25 / (1 + 0.25 + 0.1));


    canva_pixels.set_align(Align::CENTRE, Align::DEBUT);
    //canva.set_dim_fonte(config.titre_dim/10);
    canva_pixels.set_dim_fonte(FONTE_TITRE);
    canva_pixels.set_epaisseur(1);
    canva_pixels.set_couleur(couleur_avant_plan);
    canva_pixels.texte(0.5*sx_graphique, 0, config.titre, sx_graphique, HAUTEUR_TITRE);

    DBG(msg(" raxes: fin."));
    return canva_res;
  }
# endif

  // Quel est le principe :
  //  - on connait vmax, vmin
  //  -
  void calcule_tics(Axe &axe) const
  {

    if(axe.config.tic_manuels)
    {
      axe.tics_majeurs_pos.clear();
      axe.tics_mineurs_pos.clear();

      for(auto p: axe.config.tics)
        axe.tics_majeurs_pos.push_back(p.first);

      return;
    }

    bool mode_log = axe.config.echelle_logarithmique;

    if(axe.config.ecart_pixel != -1)
    {
      DBG(msg("Ecart pixel = {}", axe.config.ecart_pixel);)
      axe.config.vmin = - (axe.config.ecart_valeur * axe.dim) / axe.config.ecart_pixel;
      axe.config.vmax = 0;
      // On doit reconfigurer pour pplsb, etc.
      axe.configure(axe.config, axe.dim);
    }

    double vmin = axe.config.vmin, vmax = axe.config.vmax;

    if(mode_log)
    {
      DBG(msg("mode log: vmin {} -> {}, vmax {} -> {}", vmin, log10(vmin), vmax, log10(vmax)));
      vmin = log10(vmin);
      vmax = log10(vmax);
    }

    DBG(msg("Axe [{}]: dim={}, cfg.min={}, cfg.max={}, pixels_min={}, ecart-pixel={}, mode={}",
          axe.config.id, axe.dim, vmin, vmax,
          config.grille_majeure.pixels_par_tic_min,
          axe.config.ecart_pixel,
          mode_log ? "log" : "lin");)

    if(vmin > vmax)
    {
      msg_erreur("Vmin ({}) > vmax ({}) !", vmin, vmax);
    }

    // Calcul des tics majeurs et mineurs

    // δ = diff vmax - vmin, en unité utilisateur
    float δ = abs(vmax - vmin);

    // exemple : δ = 15
    // => unité = 1

    // => On cherche le plus petit tic majeur < au plus petit
    //     de abs(vmax) / abs(vmin)

    float pixels_par_tic_majeur_min = config.grille_majeure.pixels_par_tic_min;

    essai:

    // Imposé par la résolution
    // eg lsb_par_pixel / pixel_par_tic_majeur_min

    // Exemple : δ = 1000, axe.dim = 100, pptm_min = 50
    // tic_majeur_min = 10 * 50 = 500 LSB

    // Tic majeur mini, en unité LSB (ou LSB_LOG en mode logarithmique)
    // Au moins ... pixels => au moins ... en valeur utilisateur
    float tic_majeur_min = (δ / axe.dim) * pixels_par_tic_majeur_min;
    float tic_mineur_min = (δ / axe.dim) * config.grille_mineure.pixels_par_tic_min;




    tic_majeur = calcule_tic(tic_majeur_min);

    DBG(msg("  δ={}, tcmin={} -> essai tic majeur = {}", δ, tic_majeur_min, tic_majeur);)

    // En mode logarithmique, travaille au mimimun par puissance entière de 10
    if(mode_log)
      tic_majeur = ceil(tic_majeur);

    if(axe.config.ecart_pixel != -1)
      tic_majeur = axe.config.ecart_pixel * (δ / axe.dim);

    if((tic_majeur == 0) || (isnan(tic_majeur)))
    {
      msg_avert("Tic majeur : {}", tic_majeur);
      return;
    }

    if(tic_majeur <= numeric_limits<double>::epsilon())
    {
      tic_majeur = 10 * numeric_limits<double>::epsilon();
      DBG(msg_avert("tcm < eps -> {}...", tic_majeur);)
    }

    axe.ecart_tics_majeurs_valeur = tic_majeur;
    axe.ecart_tics_majeurs_pixels = tic_majeur * (axe.dim / δ);

    DBG(msg("  Essai tic majeur = {}, pp_par_tic_majeur_min = {} -> ecart_tic_majeurs_pix = {}...", tic_majeur, pixels_par_tic_majeur_min, axe.ecart_tics_majeurs_pixels);)

    tic_mineur = tic_majeur / config.grille_majeure.sous_graduation;
    axe.tics_mineurs_presents = tic_mineur >= tic_mineur_min;

    axe.tics_majeurs_pos.clear();
    axe.tics_mineurs_pos.clear();

    // Entre tic_majeur_min et δ, on cherche des valeurs
    // pas trop bêtes

    bool zero_present = (vmin <= 0) && (vmax >= 0);

    // Les variables tk, tic_majeur, etc. sont en logs
    auto adapt = [&](double v)
    {
      if(mode_log)
        return pow(10.0f, v);
      return v;
    };

    // Ajoute des tics mineurs au dessus de tk (tic majeur)s
    auto ajouter_tics_mineurs = [&](float tk)
    {
      if(!axe.tics_mineurs_presents)
        return;

      // En mode logarithmique :
      // tic_mineur, tic_majeur sont en logarithmique
      // d'ailleurs, tic_mineur n'a pas de sens
      // (    tic_mineur = tic_majeur / config.grille_majeure.sous_graduation;)
      // On commence par tk
      // tk est en logarithmique.
      // Par exemple tk = -10
      // disons tic_majeur = 10
      // On va y aller de manière linéaire et pas log (pas de 2 en 2 en log)
      // 10^tk, 10^tk - R, 10^tk - 2*R, ... , [10^(tk-ticmajeur)
      // => en log : log(10^tk - R), log(10^tk - 2*R), etc...
      // R = ?
      // Exemple : tic_majeur = 5, 5 pas
      // Intervalle linéaire : 10^5
      // à couper en 5 parties :
      // 0, [10^5/5, 2*10^5/5, 3*10^5/5, 4*10^5/5,] 5*10^5/5
      // => R = 10^(tic_majeur) / nb_grad

      if(mode_log)
      {
        // si sous grad =  1/10:
        // R = 10^(tk+M) / 10 = 10^{tk+M-1}
        // tk0 = log10(10^tk + 10^{tk+M-1})
        //       log10(10^tk * (1 + 10^{M-1})
        //       tk + log10(1 + 10^{M-1})
        //     ~ tk + M - 1
        // si M = 5
        //
        // Exemple : tic_majeur = 1, tk = 1 (10^1 -> 10^2)
        // R = 100 / 10 = 10;
        // tkm0 = log10(20) = 1.3



        double R = pow(10.0f, tk + tic_majeur) / config.grille_majeure.sous_graduation;
        int k = 1;
        double tkm = log10(pow(10.0, tk) + k * R);
        while(tkm < tk + tic_majeur)
        {
          if(axe.tic_ok(adapt(tkm)))
            axe.tics_mineurs_pos.push_back(adapt(tkm));
          tkm = log10(pow(10.0, tk) + k * R);
          k++;
        }
      }
      else
      {
        double tkm = tk + tic_mineur;
        while(tkm < tk + tic_majeur)
        {
          if(axe.tic_ok(adapt(tkm)))
            axe.tics_mineurs_pos.push_back(adapt(tkm));
          tkm += tic_mineur;
        }
      }
    };

    if(zero_present)
    {
      // => On a forcément un tic à zéro
      double rmax = max(-vmin, vmax);
      double tk = 0;//tic_majeur;
      while(tk < rmax)
      {
        if(axe.tic_ok(adapt(tk)))
          axe.tics_majeurs_pos.push_back(adapt(tk));
        if((tk != 0) && axe.tic_ok(adapt(-tk)))
          axe.tics_majeurs_pos.push_back(adapt(-tk));

        ajouter_tics_mineurs(tk);
        ajouter_tics_mineurs(-tk-tic_majeur);

        tk += tic_majeur;
      }
    }
    else
    {
      // Dans le cas contraire, quel est le point de référence ???
      double tk;

      // vmax > vmin > 0
      if(vmin >= 0)
      {
        // => Le plus petit multiple de tic_majeur > config.min
        DBG(msg("vmin={}, tmaj={} --> ceil={}", vmin, tic_majeur, tic_majeur * ceil(vmin / tic_majeur));)
        tk = tic_majeur * ceil(vmin / tic_majeur);

        if(!axe.tic_ok(adapt(tk)))
          tk += tic_majeur;

        DBG(msg(" tk={}, tic ok = {}", tk, axe.tic_ok(adapt(tk)));)
        DBG(msg("valeur_vers_pixel : {}", axe.valeur_vers_pixel(tk));)

        while(axe.tic_ok_souple(adapt(tk)))
        {
          if(axe.tic_ok(adapt(tk)))
            axe.tics_majeurs_pos.push_back(adapt(tk));
          ajouter_tics_mineurs(tk);
          tk += tic_majeur;
        }
      }
      // vmin < vmax < 0
      else
      {
        // => Le plus grand multiple de tic_majeur < vmax
        tk = -tic_majeur * ceil(-vmax / tic_majeur);

        DBG(msg("  -> essai premier tic = {} LSB.", tk));

        while(axe.tic_ok_souple(adapt(tk)))
        {
          if(axe.tic_ok(adapt(tk)))
            axe.tics_majeurs_pos.push_back(adapt(tk));
          ajouter_tics_mineurs(tk);
          tk -= tic_majeur;
        }
      }
    }

    //DBG(msg("  Tics majeurs :");
    //for(auto &p: axe.tics_majeurs_pos)
      //msg("   {:.5f} (log10 = {:.2f})", p, log10(p)););


    if(axe.config.ecart_pixel == -1)
    {
      if((axe.tics_majeurs_pos.size() <= 1) && (pixels_par_tic_majeur_min > 1))
      {
        pixels_par_tic_majeur_min = floor(pixels_par_tic_majeur_min * 0.8f);
        if(pixels_par_tic_majeur_min >= 2)
        {
          DBG(msg("Calcul tics : aucun tic majeur trouvé, ré-essai en réduisant l'écart / pixel..."));
          goto essai;
        }
        else
        {
          DBG(msg_avert("Calcul tics : aucun tic majeur trouvé -> abandon."));
        }
      }
    }


    DBG(msg("Calcule tics : {} tics mineurs, {} majeurs.",
        axe.tics_mineurs_pos.size(), axe.tics_majeurs_pos.size());)

  }




  Couleur affaiblir_couleur(const Couleur &c, float ratio = 0.5) const
  {
    Couleur res;
    if(mat_sur_clair)
      return c.eclaircir(ratio);
    return c.assombrir(ratio);
  }

  Couleur renforcer_couleur(const Couleur &c, float ratio = 0.5) const
  {
    if(mat_sur_clair)
      return c.assombrir(ratio);
    return c.eclaircir(ratio);
  }

  float pixel_vers_x(float x) const
  {
    return axe_horizontal.pixel_vers_valeur(x);
  }

  float pixel_vers_y(float y) const
  {
    return axe_vertical.pixel_vers_valeur(y);
  }

  Point pos(float x, float y) const
  {
    return Point{x_vers_pixel(x), y_vers_pixel(y)};
  }

  tuple<float, float> posf(float x, float y) const
  {
    return {axe_horizontal.valeur_vers_pixelf(x),
            axe_vertical.valeur_vers_pixelf(y)};
  }

  int x_vers_pixel(float v) const
  {
    return canva_res.v2c(Pointf{v,0.0f}).x;
    //return axe_horizontal.valeur_vers_pixel(v);
  }

  int y_vers_pixel(float v) const
  {
    return canva_res.v2c(Pointf{0.0f,v}).y;
    //return axe_vertical.valeur_vers_pixel(v);
  }


  void affiche_axe_vertical(Canva O) const
  {
    if(!config.axe_vertical.afficher)
      return;

    int x0 = x_grad_vert;

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    O.ligne(x0,marge_haute,x0,marge_haute+sy_graphique-1);
    O.set_dim_fonte(config.axe_vertical.valeurs_echelle);

    auto [expo, nbchiffres] = tsd::vue::unites::calc_expo_nb_chiffres_commun(axe_vertical.tics_majeurs_pos, config.axe_vertical.unite);

    if(axe_vertical.tics_majeurs_pos.empty())
      msg_avert("Axes : aucun tic vertical !");

    for(auto t: axe_vertical.tics_majeurs_pos)
    {
      string s;
      if(axe_vertical.config.echelle_logarithmique)
        s = fmt::format("{:.1e}", t);
      else
        s = tsd::vue::unites::valeur_vers_chaine(t, config.axe_vertical.unite, expo, nbchiffres);
      affiche_tic_majeur(O, t, 5, x0, s);
    }
    int lg = 1;
    for(auto t: axe_vertical.tics_mineurs_pos)
      affiche_tic_mineur(O, t, lg, x0);

    O.set_dim_fonte(FONTE_X_LABEL);
    if((!config.axe_vertical.label.empty()) && (sy > 20) && (sx > 100))
    {
      auto label_y = config.axe_vertical.label;
      if(label_y.size() > 1)
        O.set_orientation(Orientation::VERTICALE);
      O.texte(sx - 25, sy/2, label_y, 25, sy);
      O.set_orientation(Orientation::HORIZONTALE);
    }
  }


  // Ajoute un plan sur l'image en cours.
  // Ne modifie que les pixels qui sont encore à la couleur de fond
  void affiche_grille(Canva canva) const
  {
    DBG(msg(" raxes/affiche_grille : lignes h"));
    if(config.grille_majeure.afficher)
    {
      if(axe_vertical.tics_mineurs_presents && config.grille_mineure.afficher)
      {
        for(auto &t: axe_vertical.tics_mineurs_pos)
        {
          int y = canva_res.v2c(Pointf{0.0f,(float) t}).y;//y_vers_pixel(t);
          ligne_stylee_horizontale(canva, {marge_gauche, y},
                                  {sx-1-marge_droite,y},
                                   config.grille_mineure);
        }
    }

    DBG(msg(" raxes/affiche_grille : lignes v"));
    if(axe_horizontal.tics_mineurs_presents && config.grille_mineure.afficher)
    {
      for(auto &t: axe_horizontal.tics_mineurs_pos)
      {
        int x = canva_res.v2c(Pointf{(float) t,0.0f}).x;
        //int x = x_vers_pixel(t);
        ligne_stylee_verticale(canva,
                               {marge_gauche+x, marge_haute},
                               {marge_gauche+x,sy-1-marge_basse},
                               config.grille_mineure);
      }
    }

    DBG(msg(" raxes/affiche_grille : lignes V"));
    for(auto t: axe_horizontal.tics_majeurs_pos)
    {
      int x = canva_res.v2c(Pointf{(float) t,0.0f}).x;
      //int x = x_vers_pixel(t);
      ligne_stylee_verticale(canva,
                             {marge_gauche+x, marge_haute},
                             {marge_gauche+x, sy-1-marge_basse},
                             config.grille_majeure);
    }

    DBG(msg(" raxes/affiche_grille : lignes H"));
    for(auto &t: axe_vertical.tics_majeurs_pos)
    {
      int y = canva_res.v2c(Pointf{0.0f,(float)t}).y;
      //int y = y_vers_pixel(t);
      ligne_stylee_horizontale(canva, {marge_gauche, y},
                               {sx-1-marge_droite,y},
                               config.grille_majeure);
    }
    }

    DBG(msg(" raxes/affiche_grille : axes"));
    affiche_axe_vertical(canva);
    affiche_axe_horizontal(canva);
    DBG(msg(" raxes/affiche_grille : fait"));


  }

  void affiche_axe_horizontal(Canva O) const
  {
    if(!config.axe_horizontal.afficher)
      return;

    int y0 = y_grad_hor;

    DBG(msg(" affichage axe horizontal : sy={}, y_grad_hor={}", sy, y_grad_hor);)

    O.set_couleur(couleur_avant_plan);
    O.ligne(Point{marge_gauche, y0}, Point{marge_gauche+sx_graphique-1, y0});

    auto [expo, nbchiffres] = tsd::vue::unites::calc_expo_nb_chiffres_commun(axe_horizontal.tics_majeurs_pos, config.axe_horizontal.unite);

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    O.set_align(Align::CENTRE, Align::DEBUT);

    O.set_dim_fonte(config.axe_horizontal.valeurs_echelle);

    int idx = 0;
    for(auto t: axe_horizontal.tics_majeurs_pos)
    {
      int ti = x_vers_pixel(t);
      O.ligne(Point{ti+marge_gauche, y0}, Point{ti+marge_gauche, y0 - 5});
      auto s = tsd::vue::unites::valeur_vers_chaine(t, config.axe_horizontal.unite, expo, nbchiffres);

      if(config.axe_horizontal.tic_manuels && (idx < (int) config.axe_horizontal.tics.size()))
        s = config.axe_horizontal.tics[idx].second;

      O.texte(marge_gauche+ti,
              y0 + ECART_X_VALEURS, s,
              axe_horizontal.ecart_tics_majeurs_pixels,
              HAUTEUR_X_VALEURS - ECART_X_VALEURS);

      idx++;
    }

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    for(auto t: axe_horizontal.tics_mineurs_pos)
    {
      int x = x_vers_pixel(t);
      if(y0 < (int) sy/2)
        O.ligne(Point{x, y0+1}, Point{x, y0});
      else
        O.ligne(Point{x, y0-1}, Point{x, y0});
    }

    if((!config.axe_horizontal.label.empty()) && (sy > 20) && (sx > 100))
    {
      O.set_dim_fonte(FONTE_X_LABEL);
      O.texte(sx/2, sy-HAUTEUR_X_LABEL-1, config.axe_horizontal.label, 0.8*sx_graphique, HAUTEUR_X_LABEL);
    }
  }

  void affiche_tic_majeur(Canva &O, float v, int lg, int x0, const string &val) const
  {
    int y = y_vers_pixel(v);

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    if(x0 < (int) sx/2)
      O.ligne(Point{x0+lg, y}, Point{x0, y});
    else
      O.ligne(Point{x0-lg, y}, Point{x0, y});

    //auto s = valeur_vers_chaine(v, config.axe_vertical.unite);

    int pos = x0 + 6;
    Align align = Align::DEBUT;
    if(x0 < (int) sx/2)
    {
      pos = 20;//4;
      align = Align::FIN;
    }

    O.set_align(align, Align::CENTRE);

    float dxmax = 60 - 6;
    float dymax = HAUTEUR_Y_VALEUR;
    O.texte(pos, y, val, dxmax, dymax);

    /*msg("Affiche label tic majeur : {}, @[{},{}], dim=[{},{}]", val, pos, y, dxmax, dymax);
    if(pos + dxmax == sx)
      msg_avert("Presque dépassement x !");
    if(y + dymax == sy)
      msg_avert("Presque dépassement y !");
    if(pos + dxmax > sx)
      msg_avert("Dépassement x !");
    if(y + dymax > sy)
      msg_avert("Dépassement y !");*/

  }


  void affiche_tic_mineur(Canva &O, float v, int lg = 2, int x0 = 20) const
  {
    int y = y_vers_pixel(v);
    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    if(x0 < (int) sx/2)
      O.ligne(Point{x0+lg, y}, Point{x0, y});
    else
      O.ligne(Point{x0-lg, y}, Point{x0, y});
  }

  void ligne_stylee_horizontale(Canva canva,
      const Point p0, const Point p1,
      const ConfigGrille &style) const
  {
    canva.set_couleur(style.couleur);
    canva.set_epaisseur(style.epaisseur);
    canva.set_dotted(style.lg_trou != 0);
    canva.ligne(p0.x, p0.y, p1.x, p0.y);
    canva.set_dotted(false);
    //for(auto x = p0.x; x < (int) (p1.x-style.lg_tiret); x += style.lg_trou + style.lg_tiret)
    //  canva.ligne(x, p0.y, x+style.lg_tiret, p0.y);
  }

  void ligne_stylee_verticale(Canva canva,
      const Point p0, const Point p1,
      const ConfigGrille &style) const
  {
    canva.set_couleur(style.couleur);
    canva.set_epaisseur(style.epaisseur);
    canva.set_dotted(style.lg_trou != 0);
    canva.ligne(p0.x, p0.y, p0.x, p1.y);
    canva.set_dotted(false);
    //for(auto y = p0.y; y < (int) (p1.y-style.lg_tiret); y += style.lg_trou + style.lg_tiret)
    //  canva.ligne(p0.x, y, p0.x, y+style.lg_tiret);
  }


  void affiche_cartouche_legende(Canva &canva) const
  {
    DBG(msg("raxes: legende..."));
    if((sx_graphique < 100) || (sy_graphique < 100))
      // Pas assez de place
      return;

    bool a_description = false;
    for(auto c: config.series)
    {
      if(c.nom.size() > 0)
      {
        a_description = true;
        break;
      }
    }

    if(!a_description)
    {
      DBG(msg("raxes: pas série configurée --> pas de légende."));
      return;
    }

    int n = config.series.size();

    Image tmp;


    float dim_texte = config.legende.dim;

    DBG(msg("Calcul dim max légende : échelle = {}", dim_texte));

    // Calcul dimension maximale des différents textes
    int lx_max = 0, ly_max = 0;
    for(auto i = 0; i < n; i++)
    {
      auto dim = tmp.texte_dim(config.series[i].nom, dim_texte);
      DBG(msg("dim[{}] : {}", i, dim);)
      lx_max = max(dim.l, lx_max);
      ly_max = max(dim.h, ly_max);
    }

    // 187x28
    // 165x22
    //

    int hauteur_max     = sy_graphique / 2;
    int hauteur_totale  = (ly_max + 2) * n;

    int largeur_max     = sx_graphique / 3;
    int largeur_totale  = lx_max + 50;

    DBG(msg("l = {}, lmax = {}, h = {}, hmax = {}", largeur_totale, largeur_max, hauteur_totale, hauteur_max));

    largeur_totale = min(largeur_max, largeur_totale);

    if(hauteur_totale > hauteur_max)
    {
      // Diminue l'échelle de manière à avoir la même échelle finale pour tous les éléments
      dim_texte *= (1.0f * hauteur_max) / hauteur_totale;
      hauteur_totale = hauteur_max;
    }

    DBG(msg("Après seuillage : largeur totale = {}, hauteur totale = {}", largeur_totale, hauteur_totale));

    lx_max = largeur_totale - 50;
    int h1 = (int) ceil((1.0f * hauteur_totale) / n);
    hauteur_totale = h1 * n;

    // -> h1 = 30

    DBG(msg("Après arrondi : largeur totale = {}, hauteur totale = {}", largeur_totale, hauteur_totale));
    DBG(msg("h1 = {}", h1));

    // Normalement, h1 = hauteur de la ligne la plus haute + 2 pixels
    // Ex : rép : 10, Filtre : 8

    // Sauf si : legende.dim trop élevé, et seuillage à sy_graphique / (2 * n)

    if((lx_max < 0) || (h1 <= 8))
    {
      DBG(msg("raxes: espace trop petit pour la légende."));
      // Légende non affichée
      return;
    }

    int yoffset = 0;
    if(!config.titre.empty())
      yoffset = 35;

    // Nord ouest
    Rect r = Rect{marge_gauche + 2, 2+yoffset, largeur_totale, hauteur_totale};

    switch(config.legende.position.type)
    {
      case ConfigAxes::Legende::Position::Type::NORD_EST:
        r.x = marge_gauche + sx_graphique - largeur_totale - 2;
        break;
      case ConfigAxes::Legende::Position::Type::SUD_EST:
        r.x  = marge_gauche + sx_graphique - largeur_totale - 2;
        r.y  = marge_haute + sy_graphique - hauteur_totale - 2;
        break;
      case ConfigAxes::Legende::Position::Type::SUD_OUEST:
        r.y  = marge_haute + sy_graphique - hauteur_totale - 2;
      default:
        ;
    }

    Point p0 = r.tl();

    canva.set_couleur(couleur_avant_plan);
    canva.set_remplissage(true, config.couleur_arriere_plan);
    canva.rectangle(p0.x, p0.y, p0.x + largeur_totale - 1, p0.y + hauteur_totale + 2);//- 1);

    for(auto i = 0; i < n; i++)
    {
      auto &ch = config.series[i];

      canva.set_couleur(ch.couleur);
      canva.set_epaisseur(ch.largeur_de_trait + 1);

      canva.set_align(Align::DEBUT, Align::CENTRE);

      int yc = p0.y + i * h1 + h1/2;

      // Comment faire pour que tous les labels aient la même échelle ?
      canva.set_dim_fonte(dim_texte);
      canva.texte(p0.x + 5, yc+2, ch.nom, lx_max, h1-2);


      int trait_xmin = p0.x+10+lx_max,
          trait_xmax = p0.x+largeur_totale-5;

      if(ch.remplissage)
      {
        auto s = ch.couleur.eclaircir(0.3);
        canva.set_remplissage(true, s);
        canva.rectangle(Rect{trait_xmin,  yc - 7, trait_xmax - trait_xmin + 1, 14});
        canva.set_remplissage(false);
      }
      else
      {
        canva.set_dotted(ch.trait == Trait::DOTTED);
        canva.ligne(Point{trait_xmin,yc}, Point{trait_xmax, yc});
        canva.set_dotted(false);
        if(ch.marqueur != Marqueur::AUCUN)
        {
          canva.set_remplissage(true, ch.couleur);
          canva.set_couleur(ch.couleur);
          canva.marqueur({(float) (trait_xmin + trait_xmax)/2, (float) yc}, ch.marqueur, 5);
          canva.set_remplissage(false);
        }
      }
    }
  }

};



void Axes::def_echelle(const string &x, const string &y)
{
  impl->config.axe_horizontal.echelle_logarithmique = x == "log";
  impl->config.axe_vertical.echelle_logarithmique   = y == "log";
  if(impl->config.axe_vertical.echelle_logarithmique)
  {
    impl->config.grille_majeure.sous_graduation = 10;
    impl->config.grille_mineure.lg_tiret = 10;
    impl->config.grille_mineure.lg_trou = 0;
  }
}


void Axes::active_grilles(bool actives)
{
  get_config().grille_mineure.afficher  = actives;
  get_config().grille_majeure.afficher  = actives;
}

void Axes::supprime_decorations()
{
  get_config().legende.afficher         = false;
  get_config().axe_horizontal.afficher  = false;
  get_config().axe_vertical.afficher    = false;
  active_grilles(false);
}

void Axes::set_isoview(bool isoview)
{
  impl->config.mode_isoview = isoview;
  impl->config.aspect_ratio = 1.0f;
}

void Axes::fixe_aspect_ratio(float ysx)
{
  impl->config.mode_isoview = true;
  impl->config.aspect_ratio = ysx;
}

Canva Axes::rendre(Canva canva, float xmin, float xmax, float ymin, float ymax) const
{
  return impl->rendre(canva, xmin, xmax, ymin, ymax);
}

Canva Axes::rendre(Canva canva, const Rectf &rdi) const
{
  return impl->rendre(canva, rdi);
}

void Axes::post_rendu(Canva canva_) const
{
  impl->post_rendu(canva_);
}

Axes Axes::clone()
{
  Axes res;
  *(res.impl) = *impl;
  return res;
}

void Axes::def_rdi_visible_abs(float xmin, float xmax, float ymin, float ymax)
{
  impl->config.def_rdi_visible_abs(xmin, xmax, ymin, ymax);
}

void ConfigAxes::def_rdi_visible_abs(float xmin, float xmax, float ymin, float ymax)
{
  DBG(msg("axes : def_rdi_visible_abs({}x{} -> {}x{}", xmin, ymin, xmax, ymax);)
  axe_horizontal.vmin  = xmin;
  axe_horizontal.vmax  = xmax;
  axe_vertical.vmin    = ymin;
  axe_vertical.vmax    = ymax;
}

void Axes::def_rdi_visible(float xmin, float xmax, float ymin, float ymax)
{
  impl->config.def_rdi_visible(xmin, xmax, ymin, ymax);
}

Rectf Axes::get_rdi() const
{
  return {impl->config.axe_horizontal.vmin, impl->config.axe_vertical.vmin,
    impl->config.axe_horizontal.vmax - impl->config.axe_horizontal.vmin,
    impl->config.axe_vertical.vmax - impl->config.axe_vertical.vmin
  };
}


void ConfigAxes::def_rdi_visible(float xmin, float xmax, float ymin, float ymax)
{
  DBG(msg("axes : def_rdi_visible({}x{} -> {}x{}", xmin, ymin, xmax, ymax);)

  if(isnan(xmin) || isnan(ymin)  || isnan(xmax) || isnan(ymax)
    || isinf(xmin) || isinf(ymin)  || isinf(xmax) || isinf(ymax))
  {
    msg_erreur("def rdi visible : xmin = {}, ymin = {}, xmax = {}, ymax = {}", xmin, ymin, xmax, ymax);
    return;
  }

  //tsd_assert(xmax >= xmin);

  auto dx = xmax - xmin;
  auto dy = ymax - ymin;

  if(dx <= 0)
    dx = 1;


  /*if(impl->config.mode_isoview)
  {
    auto ar = impl->config.aspect_ratio;
    if(dx * ar > dy)
    {
      auto cy = (ymin + ymax) / 2;
      auto delta = dx * ar/2;
      ymin = cy - delta;
      ymax = cy + delta;
      msg("### dx = {}, ar = {}, delta = +/- {}", dy, ar, delta);
      dy = dx;
    }
    else
    {
      auto cx = (xmin + xmax) / 2;
      auto delta = (dy/ar)/2;
      xmin = cx - delta;
      xmax = cx + delta;
      msg("### dy = {}, ar = {}, delta = +/- {}", dy, ar, delta);
      dx = dy;
    }
  }*/

  if(axe_horizontal.echelle_logarithmique)
  {
    float r = xmax / xmin;
    axe_horizontal.vmin = xmin / pow(r, 0.1);
    axe_horizontal.vmax = xmax * pow(r, 0.1);
  }
  else
  {
    axe_horizontal.vmin = xmin - dx * 0.1;
    axe_horizontal.vmax = xmax + dx * 0.1;
  }


  if(axe_vertical.echelle_logarithmique)
  {
    // PB ICI !!!
    double r = ymax / ymin;
    double r2 = pow(r, 0.1);

    axe_vertical.vmin   = ymin;
    axe_vertical.vmax   = ymax;

    if(!isinf(r2))
    {
      axe_vertical.vmin /= r2;
      axe_vertical.vmax *= r2;
    }

    //msg("**** ymin = {}, ymax = {} -> r = {}, vmax = {}.", ymin, ymax, r, axe_vertical.vmax);
  }
  else
  {
    axe_vertical.vmin   = ymin - dy * 0.1;
    axe_vertical.vmax   = ymax + dy * 0.25;
  }
}


Axes::Axes()
{
  impl = make_shared<Impl>();
  impl->axe_vertical.vertical   = true;
  impl->axe_horizontal.vertical = false;
}

ConfigAxe::operator string() const
{
  return fmt::format("[afficher={}, label={}, unite={}, vmin={}, vmax={}]", afficher, label, unite, vmin, vmax);
}

ConfigAxes::operator string() const
{
  return fmt::format(
      "[axe vertical:{}, axe horizontal:{}, titre:{}, carp:{}, iso:{}, ar:{}, nbseries={}]",
      (string) axe_vertical, (string) axe_horizontal,
      titre, couleur_arriere_plan, mode_isoview, aspect_ratio, series.size());
}

void Axes::configure(const ConfigAxes &cfg)
{
  if(impl->config != cfg)
  {
    DBG(msg("Changement config axe.");)
    DBG(msg("nv cfg rdi : x=[{}x{}], y=[{}x{}]",
        cfg.axe_horizontal.vmin, cfg.axe_horizontal.vmax,
        cfg.axe_vertical.vmin, cfg.axe_vertical.vmax);)
  }
  impl->config = cfg;
}

Dim Axes::get_dim_graphique(const Dim &dim_totale) const
{
  return impl->get_dim_graphique(dim_totale);
}


ConfigAxes &Axes::get_config()
{
  return impl->config;
}


}

