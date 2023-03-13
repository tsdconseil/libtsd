#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"

#define DBG(AA)

#include <deque>

using namespace std;


namespace tsd::vue
{

  static void set_def(auto &x, auto v)
  {
    x = (x == -1) ? v : x;
  };

  // Constantes à supprimer !!!
  // Définit les dimensions
  static const auto HAUTEUR_X_LABEL   = 25;

  static const auto FONTE_X_LABEL  = 0.6;

  void ConfigGrille::from_style(const string &style)
  {
    si(style == ".")
    {
      lg_tiret = 1;
      lg_trou  = 2;
    }
    sinon si(style == "-")
    {
      lg_tiret = 10;
      lg_trou  = 0;
    }
    sinon si(style == "--")
    {
      lg_tiret = 5;
      lg_trou  = 5;
    }
    sinon
    {
      lg_tiret = -1;
      lg_trou  = -1;
    }
  }

  string ConfigGrille::vers_style() const
  {
    si((lg_tiret == 1) && (lg_trou == 2))
      retourne ".";
    sinon si((lg_tiret == 10) && (lg_trou == 0))
      retourne "-";
    sinon si((lg_tiret == 5) && (lg_trou == 5))
      retourne "--";
    sinon
      retourne "";
  }

  bouléen ConfigAxes::est_mat_sur_clair() const
  {
    entier score_grid  = grille_majeure.couleur.lumi();
    entier score_back = couleur_arriere_plan.lumi(); // 255

    si(grille_majeure.couleur.vers_rgba() == 0) // Par défaut
      score_grid = 255 - score_back;

    retourne score_back > score_grid;
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
    pour(auto &c: codes)
      si(c.code == id)
        type = c.pos;
  }

  ConfigAxes::Legende::Position::operator string() const
  {
    si(type == NORD_EST)
      retourne "ne";
    sinon si(type == SUD_EST)
      retourne "se";
    sinon si(type == NORD_OUEST)
      retourne "no";
    sinon si(type == SUD_OUEST)
      retourne "so";
    retourne "?";
  }


  void Axe::configure(const ConfigAxe &config, entier dim)
  {
    DBG(msg("Configure axe : vmin={}, vmax={}", config.vmin, config.vmax);)

    si(config.vmin > config.vmax)
      echec("Configuration axe : vmin ({}) > vmax ({})", config.vmin, config.vmax);

    this->config = config;
    this->dim    = dim;
    decalage_zero   = (config.vmax + config.vmin) / 2;
    float diff = config.vmax - config.vmin;
    si(diff == 0)
      diff = 1;
    pixels_par_lsb = dim / diff;

    si(config.echelle_logarithmique)
    {
      float diff_db  = (log10(config.vmax) - log10(config.vmin));
      pixels_par_lsb = dim / diff_db;
      decalage_zero  = (log10(config.vmax) + log10(config.vmin)) / 2;
    }
  }


  bouléen Axe::tic_ok(float v)
  {
    soit y = valeur_vers_pixel(v);
    retourne (y >= 5) && (y < dim - 5);
  }

  bouléen Axe::tic_ok_souple(float v)
  {
    soit y = valeur_vers_pixel(v);
    retourne (y >= 0) && (y < dim);
  }


  entier Axe::valeur_vers_pixel(float v) const
  {
    si(config.echelle_logarithmique)
    {
      // vmin et vmax sont en linéaire, v aussi
      si(v == 0)
        retourne -1;
      v = log10(v);
    }
    soit x = (entier) ((v  - decalage_zero) * pixels_par_lsb);
    retourne dim / 2 + (vertical ? -x : x);
  }

  float Axe::valeur_vers_pixelf(float v) const
  {
    si(config.echelle_logarithmique)
    {
      si(v == 0)
        retourne -1;
      // vmin et vmax sont en linéaire, v aussi
      v = log10(v);
    }
    soit x = ((v  - decalage_zero) * pixels_par_lsb);
    retourne dim / 2.0f + (vertical ? -x : x);
  }

  float Axe::pixel_vers_valeur(float p) const
  {
    soit tmp = vertical ? dim/2-p : p-dim/2;
    soit res = tmp / pixels_par_lsb + decalage_zero;

    si(config.echelle_logarithmique)
      retourne pow(10.0f, res);

    retourne res;
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

  si(v == 1)
    // Dans ce cas, tic_min est déjà une puissance de 10
    res = 1;
  sinon si(v <= 2)
    res = 2;
  sinon si(v <= 2.5)
    res = 2.5;
  sinon si(v <= 5.0)
    res = 5.0;
  sinon
    res = 10.0;

  retourne res * pnorm;
}





struct Axes::Impl
{
  // TODO: regrouper les mutable
  mutable ConfigAxes config;

  // Position horizontale de l'axe vertical
  mutable entier x_grad_vert = 0, y_grad_hor = 0;

  mutable double tic_majeur = 0, tic_mineur = 0;
  mutable unsigned int tic_majeur_pixels = 0, tic_mineur_pixels = 0;
  mutable Axe axe_horizontal, axe_vertical;
  /** Dimension de la vue, en pixels */
  mutable entier sx = 1, sy = 1;
  /** Dimension allouée pour les courbes */
  mutable entier sx_graphique = 1, sy_graphique = 1;
  mutable entier marge_gauche = 0, marge_droite = 0;
  mutable entier marge_basse = 0, marge_haute = 0;
  mutable   bouléen    mat_sur_clair = non;
  mutable Couleur couleur_avant_plan = {255,255,255};
  mutable Couleur fond_grille;
  mutable Canva canva_pixels, canva_res;

  void post_rendu(Canva canva_) const
  {
    si(config.legende.afficher)
    {
      soit canva = canva_.vue(Rect{0, 0, sx, sy});
      affiche_cartouche_legende(canva);
    }
  }

  Dim get_dim_graphique(const Dim &dim_totale) const
  {
    maj_dimensions(dim_totale);
    retourne {sx_graphique, sy_graphique};
  }

  mutable int marge_valeurs = 60;

  // dim = allocation
  void maj_dimensions(const Dim &dim) const
  {
    si(get_mode_impression())
    {
      marge_valeurs = 100;
    }


    soit &cah = config.axe_horizontal,
         &cav = config.axe_vertical;
    soit &cM  = config.grille_majeure,
         &cm  = config.grille_mineure;
    sx = dim.l;
    sy = dim.h;
    marge_haute = marge_basse  = marge_gauche = marge_droite = 0;

    si(cah.afficher)
    {
      marge_basse  = cah.valeurs_hauteur_max + cah.valeurs_décalage + 2;

      si(!cah.label.empty())
        marge_basse += HAUTEUR_X_LABEL;
    }

    si(cav.afficher)
    {
      si(config.position_axe_vertical == ConfigAxes::AXE_VERTICAL_GAUCHE)
      {
        marge_gauche = marge_valeurs;
        si(!cav.label.empty())
          marge_gauche += 25;
      }
      sinon
      {
        marge_droite = marge_valeurs;
        si(!cav.label.empty())
          marge_droite += 25;
      }
    }

    soit marge = marge_gauche + marge_droite;
    soit margev = marge_basse + marge_haute;
    sx_graphique = sx > marge ? sx - marge : sx;
    sy_graphique = sy > margev ? sy - margev : sy;

    x_grad_vert = marge_gauche;

    si(config.position_axe_vertical == ConfigAxes::AXE_VERTICAL_DROITE)
      x_grad_vert = sx - marge_droite;

    y_grad_hor = sy - marge_basse;







    if(get_mode_impression())
    {
      set_def(cM.pixels_par_tic_min_x, sx / 8);
      set_def(cM.pixels_par_tic_min_y, sx / 12);
    }
    else
    {
      set_def(cM.pixels_par_tic_min_x, sx / 12);
      set_def(cM.pixels_par_tic_min_y, sy / 12);
    }

    set_def(cm.pixels_par_tic_min_y, cM.pixels_par_tic_min_y / 20);
    set_def(cm.pixels_par_tic_min_x, cM.pixels_par_tic_min_x / 20);

    cah.id = "X";
    cav.id = "Y";
  }



  Canva rendre(Canva canva_, const Rectf &rdi) const
  {
    retourne rendre(canva_, rdi.x, rdi.x + rdi.l, rdi.y + rdi.h, rdi.y);// !!
  }

  Canva rendre(Canva canva_, float xmin, float xmax, float ymin, float ymax) const
  {
    soit dim = canva_.get_allocation(); // Allocation en pixels

    maj_dimensions(dim);

    // Vue en pixels
    canva_pixels = canva_.vue({0, 0, sx, sy});


    si(isinf(xmin) || isinf(xmax) || isinf(ymin) || isinf(ymax))
    {
      msg_erreur("Canva::rendre({},{},{},{}) : borne infinie.", xmin, xmax, ymin, ymax);
      retourne canva_res;
    }


    soit &cah = config.axe_horizontal,
         &cav = config.axe_vertical;
    soit &cM  = config.grille_majeure,
         &cm  = config.grille_mineure;

    cah.vmin = xmin;
    cah.vmax = xmax;
    cav.vmax = ymax;
    cav.vmin = ymin;


    DBG(msg("rendu axes : hmin={}, hmax={}, vmin={}, vmax={}, ecart pix={}, sx={}, sy={}",
        cah.vmin, cah.vmax,
        cav.vmin, cav.vmax,
        cah.ecart_pixel,
        sx, sy);)


    si(config.mode_isoview)
    {
      // Il faut se débrouiller pour que la résolution soit identique
      // en x et en y. soit LSB / pixel identique.
      soit dh = cah.vmax - cah.vmin, dv = cav.vmax - cav.vmin;
      soit cth = (cah.vmax + cah.vmin) / 2, ctv = (cav.vmax + cav.vmin) / 2;

      // Résolution = nb lsb / pixel, en horizontal et en vertical
      soit reso_x = dh / sx_graphique, reso_y = dv / sy_graphique;
      soit ratio = (reso_x + 1e-200) / (reso_y + 1e-200);

      // Comment forcer les résolutions à être identique ?
      // On augmente l'intervalle pour la résolution la plus haute
      soit ar = config.aspect_ratio;

      // ar = ratio résolution y / x

      DBG(msg("RDI avant AR : extension_X/extension_Y = {} VS ratio d'aspet demandé = {}", ratio, ar);)
      si(ratio <  1.0 / ar)
      {
        soit δ = sx_graphique * reso_y / (2 * ar);
        cah.vmin = cth - δ;
        cah.vmax = cth + δ;
        DBG(float nv_reso_x = 2 * δ / sx_graphique;)
        DBG(msg("cah : δ=+/- {}, nv ratio = {}", δ, nv_reso_x / reso_y);)
        // 2 (delta / sxg) / reso_y = 1/ar
        // delta = sxg * resoy / (2 * ar)
      }
      sinon
      {
        soit δ = sy_graphique * reso_x / (2 * ar);
        cav.vmin = ctv - δ;
        cav.vmax = ctv + δ;
        DBG(soit nv_reso_y = 2 * δ / sy_graphique;)
        DBG(msg("cav : δ=+/- {}, ratio d'aspet final = {}", δ, reso_x / nv_reso_y);)
      }
    }

    // TODO: la même chose pour l'axe horizontal
    si(cav.echelle_logarithmique)
    {
      //DBG(msg("mode log: vmin {} -> {}, vmax {} -> {}", vmin, log10(vmin), vmax, log10(vmax)));
      soit vmin = log10(cav.vmin),
           vmax = log10(cav.vmax);

      si(floor(vmin) == floor(vmax))
      {
        // Augmente un peu l'intervalle,
        // car les tics majeurs sont au minimum par puissances entières de 10
        vmin = floor(vmin) - 0.1f;
        vmax = vmin + 1.2f;
        cav.vmin = pow(10.0f, vmin);
        cav.vmax = pow(10.0f, vmax);
      }
    }

    xmin = cah.vmin;
    xmax = cah.vmax;
    ymin = cav.vmin;
    ymax = cav.vmax;

    Rectf rdi{xmin,ymax, xmax - xmin, ymin - ymax};
    canva_res = canva_pixels.clip_alt(
        Rect{0, marge_haute, sx_graphique, sy_graphique}, xmin, xmax, ymax, ymin,
        cah.echelle_logarithmique, cav.echelle_logarithmique);


    DBG(msg("raxes: cfg axes..."));
    axe_horizontal.configure(cah, sx_graphique);
    axe_vertical.configure(cav, sy_graphique);



    DBG(msg("raxes: cpt tics..."));
    calcule_tics(axe_vertical);
    DBG(msg("Calcul tics horizontaux...");)

    si(config.mode_isoview)
    {
      axe_horizontal.ecart_tics_majeurs_valeur = axe_vertical.ecart_tics_majeurs_valeur;
      axe_horizontal.ecart_tics_majeurs_pixels = axe_vertical.ecart_tics_majeurs_pixels;
      axe_horizontal.ecart_tics_mineurs_valeur = axe_vertical.ecart_tics_mineurs_valeur;
      axe_horizontal.ecart_tics_mineurs_pixels = axe_vertical.ecart_tics_mineurs_pixels;
    }

    calcule_tics(axe_horizontal, config.mode_isoview);
    DBG(msg("ok.");)

    DBG(msg("Axes / couleurs :"));
    DBG(msg(" - ARP    = {}", config.couleur_arriere_plan));
    DBG(msg(" - grille (majeure) = {}", config.grille_majeure.couleur));

    mat_sur_clair = config.est_mat_sur_clair();
    couleur_avant_plan = mat_sur_clair ? Couleur::Noir : Couleur::Blanc;

    DBG(msg(" - mat sur clair = {}", mat_sur_clair));
    DBG(msg(" - => avant plan = {}", couleur_avant_plan));

    set_def(cM.afficher, 1);
    set_def(cm.afficher, 1);

    si(cM.lg_tiret == -1)
    {
      cM.lg_tiret = 10;
      cM.lg_trou  = 0;
    }
    si(cm.lg_tiret == -1)
    {
      cm.lg_tiret = 1;
      cm.lg_trou  = 3;
    }

    si(cM.couleur.vers_rgba() == 0)
      cM.couleur = affaiblir_couleur(couleur_avant_plan, 0.3);

    si(cm.couleur.vers_rgba() == 0)
     cm.couleur = affaiblir_couleur(cM.couleur, 0.3);


    DBG(msg(" - grille majeure = {}", cM.couleur));
    DBG(msg(" - grille mineure = {}", cm.couleur));

    DBG(msg("raxes: ok, grille..."));

    affiche_grille(canva_pixels);

    DBG(msg("raxes: ok, titre..."));

    canva_pixels.set_align(Align::CENTRE, Align::DEBUT);
    canva_pixels.set_dim_fonte(config.titre_echelle);
    canva_pixels.set_epaisseur(1);
    canva_pixels.set_couleur(couleur_avant_plan);
    canva_pixels.texte(0.5*sx_graphique, 0, config.titre, sx_graphique, sy_graphique / 4);

    DBG(msg(" raxes: fin."));
    retourne canva_res;
  }




  // Quel est le principe :
  //  - on connait vmax, vmin
  //  -
  void calcule_tics(Axe &axe, bool tics_forcés = non) const
  {
    si(axe.config.tic_manuels)
    {
      axe.tics_majeurs_pos.clear();
      axe.tics_mineurs_pos.clear();

      pour(auto p: axe.config.tics)
        axe.tics_majeurs_pos.push_back(p.first);

      retourne;
    }

    bouléen mode_log = axe.config.echelle_logarithmique;

    si(axe.config.ecart_pixel != -1)
    {
      DBG(msg("Ecart pixel = {}", axe.config.ecart_pixel);)
      axe.config.vmin = - (axe.config.ecart_valeur * axe.dim) / axe.config.ecart_pixel;
      axe.config.vmax = 0;
      // On doit reconfigurer pour pplsb, etc.
      axe.configure(axe.config, axe.dim);
    }

    double vmin = axe.config.vmin,
           vmax = axe.config.vmax;

    si(mode_log)
    {
      DBG(msg("mode log: vmin {} -> {}, vmax {} -> {}", vmin, log10(vmin), vmax, log10(vmax)));
      vmin = log10(vmin);
      vmax = log10(vmax);

      /*si(floor(vmin) == floor(vmax))
      {
        // Augmente un peu l'intervalle,
        // car les tics majeurs sont au minimum par puissances entières de 10
        vmin = floor(vmin) - 1.1f;
        vmax = vmin + 2.2f;
        axe.config.vmin = pow(10.0f, vmin);
        axe.config.vmax = pow(10.0f, vmax);
        axe.configure(axe.config, axe.dim);
      }*/
    }

    DBG(msg("Axe [{}]: dim={}, cfg.min={}, cfg.max={}, pixels_min:x={},y={}, ecart-pixel={}, mode={}",
          axe.config.id, axe.dim, vmin, vmax,
          config.grille_majeure.pixels_par_tic_min_x,
          config.grille_majeure.pixels_par_tic_min_y,
          axe.config.ecart_pixel,
          mode_log ? "log" : "lin");)

    si(vmin > vmax)
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

    float pixels_par_tic_majeur_min =
        axe.vertical ? config.grille_majeure.pixels_par_tic_min_y
            : config.grille_majeure.pixels_par_tic_min_x;

    essai:

    // Imposé par la résolution
    // eg lsb_par_pixel / pixel_par_tic_majeur_min

    // Exemple : δ = 1000, axe.dim = 100, pptm_min = 50
    // tic_majeur_min = 10 * 50 = 500 LSB



    si(tics_forcés)
    {
      tic_majeur = axe.ecart_tics_majeurs_valeur;
      tic_mineur = axe.ecart_tics_mineurs_valeur;
    }
    sinon
    {
      // Tic majeur mini, en unité LSB (ou LSB_LOG en mode logarithmique)
      // Au moins ... pixels => au moins ... en valeur utilisateur
      float tic_majeur_min = (δ / axe.dim) * pixels_par_tic_majeur_min;
      float tic_mineur_min = (δ / axe.dim) *
          (axe.vertical ? config.grille_mineure.pixels_par_tic_min_y
          : config.grille_mineure.pixels_par_tic_min_x);

      // TODO : pas une variable de classe !
      tic_majeur = calcule_tic(tic_majeur_min);



      DBG(msg("  δ={}, tcmin={} -> essai tic majeur = {}", δ, tic_majeur_min, tic_majeur);)

      // En mode logarithmique, travaille au mimimun par puissance entière de 10
      si(mode_log)
      {
        tic_majeur = ceil(tic_majeur);
      }

      si(axe.config.ecart_pixel != -1)
        tic_majeur = axe.config.ecart_pixel * (δ / axe.dim);

      si((tic_majeur == 0) || (isnan(tic_majeur)))
      {
        msg_avert("Tic majeur : {}", tic_majeur);
        retourne;
      }

      si(tic_majeur <= numeric_limits<double>::epsilon())
      {
        tic_majeur = 10 * numeric_limits<double>::epsilon();
        DBG(msg_avert("tcm < eps -> {}...", tic_majeur);)
      }


      axe.ecart_tics_majeurs_valeur = tic_majeur;
      axe.ecart_tics_majeurs_pixels = tic_majeur * (axe.dim / δ);

      tic_mineur = tic_majeur / config.grille_majeure.sous_graduation;
      axe.ecart_tics_mineurs_valeur = tic_mineur;

      axe.tics_mineurs_presents = tic_mineur >= tic_mineur_min;

      DBG(msg("  Essai tic majeur = {}, pp_par_tic_majeur_min = {} -> ecart_tic_majeurs_pix = {}...", tic_majeur, pixels_par_tic_majeur_min, axe.ecart_tics_majeurs_pixels);)
    }

    axe.tics_majeurs_pos.clear();
    axe.tics_mineurs_pos.clear();

    // Entre tic_majeur_min et δ, on cherche des valeurs
    // pas trop bêtes

    bouléen zero_present = (vmin <= 0) && (vmax >= 0);

    // Les variables tk, tic_majeur, etc. sont en logs
    soit adapt = [&](double v)
    {
      si(mode_log)
        retourne pow(10.0f, v);
      retourne v;
    };

    // Ajoute des tics mineurs au dessus de tk (tic majeur)s
    soit ajouter_tics_mineurs = [&](float tk)
    {
      si(!axe.tics_mineurs_presents)
        retourne;

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

      si(mode_log)
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
        entier k = 1;
        double tkm = log10(pow(10.0, tk) + k * R);
        tantque(tkm < tk + tic_majeur)
        {
          si(axe.tic_ok(adapt(tkm)))
            axe.tics_mineurs_pos.push_back(adapt(tkm));
          tkm = log10(pow(10.0, tk) + k * R);
          k++;
        }
      }
      sinon
      {
        double tkm = tk + tic_mineur;
        tantque(tkm < tk + tic_majeur)
        {
          si(axe.tic_ok(adapt(tkm)))
            axe.tics_mineurs_pos.push_back(adapt(tkm));
          tkm += tic_mineur;
        }
      }
    };

    si(zero_present)
    {
      // => On a forcément un tic à zéro
      double rmax = max(-vmin, vmax);
      double tk = 0;//tic_majeur;
      tantque(tk < rmax)
      {
        si(axe.tic_ok(adapt(tk)))
          axe.tics_majeurs_pos.push_back(adapt(tk));
        si((tk != 0) && axe.tic_ok(adapt(-tk)))
          axe.tics_majeurs_pos.push_back(adapt(-tk));

        ajouter_tics_mineurs(tk);
        ajouter_tics_mineurs(-tk-tic_majeur);

        tk += tic_majeur;
      }
    }
    sinon
    {
      // Dans le cas contraire, quel est le point de référence ???
      double tk;

      // vmax > vmin > 0
      si(vmin >= 0)
      {
        // => Le plus petit multiple de tic_majeur > config.min
        DBG(msg("vmin={}, tmaj={} --> ceil={}", vmin, tic_majeur, tic_majeur * ceil(vmin / tic_majeur));)
        tk = tic_majeur * ceil(vmin / tic_majeur);

        si(!axe.tic_ok(adapt(tk)))
          tk += tic_majeur;

        DBG(msg(" tk={}, tic ok = {}", tk, axe.tic_ok(adapt(tk)));)
        DBG(msg("valeur_vers_pixel : {}", axe.valeur_vers_pixel(tk));)

        tantque(axe.tic_ok_souple(adapt(tk)))
        {
          si(axe.tic_ok(adapt(tk)))
            axe.tics_majeurs_pos.push_back(adapt(tk));
          ajouter_tics_mineurs(tk);
          tk += tic_majeur;
        }
      }
      // vmin < vmax < 0
      sinon
      {
        // => Le plus grand multiple de tic_majeur < vmax
        tk = -tic_majeur * ceil(-vmax / tic_majeur);

        DBG(msg("  -> essai premier tic = {} LSB.", tk));

        tantque(axe.tic_ok_souple(adapt(tk)))
        {
          si(axe.tic_ok(adapt(tk)))
            axe.tics_majeurs_pos.push_back(adapt(tk));
          ajouter_tics_mineurs(tk);
          tk -= tic_majeur;
        }
      }
    }

    //DBG(msg("  Tics majeurs :");
    //pour(auto &p: axe.tics_majeurs_pos)
      //msg("   {:.5f} (log10 = {:.2f})", p, log10(p)););


    si((axe.config.ecart_pixel == -1) && !tics_forcés)
    {
      si((axe.tics_majeurs_pos.size() <= 1) && (pixels_par_tic_majeur_min > 1))
      {
        pixels_par_tic_majeur_min = floor(pixels_par_tic_majeur_min * 0.8f);
        si(pixels_par_tic_majeur_min >= 2)
        {
          DBG(msg("Calcul tics : aucun tic majeur trouvé, ré-essai en réduisant l'écart / pixel..."));
          goto essai;
        }
        sinon
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
    retourne mat_sur_clair ? c.eclaircir(ratio) : c.assombrir(ratio);
  }

  Couleur renforcer_couleur(const Couleur &c, float ratio = 0.5) const
  {
    retourne mat_sur_clair ? c.assombrir(ratio) : c.eclaircir(ratio);
  }

  float pixel_vers_x(float x) const
  {
    retourne axe_horizontal.pixel_vers_valeur(x);
  }

  float pixel_vers_y(float y) const
  {
    retourne axe_vertical.pixel_vers_valeur(y);
  }

  Point pos(float x, float y) const
  {
    retourne {x_vers_pixel(x), y_vers_pixel(y)};
  }

  tuple<float, float> posf(float x, float y) const
  {
    retourne {axe_horizontal.valeur_vers_pixelf(x),
              axe_vertical.valeur_vers_pixelf(y)};
  }

  entier x_vers_pixel(float v) const
  {
    retourne canva_res.v2c({v,0.0f}).x;
  }

  entier y_vers_pixel(float v) const
  {
    retourne canva_res.v2c({0.0f,v}).y;
  }


  void affiche_axe_vertical(Canva O, sptr<Canva::GroupeTextes> grp) const
  {
    soit &av  = axe_vertical;
    soit &cav = config.axe_vertical;

    si(!cav.afficher)
      retourne;

    soit x0 = x_grad_vert;

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    O.ligne(x0, marge_haute, x0, marge_haute+sy_graphique-1);
    O.set_dim_fonte(cav.valeurs_echelle);

    soit [expo, nbchiffres] = tsd::vue::unites::calc_expo_nb_chiffres_commun(
        av.tics_majeurs_pos, cav.unité);

    si(av.tics_majeurs_pos.empty())
      msg_avert("Axes : aucun tic vertical !");

    pour(auto t: av.tics_majeurs_pos)
    {
      string s;
      si(av.config.echelle_logarithmique)
        s = fmt::format("{:.1e}", t);
      sinon
        s = tsd::vue::unites::valeur_vers_chaine(t, cav.unité, expo, nbchiffres);
      affiche_tic_majeur(O, t, 5, x0, s, grp);
    }
    soit lg = 1;
    pour(auto t: av.tics_mineurs_pos)
      affiche_tic_mineur(O, t, lg, x0);

    O.set_dim_fonte(FONTE_X_LABEL);
    si((!cav.label.empty()) && (sy > 20) && (sx > 100))
    {
      soit label_y = config.axe_vertical.label;
      si(label_y.size() > 1)
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
    si(config.grille_majeure.afficher)
    {
      si(axe_vertical.tics_mineurs_presents && config.grille_mineure.afficher)
      {
        pour(auto &t: axe_vertical.tics_mineurs_pos)
        {
          entier y = canva_res.v2c(Pointf{0, t}).y;
          ligne_stylee_horizontale(canva, {marge_gauche, y},
                                  {sx-1-marge_droite,y},
                                   config.grille_mineure);
        }
      }

      DBG(msg(" raxes/affiche_grille : lignes v"));
      si(axe_horizontal.tics_mineurs_presents && config.grille_mineure.afficher)
      {
        pour(auto &t: axe_horizontal.tics_mineurs_pos)
        {
          entier x = canva_res.v2c(Pointf{t,0}).x;
          ligne_stylee_verticale(canva,
                                 {marge_gauche+x, marge_haute},
                                 {marge_gauche+x,sy-1-marge_basse},
                                 config.grille_mineure);
        }
      }

      DBG(msg(" raxes/affiche_grille : lignes V"));
      pour(auto t: axe_horizontal.tics_majeurs_pos)
      {
        entier x = canva_res.v2c(Pointf{t,0}).x;
        ligne_stylee_verticale(canva,
                               {marge_gauche+x, marge_haute},
                               {marge_gauche+x, sy-1-marge_basse},
                               config.grille_majeure);
      }

      DBG(msg(" raxes/affiche_grille : lignes H"));
      pour(auto &t: axe_vertical.tics_majeurs_pos)
      {
        entier y = canva_res.v2c(Pointf{0.0f,(float)t}).y;
        ligne_stylee_horizontale(canva, {marge_gauche, y},
                                 {sx-1-marge_droite,y},
                                 config.grille_majeure);
      }
    }

    DBG(msg(" raxes/affiche_grille : axes"));
    auto grp = canva.groupe_textes();
    affiche_axe_vertical(canva, grp);
    affiche_axe_horizontal(canva, grp);
    DBG(msg(" raxes/affiche_grille : fait"));
  }

  void affiche_axe_horizontal(Canva O, sptr<Canva::GroupeTextes> grp) const
  {
    soit &ah  = axe_horizontal;
    soit &cah = config.axe_horizontal;

    si(!cah.afficher)
      retourne;

    soit y0 = y_grad_hor;

    DBG(msg(" affichage axe horizontal : sy={}, y_grad_hor={}", sy, y_grad_hor);)

    O.set_couleur(couleur_avant_plan);
    O.ligne(Point{marge_gauche, y0}, Point{marge_gauche+sx_graphique-1, y0});

    soit [expo, nbchiffres] = tsd::vue::unites::calc_expo_nb_chiffres_commun(
        ah.tics_majeurs_pos, cah.unité);

    DBG(msg(" axe horizontal: nb chiffres = {}, expo = {}", nbchiffres, expo));

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    O.set_align(Align::CENTRE, Align::DEBUT);

    O.set_dim_fonte(cah.valeurs_echelle);

    // Pb : à ce moment ci, on ne connait pas encore la dimension du bitmap,
    // et donc on ne connait pas la dimension finale des textes...
    // Conséquence : tous les tics n'auront pas la même fonte !

    // Solutions possibles :
    //   - (1) lier les textes ici, pour que au moment du rendu du canva, la plus petite fonte soit appliquée ?
    //   - (2) ?

    soit idx = 0;

    pour(auto t: ah.tics_majeurs_pos)
    {
      soit ti = x_vers_pixel(t);
      O.ligne(Point{ti+marge_gauche, y0}, Point{ti+marge_gauche, y0 - 5});
      soit s = tsd::vue::unites::valeur_vers_chaine(t, cah.unité, expo, nbchiffres);

      si(cah.tic_manuels && (idx < (entier) cah.tics.size()))
        s = cah.tics[idx].second;

      // * 0.8 pour avoir un miminum d'espace
      soit dxmax = ah.ecart_tics_majeurs_pixels * 0.8;

      si(marge_gauche+ti-dxmax/2 >= 0)
        O.texte(marge_gauche+ti, y0 + cah.valeurs_décalage, s,
                dxmax, cah.valeurs_hauteur_max,
                grp);

      idx++;
    }

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    pour(auto t: axe_horizontal.tics_mineurs_pos)
    {
      soit x = x_vers_pixel(t);
      si(y0 < sy/2)
        O.ligne(Point{x, y0+1}, Point{x, y0});
      sinon
        O.ligne(Point{x, y0-1}, Point{x, y0});
    }

    si((!cah.label.empty()) && (sy > 20) && (sx > 100))
    {
      O.set_dim_fonte(FONTE_X_LABEL);
      O.texte(sx/2, sy-HAUTEUR_X_LABEL-1, cah.label, 0.8*sx_graphique, HAUTEUR_X_LABEL);
    }
  }

  void affiche_tic_majeur(Canva &O, float v, entier lg, entier x0, const string &val, sptr<Canva::GroupeTextes> grp) const
  {
    soit y = y_vers_pixel(v);

    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    si(x0 < sx/2)
      O.ligne(Point{x0+lg, y}, Point{x0, y});
    sinon
      O.ligne(Point{x0-lg, y}, Point{x0, y});

    //soit s = valeur_vers_chaine(v, config.axe_vertical.unite);

    entier pos = x0 + 6;
    Align align = Align::DEBUT;
    si(x0 < sx/2)
    {
      pos = 20;
      align = Align::FIN;
    }

    O.set_align(align, Align::CENTRE);

    soit dxmax = marge_valeurs - 6,
         dymax = config.axe_vertical.valeurs_hauteur_max;

    si(y - dymax / 2 >= 0)
      O.texte(pos, y, val, dxmax, dymax, grp);
  }


  void affiche_tic_mineur(Canva &O, float v, entier lg = 2, entier x0 = 20) const
  {
    soit y = y_vers_pixel(v);
    O.set_couleur(couleur_avant_plan);
    O.set_epaisseur(1);
    si(x0 < sx/2)
      O.ligne(Point{x0+lg, y}, Point{x0, y});
    sinon
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
    canva.set_dotted(non);
  }

  void ligne_stylee_verticale(Canva canva,
      const Point p0, const Point p1,
      const ConfigGrille &style) const
  {
    canva.set_couleur(style.couleur);
    canva.set_epaisseur(style.epaisseur);
    canva.set_dotted(style.lg_trou != 0);
    canva.ligne(p0.x, p0.y, p0.x, p1.y);
    canva.set_dotted(non);
  }


  void affiche_cartouche_legende(Canva &canva) const
  {
    DBG(msg("raxes: legende..."));
    si((sx_graphique < 100) || (sy_graphique < 100))
      // Pas assez de place
      retourne;

    soit a_description = non;
    pour(auto c: config.series)
    {
      si(c.nom.size() > 0)
      {
        a_description = oui;
        break;
      }
    }

    si(!a_description)
    {
      DBG(msg("raxes: pas série configurée --> pas de légende."));
      retourne;
    }

    entier n = config.series.size();

    Image tmp;
    soit dim_texte = config.legende.dim;

    DBG(msg("Calcul dim max légende : échelle = {}", dim_texte));

    // Calcul dimension maximale des différents textes
    soit lx_max = 0, ly_max = 0;
    pour(auto i = 0; i < n; i++)
    {
      soit dim = tmp.texte_dim(config.series[i].nom, dim_texte);
      DBG(msg("dim[{}] : {}", i, dim);)
      lx_max = max(dim.l, lx_max);
      ly_max = max(dim.h, ly_max);
    }

    soit afficher_traits = non;
    pour(auto i = 0; i < n; i++)
    {
      soit &ch = config.series[i];
      si(ch.remplissage || (ch.marqueur != Marqueur::AUCUN))
      {
        afficher_traits = oui;
        break;
      }
    }

    soit hauteur_max     = sy_graphique / 2,
         hauteur_totale  = (ly_max + 2) * n,
         largeur_max     = sx_graphique / 3,
         largeur_totale  = lx_max + (afficher_traits ? 50 : 0);

    DBG(msg("l = {}, lmax = {}, h = {}, hmax = {}", largeur_totale, largeur_max, hauteur_totale, hauteur_max));

    largeur_totale = min(largeur_max, largeur_totale);

    si(hauteur_totale > hauteur_max)
    {
      // Diminue l'échelle de manière à avoir la même échelle finale pour tous les éléments
      dim_texte *= (1.0f * hauteur_max) / hauteur_totale;
      hauteur_totale = hauteur_max;
    }

    DBG(msg("Après seuillage : largeur totale = {}, hauteur totale = {}", largeur_totale, hauteur_totale));



    lx_max = largeur_totale;// - 50;

    si(afficher_traits)
      lx_max -= 50;

    soit h1 = (entier) ceil((1.0f * hauteur_totale) / n);
    hauteur_totale = h1 * n;

    // -> h1 = 30

    DBG(msg("Après arrondi : largeur totale = {}, hauteur totale = {}", largeur_totale, hauteur_totale));
    DBG(msg("h1 = {}", h1));

    // Normalement, h1 = hauteur de la ligne la plus haute + 2 pixels
    // Ex : rép : 10, Filtre : 8

    // Sauf si : legende.dim trop élevé, et seuillage à sy_graphique / (2 * n)

    si((lx_max < 0) || (h1 <= 8))
    {
      DBG(msg("raxes: espace trop petit pour la légende."));
      // Légende non affichée
      retourne;
    }

    soit yoffset = 0;
    si(!config.titre.empty())
      yoffset = 35;

    soit dx = 8, dy = 8;

    // Nord ouest
    Rect r{marge_gauche + dx, dy + yoffset,
           largeur_totale, hauteur_totale};

    switch(config.legende.position.type)
    {
      case ConfigAxes::Legende::Position::Type::NORD_EST:
        r.x = marge_gauche + sx_graphique - largeur_totale - dx;
        break;
      case ConfigAxes::Legende::Position::Type::SUD_EST:
        r.x  = marge_gauche + sx_graphique - largeur_totale - dx;
        r.y  = marge_haute + sy_graphique - hauteur_totale - dy;
        break;
      case ConfigAxes::Legende::Position::Type::SUD_OUEST:
        r.y  = marge_haute + sy_graphique - hauteur_totale - dy;
      default:
        ;
    }

    //msg("rendu axes : pos légende = {}", (entier) config.legende.position.type);

    soit p0 = r.tl();

    canva.set_couleur(afficher_traits ? couleur_avant_plan : config.couleur_arriere_plan);
    canva.set_remplissage(oui, config.couleur_arriere_plan);

    canva.rectangle(p0.x, p0.y,
                    p0.x + largeur_totale - 1,
                    p0.y + hauteur_totale + 2);

    sptr<Canva::GroupeTextes> grp = canva.groupe_textes();




    pour(auto i = 0; i < n; i++)
    {
      soit &ch = config.series[i];

      canva.set_couleur(ch.couleur);
      canva.set_epaisseur(ch.largeur_de_trait + 1);

      canva.set_align(Align::DEBUT, Align::CENTRE);

      entier yc = p0.y + i * h1 + h1/2;

      // Comment faire pour que tous les labels aient la même échelle ?
      canva.set_dim_fonte(dim_texte);
      canva.texte(p0.x + 5, yc+2, ch.nom, lx_max, h1-6, grp);

      entier trait_xmin = p0.x+10+lx_max,
             trait_xmax = p0.x+largeur_totale-5;

      si(afficher_traits)
      {
        si(ch.remplissage)
        {
          soit s = ch.couleur.eclaircir(0.3);
          canva.set_remplissage(oui, s);
          canva.rectangle(Rect{trait_xmin,  yc - 7, trait_xmax - trait_xmin + 1, 14});
          canva.set_remplissage(non);
        }
        sinon
        {
          si(ch.trait != Trait::AUCUN)
          {
            canva.set_dotted(ch.trait == Trait::DOTTED);
            canva.ligne(Point{trait_xmin,yc}, Point{trait_xmax, yc});
            canva.set_dotted(non);
          }
          si(ch.marqueur != Marqueur::AUCUN)
          {
            canva.set_remplissage(oui, ch.couleur);
            canva.set_couleur(ch.couleur);
            canva.marqueur({(trait_xmin + trait_xmax)/2, yc}, ch.marqueur, 5);
            canva.set_remplissage(non);
          }
        }
      }
    }
  }

};



void Axes::def_echelle(const string &x, const string &y)
{
  impl->config.axe_horizontal.echelle_logarithmique = x == "log";
  impl->config.axe_vertical.echelle_logarithmique   = y == "log";
  si(impl->config.axe_vertical.echelle_logarithmique)
  {
    impl->config.grille_majeure.sous_graduation = 10;
    impl->config.grille_mineure.lg_tiret        = 10;
    impl->config.grille_mineure.lg_trou         = 0;
  }
}


void Axes::active_grilles(bouléen actives)
{
  get_config().grille_mineure.afficher  = actives;
  get_config().grille_majeure.afficher  = actives;
}

void Axes::supprime_decorations()
{
  get_config().legende.afficher         = non;
  get_config().axe_horizontal.afficher  = non;
  get_config().axe_vertical.afficher    = non;
  active_grilles(non);
}

void Axes::set_isoview(bouléen isoview)
{
  impl->config.mode_isoview = isoview;
  impl->config.aspect_ratio = 1.0f;
}

void Axes::fixe_aspect_ratio(float ysx)
{
  impl->config.mode_isoview = oui;
  impl->config.aspect_ratio = ysx;
}

Canva Axes::rendre(Canva canva, float xmin, float xmax, float ymin, float ymax) const
{
  retourne impl->rendre(canva, xmin, xmax, ymin, ymax);
}

Canva Axes::rendre(Canva canva, const Rectf &rdi) const
{
  retourne impl->rendre(canva, rdi);
}

void Axes::post_rendu(Canva canva_) const
{
  impl->post_rendu(canva_);
}

Axes Axes::clone()
{
  Axes res;
  *(res.impl) = *impl;
  retourne res;
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
  retourne {impl->config.axe_horizontal.vmin, impl->config.axe_vertical.vmin,
    impl->config.axe_horizontal.vmax - impl->config.axe_horizontal.vmin,
    impl->config.axe_vertical.vmax - impl->config.axe_vertical.vmin
  };
}


void ConfigAxes::def_rdi_visible(float xmin, float xmax, float ymin, float ymax)
{
  DBG(msg("axes : def_rdi_visible({}x{} -> {}x{}", xmin, ymin, xmax, ymax);)

  si(isnan(xmin) || isnan(ymin)  || isnan(xmax) || isnan(ymax)
    || isinf(xmin) || isinf(ymin)  || isinf(xmax) || isinf(ymax))
  {
    msg_erreur("def rdi visible : xmin = {}, ymin = {}, xmax = {}, ymax = {}", xmin, ymin, xmax, ymax);
    retourne;
  }

  //tsd_assert(xmax >= xmin);

  soit dx = xmax - xmin;
  soit dy = ymax - ymin;

  si(dx <= 0)
    dx = 1;

  si(axe_horizontal.echelle_logarithmique)
  {
    float r = xmax / xmin;
    axe_horizontal.vmin = xmin / pow(r, 0.1);
    axe_horizontal.vmax = xmax * pow(r, 0.1);
  }
  sinon
  {
    axe_horizontal.vmin = xmin - dx * 0.1;
    axe_horizontal.vmax = xmax + dx * 0.1;
  }


  si(axe_vertical.echelle_logarithmique)
  {
    // PB ICI !!!
    double r = ymax / ymin;
    double r2 = pow(r, 0.1);

    axe_vertical.vmin   = ymin;
    axe_vertical.vmax   = ymax;

    si(!isinf(r2))
    {
      axe_vertical.vmin /= r2;
      axe_vertical.vmax *= r2;
    }

    //msg("**** ymin = {}, ymax = {} -> r = {}, vmax = {}.", ymin, ymax, r, axe_vertical.vmax);
  }
  sinon
  {
    axe_vertical.vmin   = ymin - dy * 0.1;
    axe_vertical.vmax   = ymax + dy * 0.25;
  }
}


Axes::Axes()
{
  impl = make_shared<Impl>();
  impl->axe_vertical.vertical   = oui;
  impl->axe_horizontal.vertical = non;
}

ConfigAxe::operator string() const
{
  retourne fmt::format("[afficher={}, label={}, unite={}, vmin={}, vmax={}]", afficher, label, unité, vmin, vmax);
}

ConfigAxes::operator string() const
{
  retourne fmt::format(
      "[axe vertical:{}, axe horizontal:{}, titre:{}, carp:{}, iso:{}, ar:{}, nbseries={}]",
      (string) axe_vertical, (string) axe_horizontal,
      titre, couleur_arriere_plan, mode_isoview, aspect_ratio, series.size());
}

void Axes::configure(const ConfigAxes &cfg)
{
  si(impl->config != cfg)
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
  retourne impl->get_dim_graphique(dim_totale);
}


ConfigAxes &Axes::get_config()
{
  retourne impl->config;
}


}

