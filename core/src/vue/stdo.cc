#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include <cstdarg>
#include <thread>
#include <functional>
#include <filesystem>
#include <mutex>
#include <cstdio>

namespace tsd::vue {

Stdo stdo;

fonction<void ()> stdo_fin;
fonction<void (sptr<const Rendable> fig, cstring titre)> stdo_ajoute_figure;

struct Element
{
  enum Type
  {
    CHAINE,
    FIGURE
  } type;
  string chaine, fichier_image;
};

static struct StdoPrivA
{
  string dossier_sortie;
  bouléen mode_rt = oui;
  entier cnt_racine = 0;
  entier cnt_fig = 0;
  vector<Element> elems;
  std::mutex mtx;
} priv;


void Stdo::flush()
{

  //msg("flush stdo: {} éléms.", priv.elems.size());

  si(priv.elems.empty())
    retourne;

  //si(!std::filesystem::exists(dossier_sortie))
    //std::filesystem::create_directory(dossier_sortie);
  //cutils::fichiers::creation_dossier_si_inexistant(dossier_sortie);

  si(!priv.dossier_sortie.empty() && !priv.elems.empty())
  {
    soit fn = priv.dossier_sortie + "/index.html";
    string s;

    //msg("Flush stdo : {} elems...", priv.elems.size());

    pour(auto &e: priv.elems)
    {
      si(e.type == Element::CHAINE)
      {
        s += "<p>" + e.chaine + "</p>\n";
      }
      sinon si(e.type == Element::FIGURE)
      {
        s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" width=\"1000\"></p>";
        //s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" width=\"600\"></p>";
        //s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" ></p>";
      }
    }

    soit f = fopen(fn.c_str(), "wt");
    si(!f)
    {
      msg_erreur("Erreur lors de l'ouverture du fichier [{}].\n", fn);
      retourne;
    }
    fprintf(f, "%s", s.c_str());
    fclose(f);
  }

  //cutils::fichiers::enregistre_fichier_texte(fn, s);
  priv.elems.clear();
  priv.cnt_racine++;
  priv.cnt_fig = 0;

  priv.dossier_sortie = "";
}


void Stdo::def_dossier_sortie(cstring chemin)
{
  flush();

  priv.mode_rt        = non;
  priv.dossier_sortie = chemin;

  si(!std::filesystem::exists(chemin))
    std::filesystem::create_directory(chemin);
  si(!std::filesystem::exists(chemin + "/img"))
    std::filesystem::create_directory(chemin + "/img");

  //msg("Stdo: dossier de sortie = {}", chemin);
}

void Stdo::affiche(sptr<const Rendable> s, cstring titre, const Dim &dim)
{
  si(priv.mode_rt)
  {
    si(stdo_ajoute_figure)
      stdo_ajoute_figure(s, titre);
  }
  sinon
  {
    Element elmt;
    elmt.type = Element::FIGURE;
    soit nom_fichier = sformat("img-{}-{}.png", priv.cnt_racine, priv.cnt_fig++);
    elmt.fichier_image = nom_fichier;
    elmt.chaine        = titre;
    s->enregistrer(priv.dossier_sortie + "/img/" + nom_fichier, dim);
    priv.elems.push_back(elmt);
  }
}

void Stdo::printf(cstring s)
{
  priv.mtx.lock();
  Element elmt;
  elmt.type = Element::CHAINE;
  elmt.chaine = s;
  priv.elems.push_back(elmt);
  //::printf("priv : %d éléms.\n", (entier) priv.elems.size());
  priv.mtx.unlock();
}




void Stdo::fin()
{
  si(priv.mode_rt)
  {
    //msg("tsd::stdo : attente fin ihm.");
    si(stdo_fin)
      stdo_fin();
  }
}

}


