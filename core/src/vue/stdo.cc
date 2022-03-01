#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include <cstdarg>
#include <thread>
#include <functional>
#include <filesystem>
#include <mutex>
#include <cstdio>

namespace tsd::vue {

Stdo stdo;

std::function<void ()> stdo_fin;
std::function<void (sptr<const Rendable> fig, const std::string &titre)> stdo_ajoute_figure;

struct Element
{
  enum Type
  {
    CHAINE,
    FIGURE
  } type;
  std::string chaine, fichier_image;
};

static struct StdoPrivA
{
  std::string dossier_sortie;
  bool mode_rt = true;
  int cnt_racine = 0;
  int cnt_fig = 0;
  std::vector<Element> elems;
  std::mutex mtx;
} priv;


void Stdo::flush()
{

  if(priv.elems.empty())
    return;

  //if(!std::filesystem::exists(dossier_sortie))
    //std::filesystem::create_directory(dossier_sortie);
  //cutils::fichiers::creation_dossier_si_inexistant(dossier_sortie);

  if(!priv.dossier_sortie.empty() && !priv.elems.empty())
  {
    auto fn = priv.dossier_sortie + "/index.html";
    std::string s;

    //msg("Flush stdo : {} elems...", priv.elems.size());

    for(auto &e: priv.elems)
    {
      if(e.type == Element::CHAINE)
      {
        s += "<p>" + e.chaine + "</p>\n";
      }
      else if(e.type == Element::FIGURE)
      {
        s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" width=\"1000\"></p>";
        //s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" width=\"600\"></p>";
        //s += "<p><h3>" + e.chaine + "</h3><img src=\"img/" + e.fichier_image + "\" ></p>";
      }
    }

    auto f = fopen(fn.c_str(), "wt");
    if(f == nullptr)
    {
      msg_erreur("Erreur lors de l'ouverture du fichier [{}].\n", fn);
      return;
    }
    fprintf(f, "%s", s.c_str());
    fclose(f);
  }

  //cutils::fichiers::enregistre_fichier_texte(fn, s);
  priv.elems.clear();
  priv.cnt_racine++;
  priv.cnt_fig = 0;
}


void Stdo::def_dossier_sortie(const std::string &chemin)
{
  flush();

  priv.mode_rt        = false;
  priv.dossier_sortie = chemin;

  if(!std::filesystem::exists(chemin))
    std::filesystem::create_directory(chemin);
  if(!std::filesystem::exists(chemin + "/img"))
    std::filesystem::create_directory(chemin + "/img");

  //msg("Stdo: dossier de sortie = {}", chemin);
}

void Stdo::affiche(sptr<const Rendable> s, const std::string &titre, const Dim &dim)
{
  if(priv.mode_rt)
  {
    if(stdo_ajoute_figure)
      stdo_ajoute_figure(s, titre);
  }
  else
  {
    Element elmt;
    elmt.type = Element::FIGURE;
    auto nom_fichier = fmt::format("img-{}-{}.png", priv.cnt_racine, priv.cnt_fig++);
    elmt.fichier_image = nom_fichier;
    elmt.chaine        = titre;
    s->enregistrer(priv.dossier_sortie + "/img/" + nom_fichier, dim);
    priv.elems.push_back(elmt);
  }
}

void Stdo::printf(const std::string &s)
{
  priv.mtx.lock();
  Element elmt;
  elmt.type = Element::CHAINE;
  elmt.chaine = s;
  priv.elems.push_back(elmt);
  //::printf("priv : %d éléms.\n", (int) priv.elems.size());
  priv.mtx.unlock();
}




void Stdo::fin()
{
  if(priv.mode_rt)
  {
    //msg("tsd::stdo : attente fin ihm.");
    if(stdo_fin)
      stdo_fin();
  }
}

}


