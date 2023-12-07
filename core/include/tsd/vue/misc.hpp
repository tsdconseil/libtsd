#pragma once
#include "tsd/vue.hpp"
#include "tsd/temps.hpp"

namespace tsd::vue {

struct FigureEvenementModele
{
  struct Evt
  {
    tsd::temps::DateHeure t0, t1;

    string label, titre, description;
    Couleur couleur = Couleur::Vert;
    bouléen sélectionné = non,
            hide_label  = non;
    entier user_id = -1;
  };

  struct Ligne
  {
    string nom;
    Couleur couleur_fond = Couleur::Blanc;

    vector<Evt> evts;

    void concaténation();
  };
  vector<Ligne> lignes;
  tsd::temps::DateHeure t0, t1;

  bouléen orientation_verticale = non;

  fonction<void(const Evt *evt, entier idl, entier ide)> callback_chgt_selection;
};



struct FigureEvenements: Rendable
{
  FigureEvenementModele content;

  mutable struct PlotInfos
  {
    struct Evt
    {
      Rectf rdi;
    };
    struct Ligne
    {
      Rectf rdi;
      vector<Evt> evts;
    };
    vector<Ligne> lignes;

  } infos;

  void maj(const FigureEvenementModele &content);
  void rendre(Canva canva) const;
};

}

