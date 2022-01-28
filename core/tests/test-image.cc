#include "tsd/vue/image.hpp"
#include "tsd/figure.hpp"


using namespace tsd;
using namespace tsd::vue;

int test_image()
{

  Image image(550, 350);

  //image->remplir(Couleur::Rouge);

  image.def_couleur_dessin(Couleur::Bleu);
  image.puts({0, 0}, "Hello", 1);

  image.def_couleur_dessin(Couleur::Bleu);
  image.ligne({50+5,50+5}, {100+5,50+5});
  image.ligne({50+5,50+5}, {50+5,100+5});
  image.ligne({50+5,50+5}, {100+5,100+5});
  image.ligne({50+5,50+5}, {0+5,50+5});
  image.ligne({50+5,50+5}, {50+5,0+5});
  image.ligne({50+5,50+5}, {0+5,0+5});
  image.ligne({50+5,50+5}, {0+5,100+5});
  image.ligne({50+5,50+5}, {100+5,0+5});
  image.ligne({50+5,50+5}, {0+5,0+5});

  image.ligne_aa(20, 0, 50, 100, Image::POINTILLEE);

  image.def_couleur_dessin(Couleur::Rouge);
  image.ligne({100+5,0+5}, {200+5,50+5});
  image.ligne({100+5,50+5}, {200+5,0+5});

  image.def_couleur_dessin(Couleur::Orange);
  image.ligne_aa(100+5,54+5, 200+5,104+5);
  image.ligne_aa(100+5,104+5, 200+5,54+5);

  image.def_couleur_dessin(Couleur::Noir);
  image.ligne({200+5,0+5}, {250+5,200+5});
  image.ligne({250+5,0+5}, {200+5,200+5});


  image.cercle({200,200}, 100);
  image.def_couleur_dessin(Couleur::Rouge);
  image.cercle({200,200}, 80);


  image.def_couleur_remplissage(Couleur::Violet);
  image.cercle_plein({100,200}, 50);

  image.def_couleur_remplissage(Couleur::Bleu);
  image.cercle_plein({150,200}, 50);

  image.def_couleur_dessin(Couleur::Violet);
  image.ellipse({300,200}, {400,300});
  image.def_couleur_dessin(Couleur::Cyan);
  image.ellipse({300,200}, {400,250});

  image.enregister("./build/test-log/test-image");


  auto im2 = image.redim(100, 100);
  im2.enregister("./build/test-log/test-image-r1");

  auto im3 = image.redim(1000, 1000);
  im3.enregister("./build/test-log/test-image-r2");

  //Figure f;
  //f.plot(linspace(0,1,50), linspace(0,1,50));
  //f.enregistrer("./test-log/test-image");


  //f.afficher();
  //Figure::attente_ihm();

  return 0;
}
