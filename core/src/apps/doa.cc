#include "tsd/apps/doa.hpp"
#include "tsd/stats.hpp"

using namespace tsd::stats;




namespace tsd::apps::doa {


  // Nb pos    = nb senseurs
  // Unité pos = longueur d'onde
  // Nb angles = nb sources
  // Angle = 0 <=> tangent, π/2 <=> perpendiculaire
  // Renvoie une matrice Nr par Ns
  Tabcf steervec_1d(const Vecf &pos, const Vecf &angle)
  {
    soit Ns = angle.rows(), Nr = pos.rows();
    Tabcf A(Nr, Ns);
    pour(auto i = 0; i < Nr; i++)
      pour(auto j = 0; j < Ns; j++)
        A(i,j) = std::polar(1.0f, pos(i) * cos(angle(j)));
    retourne A;
  }


  Tabcf sensorcov_1d(const Vecf &pos, const Vecf &angle, float SNR_db)
  {
    // x = A s
    // R = E[ x x*] = A E[s s*] A* = k A A*
    // Avec colonnes de A = steering vectors
    soit Nr = pos.rows();
    Tabcf A = steervec_1d(pos, angle);
    retourne A * A.transpose().conjugate() + db2pow(SNR_db) * Tabf::eye(Nr);
  }


  /** @brief DOA linéaire 1d, avec micro / antennes équidistantes
   *  @param Ns Nombre de sources (ou -1 si détermination automatique)
   *  @param d  Distance entre deux antennes (en unité de longueur d'onde de la fréquence porteuse)
   *  @returns Vecteur des angles d'arrivée des signaux sources */
  Vecf musicdoa_1d(const Tabcf &R, float d, entier Ns)
  {
    SubSpaceSpectrumConfig config;
    config.debug_actif = oui;
    config.balayage = [&](entier i, entier n, entier m) -> tuple<Vecf, Veccf>
    {
      // Calcul du "vecteur de steering"
      // (ici simple exponentielle complexe, de fréquence dans l'intervalle [-0.5,0.5[,
      // c'est à dire adaptée à la recherche de signaux périodiques)
      Vecf phi(1);
      // Balayage entre -π / 2 et π / 2
      phi(0) = π * (1.0f*i) /n-1;

      // si phi = 0, fréquence = d

      retourne {phi, sigexp(cos(phi(0)) * d, m)};
    };
    soit res = subspace_spectrum(R, config);

    Vecf phi = res.var.col(0);//, res.spectrum};

    // Recherche des Ns signaux les plus forts

    Vecf lst(res.Ns);
    Vecf sp = res.spectrum;
    pour(auto i = 0; i < res.Ns; i++)
    {
      soit index = sp.index_max();
      sp(index) = 0;
      lst(i) = res.var(index, 0);
    }

    //inline void sort(){
    //   std::sort(derived().data(), derived().data()+size());
    //}

    retourne lst;
  }



}
