#include "tsd/apps/doa.hpp"
#include "tsd/stats.hpp"

using namespace tsd::stats;



namespace tsd::apps::doa {



  // Nb pos    = nb senseurs
  // Unité pos = longueur d'onde
  // Nb angles = nb sources
  // Angle = 0 <=> tangent, π/2 <=> perpendiculaire
  // Renvoie une matrice Nr par Ns
  MatrixXcf steervec_1d(const ArrayXf &pos, const ArrayXf &angle)
  {
    int Ns = angle.rows(), Nr = pos.rows();
    MatrixXcf A(Nr, Ns);
    for(auto i = 0; i < Nr; i++)
      for(auto j = 0; j < Ns; j++)
        A(i,j) = std::polar(1.0f, pos(i) * std::cos(angle(j)));
    return A;
  }


  MatrixXcf sensorcov_1d(const ArrayXf &pos, const ArrayXf &angle, float SNR_db)
  {
    // x = A s
    // R = E[ x x*] = A E[s s*] A* = k A A*
    // Avec colonnes de A = steering vectors
    int /*Ns = angle.rows(),*/ Nr = pos.rows();
    MatrixXcf A = steervec_1d(pos, angle);
    return A * A.adjoint() + db2pow(SNR_db) * MatrixXf::Identity(Nr, Nr);
  }


  /** @brief DOA linéaire 1d, avec micro / antennes équidistantes
   *  @param Ns Nombre de sources (ou -1 si détermination automatique)
   *  @param d  Distance entre deux antennes (en unité de longueur d'onde de la fréquence porteuse)
   *  @returns Vecteur des angles d'arrivée des signaux sources */
  ArrayXf musicdoa_1d(const MatrixXcf &R, float d, int Ns)
  {
    SubSpaceSpectrumConfig config;
    config.debug_actif = true;
    config.balayage = [&](int i, int n, int m) -> std::tuple<ArrayXf, ArrayXcf>
    {
      // Calcul du "vecteur de steering"
      // (ici simple exponentielle complexe, de fréquence dans l'intervalle [-0.5,0.5[,
      // c'est à dire adaptée à la recherche de signaux périodiques)
      ArrayXf phi(1);
      // Balayage entre -π / 2 et π / 2
      phi(0) = π * (1.0f*i) /n-1;

      // si phi = 0, fréquence = d

      return {phi, sigexp(cos(phi(0)) * d, m)};
    };
    auto res = subspace_spectrum(R, config);

    ArrayXf phi = res.var.col(0);//, res.spectrum};

    // Recherche des Ns signaux les plus forts

    ArrayXf lst(res.Ns);
    ArrayXf sp = res.spectrum;
    for(auto i = 0; i < res.Ns; i++)
    {
      int index;
      sp.maxCoeff(&index);
      sp(index) = 0;
      lst(i) = res.var(index, 0);
    }

    //inline void sort(){
    //   std::sort(derived().data(), derived().data()+size());
    //}

    return lst;
  }



}
