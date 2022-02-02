#include "tsd/fourier.hpp"
#include "tsd/figure.hpp"
#include "tsd/filtrage.hpp"

using namespace std;
using namespace tsd::vue;

#define VERB(AA)

namespace tsd::fourier {



// Interpolation quadratique,
// d'après https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
static float qint_loc(float ym1, float y0, float yp1)
{
  return (yp1 - ym1) / (2 * (2*y0 - yp1 - ym1));
}

// δ = position fractionnaire du max.
template<typename T>
static T qint_val(T ym1, T y0, T yp1, float δ)
{
  return y0 - (ym1 - yp1) * δ * ((T) 0.25);
}


// Ligne à retard spéciale, capable de retrouver les K derniers échantillons à tout moment
template<typename T>
struct LigneARetardExt
{
  int K = 0;
  Vecteur<T> mem;

  void configure(int K)
  {
    this->K = K;
    mem.setZero(K);
  }

  void step(const Vecteur<T> &x)
  {
    int n = x.rows();
    if(n == K)
      mem = x;
    else if(n > K)
      mem = x.tail(K);
    else
    {
      mem.head(K-n) = mem.tail(K-n);
      mem.tail(n) = x;
    }
  }
  const Vecteur<T> &get_last_K() const
  {
    return mem;
  }
  Vecteur<T> derniers(int n) const
  {
    return mem.tail(n);
  }
  Vecteur<T> get(int delais, int n) const
  {
    if((delais >= K) || (n >= delais))
    {
      echec("Ligne à retard ext : dépassement, K={}, delais={}, n={}", K, delais, n);
    }
    return mem.segment(K - delais - 1, n);
  }
};


/** Détecteur de motif fixe */
struct DetecteurImpl: Detecteur
{
  /** Ligne à retard 1 */
  LigneARetardExt<cfloat> lar;

  /** Ligne à retard 2 */
  sptr<FiltreGen<float, float>> lar1;

  /** Moniteur utilisation CPU */
  vector<MoniteurCpu> mon;

  /** Corrélateur (par filtrage OLA ou RIF) */
  sptr<FiltreGen<cfloat, cfloat>> ola;

  /** Filtre à moyenne glissante pour l'énergie du signal sur les M derniers échantillons */
  sptr<FiltreGen<float>> filtre_energie;

  int itr = 0;

  /** Ne = paquets d'entrée
   *  N  = dimension FFT (Ne + Nz), valable seulement en mode OLA
   *  M  = dimension motif
   */
  int Ne = 0, N = 1, M = 0;

  /** TFD du motif, pré-calculée (seulement en mode OLA). */
  ArrayXcf T_motif;

  /** Norme 2 du motif */
  float norme_motif;

  /** Motif, normalisé */
  ArrayXcf motif;

  bool pic_final_a_traiter = false;
  Detection pic_final;
  int dernier_n;

  MoniteursStats moniteurs()
  {
    return {{mon[0].stats(), mon[1].stats(), mon[2].stats()}};
  }

  DetecteurImpl(const DetecteurConfig &config): mon(3)
  {
    mon[0].nom() = "fft-corr/energie";
    mon[1].nom() = "fft-corr/ola";
    mon[2].nom() = "fft-corr/norm";
    this->config = config;
    configure_impl(config);
  }

  int configure_impl(const DetecteurConfig &config)
  {
    pic_final_a_traiter = false;
    itr = 0;

    M = config.motif.rows();

    filtre_energie = tsd::filtrage::filtre_mg<float,double>(M);

    Ne = config.Ne;

    if(config.Ne == 0)
    {
      float C;
      int Nf, Nz;
      tsd::fourier::ola_complexite_optimise(M, C, Nf, Nz, Ne);
      msg("FFT corrélateur : calcul auto Ne optimal : M={} --> Ne={},Nz={},Nf=Ne+Nz={},C={} FLOPS/ech", M, Ne, Nz, Nf, C);
    }

    norme_motif = sqrt(config.motif.abs2().sum());

    // Normalisation énergie du motif
    motif = config.motif / norme_motif;


    if(config.mode == DetecteurConfig::MODE_OLA)
    {
      FiltreFFTConfig ola_config;
      ola_config.nb_zeros_min             = M - 1;
      ola_config.dim_blocs_temporel       = Ne;
      ola_config.traitement_freq =
          [&](ArrayXcf &X)
          {
            tsd_assert(X.rows() == T_motif.rows());
            tsd_assert(X.rows() == N);
            X *= T_motif.conjugate();
          };

      tie(ola, N) = filtre_fft(ola_config);

      ArrayXcf tmp;
      tmp.setZero(N);
      tsd_assert(M * 2 <= N);

      tmp.head(M) = motif;

      delais_corr = Ne;

      T_motif = fft(tmp);
      //lar1 = tsd::filtrage::ligne_a_retard<float>(delais_corr - M - 1);
      lar1 = tsd::filtrage::ligne_a_retard<float>(delais_corr - M + 1);
    }
    else
    {
      msg("Détecteur : mode filtre RIF.");
      // TODO: détecter si le motif est réel -> filtre réel

      if(config.motif.imag().abs().sum() / config.motif.abs().sum() < 1e-7)
      {
        msg("  Motif réel détecté.");
        ola = tsd::filtrage::filtre_rif<float, cfloat>(config.motif.reverse().real());
      }
      else
      {
        msg("  Motif complexe détecté.");
        ola = tsd::filtrage::filtre_rif<cfloat, cfloat>(config.motif.reverse().conjugate());
      }

      delais_corr = M - 1;
    }

    lar.configure(delais_corr+1);

    msg("Corrélateur FFT: Ne (dim blocs) = {}, M (dim motif) = {}, N (dim fft) = {}.", Ne, M, N);
    msg("  norme motif = {}", norme_motif);
    return 0;
  }


  cfloat lc = 0.0f, lc0 = 0.0f;
  float alc = 0, alc0 = 0;

  int delais_corr = 0;

  void step(IArrayXcf x, ArrayXf &y)
  {
    auto &config = Configurable<DetecteurConfig>::lis_config();

    int n = x.rows();

    mon[0].commence_op();
    ArrayXf en = filtre_energie->step(x.abs2());
    // TODO : ligne à retard sur le signal d'énergie, pour synchroniser avec corr
    mon[0].fin_op();

    if(config.mode == DetecteurConfig::MODE_OLA)
      en = lar1->step(en);

    // corr = résultat du produit de corrélation entre le motif et le signal d'entrée
    mon[1].commence_op();
    ArrayXcf corr;
    ola->step(x, corr); // Retard = Ne échantillon avec un OLA, M échantillons avec un filtre
    mon[1].fin_op();


    // Pb : rien ne garantit que corr.rows() == n.

    tsd_assert_msg(corr.rows() == n, "Sortie OLA (corr) devrait faire {} échantillons, mais {}.", n, corr.rows());
    tsd_assert(en.rows() == n);

    mon[2].commence_op();

    float ratio = sqrt(1.0f * N) / sqrt(1.0f * M);

    if(config.mode == DetecteurConfig::MODE_RIF)
      ratio /= sqrt(1.0f * M);


    // Supprime les valeurs trop faibles
    float seuil = 1e-12f;

    corr = (seuil < corr.abs2()).select(corr, 0.0f).eval();

    y = ratio * (corr.abs2() / (en + 1e-20f)).sqrt();


    // TODO: Pb à régler plus proprement (premier buffer nul)
    //if((itr == 0) && (config.mode == DetecteurConfig::MODE_OLA))
      //y.head(Ne).setZero();

    /*ArrayXf y1(n);
    y1(0) = y(0);     // TODO
    y1(n-1) = y(n-1); // TODO
    for(auto i = 1; i + 1 < n; i++)
    {
      y1(i) = y(i-1) + y(i) + y(i+1);
    }*/

    // Pré-calcul pour l'érosion : sur chaque segment de largeur M, on
    // ne considère que le plus grand
    ArrayXf y2 = ArrayXf::Zero(n);
    for(auto i = 0; i < n; i += M)
    {
      int idx;
      auto mx = y.segment(i, min(M, n - i)).maxCoeff(&idx);
      y2(idx+i) = mx;
    }

    vector<int> lst = trouve(y2 > config.seuil), lst2;

    if(pic_final_a_traiter)
    {
      pic_final_a_traiter = false;
      lst2.push_back(-1);
    }

    // Attention : ci-dessous, coût quadratique en fonction du nombre de détections
    // (pb si seuil faible)
    // Il faut faire une sorte d'érosion
    for(auto idx: lst)
    {
      //msg("Cand : {} ({})", idx, y(idx));
      bool ok = true;
      for(auto idx2: lst)
      {
        if((y(idx2) > y(idx)) && (abs(idx - idx2) < M))
        {
          ok = false;
          break;
        }
      }
      if(ok)
        lst2.push_back(idx);
    }

    if(config.debug_actif)
    {
      Figures f;
      //f.subplot().plot(x1, "", "Buffer d'avant");
      f.subplot().plot(x, "", "x");
      f.subplot().plot(en.sqrt(), "", "Energie (sqrt)");
      f.subplot().plot(corr.abs(), "", "Corrélation");

      if((lst2.size() > 0) && (lst2[0] > 10) && (lst2[0] + 10 < n))
      {
        f.subplot().plot(y.segment(lst2[0]-10, 20), "", "Corrélation (zoom)");
      }


      auto s = f.subplot();
      s.plot(y, "", "Corrélation norm.");
      for(auto idx: lst2)
      {
        if(idx != -1)
          s.plot(idx, y(idx), "ro");
      }
      f.subplot().plot(y2, "", "Corrélation norm (érosion).");
      f.afficher(fmt::format("DBG Corr, itr {}", itr));
    }



    for(auto idx: lst2)
    {

      tsd_assert(idx <= (int) (n-1));
      tsd_assert(idx >= -1);

      //infos("Detection motif ok, score = %.2f, seuil = %.2f (energie = %f, row = %d).", score_max, config.seuil, en(maxRow), maxRow);
      Detection det;

      // Position fractionnaire
      float δ = 0;
      cfloat c0, c1, c2;
      float ac0, ac1, ac2;

      if(idx == -1)
      {
        det                = pic_final;
        det.position      -= dernier_n;
        det.position_prec -= dernier_n;
      }
      else
      {
        det.score     = y(idx);
        det.position  = idx - delais_corr;
        // position entre -Ne et -1 (sauf si cnt_ech != 0)
        det.θ         = arg(corr(idx));

        // Solution non optimale qui fait rentrer aussi le bruit
        //det.gain          = sqrt(en(idx)) / norme_motif;

        // PB si position fractionnaire : le gain n'est pas exact, e.g. erreur de 1 %
        // -> se traduit par une erreur d'estimation du bruit
        det.gain           = abs(corr(idx)) / (norme_motif / sqrt(N));
        det.position_prec  = det.position;

        // Position fractionnaire
        c1 = corr(idx);
        ac1 = y(idx);
      }

      if(idx == -1)
      {
        ac0 = alc0;
        c0  = lc0;
        ac1 = alc;
        c1  = lc;
        ac2 = y(0);
        c2  = corr(0);
        if(ac1 < ac2)
        {
          VERB(msg("Dernier échan : pas le pic.");)
          goto suite;
        }
      }
      else if(idx == 0)
      {
        ac0 = alc;
        c0  = lc;
        ac1 = y(idx);
        c1 = corr(idx);
        ac2 = y(1);
        c2  = corr(1);
        // Si le dernier échantillons du buffer précédant était le pic.
        if(ac1 < ac0)
        {
          VERB(msg("Premier échan : pas le pic.");)
          goto suite;
        }
      }
      else if(idx == n-1)
      {
        // Ceci sera traité au début de l'itération suivante,
        // car on a besoin d'un échantillon de plus pour faire l'interpolation quadratique
        VERB(msg("idx = n-1 -> pic final");)
        pic_final_a_traiter = true;
        pic_final = det;
        goto suite;
      }
      else
      {
        ac0 = y(idx-1);
        c0  = corr(idx-1);
        ac1 = y(idx);
        c1 = corr(idx);
        ac2 = y(idx+1);
        c2  = corr(idx+1);
      }

      if(config.mode == DetecteurConfig::MODE_RIF)
      {
        c0 /= sqrt(1.0f * M);
        c1 /= sqrt(1.0f * M);
        c2 /= sqrt(1.0f * M);
      }

      // Interpolation quadratique pour améliorer la précision
      δ = qint_loc(ac0, ac1, ac2);
      if((δ < -0.5) || (δ > 0.5))
      {
        msg_avert("Echec interpolation quadratique: pics = {} {} {}, δ={}, idx={}, Ne={}",
            ac0, ac1, ac2, δ, idx, Ne);
        δ = clamp(δ, -0.5f, 0.5f);
      }
      det.position_prec += δ;

      {
        float r =  sqrt((float) N) / norme_motif;
        cfloat g2 = qint_val(c0, c1, c2, δ) * r;
        det.gain  = abs(g2);
        det.θ     = arg(g2);
      }

      ArrayXcf recu_theo = std::polar(det.gain, det.θ) * config.motif;

      // TODO: un peu arbitraire
      //if(abs(δ) > 0.01)
      //{
      recu_theo = delais(recu_theo, δ);
      //}


      // On a besoin ici d'une sorte de délais variable
      // Donc, corr = [....... idx ...... ]
      //       x    = [.idx-D ........... ]

      // Délais maximum de la ligne à retard : dépends de n !

      ArrayXcf recu_vraiment(M);

      int id = idx - delais_corr;

      // Note : on suppose delais_corr >= M (pas besoin d'échantillon dans le futur)

      // Nb échantillons dispo dans le x en cours
      int nspl_in_current = M;
      if(id >= 0)
        nspl_in_current = M;
      else if((id < 0) && (id > -M))
        nspl_in_current = M + id;
      else
        nspl_in_current = 0;
      // Nb échantillons à prendre dans la ligne à retard
      int nspl_avant = M - nspl_in_current;

      //msg("nspl_in_curr={}, nspl_avant={}, M={}, lar.K={}, id={}", nspl_in_current, nspl_avant, M, lar.K, id);

      recu_vraiment.tail(nspl_in_current) = x.segment(max(id,0), nspl_in_current);
      if(nspl_avant < M)
        recu_vraiment.head(nspl_avant)      = lar.derniers(nspl_avant);
      else
        recu_vraiment.head(nspl_avant)      = lar.mem.segment(id + lar.K, nspl_avant); // Erreur ici
      ArrayXcf bruit     = recu_vraiment - recu_theo;

      auto var_bruit  = bruit.segment(1, M-2).abs2().mean();
      auto var_signal = pow(det.gain * norme_motif, 2.0f) / M;
      det.σ_noise = sqrt(var_bruit);
      det.SNR_dB  = pow2db(var_signal / var_bruit);


      if(config.debug_actif)
      {
        Figures f;
        f.subplot().plot(config.motif, "", "motif");
        f.subplot().plot(recu_theo, "", "recu_theo");
        f.subplot().plot(recu_vraiment, "", "recu_vraiment");
        f.subplot().plot(bruit, "r-", "bruit");
        f.afficher(fmt::format("DBG Det, itr {}", itr));
        msg("Détection : {}", det);
      }


      if(config.gere_detection)
        config.gere_detection(det);
    }

    suite:

    lar.step(x);

    mon[2].fin_op();

    alc0 = y(n-2);
    lc0  = corr(n-2);

    alc = y(n-1);
    lc  = corr(n-1);
    itr++;

    dernier_n = n;
  }


};


ostream &operator <<(ostream &os, const Detection &det)
{
  os << fmt::format("Détection : score={:.3f}, pos={} ({:.3f}), gain={:.5e}, θ={:.1f}°, σ={:.2e}, SNR={:.1f} dB.",
      det.score, det.position, det.position_prec, det.gain, rad2deg(det.θ), det.σ_noise, det.SNR_dB);
  return os;
}

sptr<Detecteur> détecteur_création(const DetecteurConfig &config)
{
  return make_shared<DetecteurImpl>(config);
}


}

