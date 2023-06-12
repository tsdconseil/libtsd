#include "tsd/tsd-all.hpp"

#define VERB(AA)

namespace tsd::fourier {

using namespace std;

// Interpolation quadratique,
// d'après https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
static float qint_loc(float ym1, float y0, float yp1)
{
  retourne (yp1 - ym1) / (2 * (2*y0 - yp1 - ym1));
}

// δ = position fractionnaire du max.
template<typename T>
static T qint_val(T ym1, T y0, T yp1, float δ)
{
  retourne y0 - (ym1 - yp1) * δ * ((T) 0.25);
}


// Ligne à retard spéciale, capable de retrouver les K derniers échantillons à tout moment
template<typename T>
struct LigneARetardExt
{
  entier K = 0;
  Vecteur<T> mem;

  void configure(entier K)
  {
    this->K = K;
    mem.setZero(K);
  }

  void step(const Vecteur<T> &x)
  {
    entier n = x.rows();
    si(n == K)
      mem = x;
    sinon si(n > K)
      mem = x.tail(K);
    sinon
    {
      mem.head(K-n) = mem.tail(K-n);
      mem.tail(n) = x;
    }
  }
  const Vecteur<T> &get_last_K() const
  {
    retourne mem;
  }
  Vecteur<T> derniers(entier n) const
  {
    retourne mem.tail(n);
  }
  Vecteur<T> get(entier delais, entier n) const
  {
    si((delais >= K) || (n >= delais))
      échec("Ligne à retard ext : dépassement, K={}, delais={}, n={}", K, delais, n);
    retourne mem.segment(K - delais - 1, n);
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

  entier itr = 0, dernier_n = 0;

  /** Ne = paquets d'entrée
   *  N  = dimension FFT (Ne + Nz), valable seulement en mode OLA
   *  M  = dimension motif
   */
  entier Ne = 0, N = 1, M = 0;

  /** TFD du motif, pré-calculée (seulement en mode OLA). */
  Veccf T_motif;

  /** Norme 2 du motif */
  float norme_motif;

  /** Motif, normalisé */
  Veccf motif;

  bouléen pic_final_a_traiter = non;

  Detection pic_final;

  MoniteursStats moniteurs()
  {
    retourne {{mon[0].stats(), mon[1].stats(), mon[2].stats()}};
  }

  DetecteurImpl(const DetecteurConfig &config): mon(3)
  {
    this->config = config;
    mon[0].nom() = "fft-corr/energie";
    mon[1].nom() = "fft-corr/ola";
    mon[2].nom() = "fft-corr/norm";
    configure_impl(config);
  }

  void configure_impl(const DetecteurConfig &config)
  {
    pic_final_a_traiter = non;
    itr                 = 0;

    norme_motif = config.motif.norme(); // sqrt(abs2(config.motif).somme());

    // Normalisation énergie du motif
    motif = config.motif / norme_motif;

    M = motif.rows();

    filtre_energie = filtre_mg<float,double>(M);

    Ne = config.Ne;

    si(config.Ne == 0)
    {
      float C;
      entier Nf, Nz;
      // TODO: retour multiples
      ola_complexité_optimise(M, C, Nf, Nz, Ne);
      msg("FFT corrélateur : calcul auto Ne optimal : M={} --> Ne={},Nz={},Nf=Ne+Nz={},C={} FLOPS/ech",
          M, Ne, Nz, Nf, C);
    }

    si(config.mode == DetecteurConfig::MODE_OLA)
    {
      FiltreFFTConfig ola_config;
      ola_config.nb_zeros_min             = M - 1;
      ola_config.dim_blocs_temporel       = Ne;
      ola_config.traitement_freq =
          [&](Veccf &X)
          {
            assertion(X.rows() == T_motif.rows());
            assertion(X.rows() == N);
            X *= T_motif.conjugate();
          };

      tie(ola, N) = filtre_fft(ola_config);

      Veccf tmp;
      tmp.setZero(N);
      assertion(M * 2 <= N);

      tmp.head(M) = motif;

      delais_corr = Ne;

      T_motif = fft(tmp);
      lar1 = ligne_a_retard<float>(delais_corr - M + 1);
    }
    sinon
    {
      msg("Détecteur : mode filtre RIF.");
      // Détection si le motif est réel -> filtre réel
      // (pour plus d'efficacité)

      // TODO: définir et utiliser une méthode .est_réel()
      si(abs(imag(motif)).somme() / abs(motif).somme() < 1e-7)
      {
        msg("  Motif réel détecté.");
        ola = filtre_rif<float, cfloat>(real(motif.reverse()));
      }
      sinon
      {
        msg("  Motif complexe détecté.");
        ola = filtre_rif<cfloat, cfloat>(motif.reverse().conjugate());
      }

      delais_corr = M - 1;
    }

    lar.configure(delais_corr+1);

    msg("Corrélateur FFT: Ne (dim blocs) = {}, M (dim motif) = {}, N (dim fft) = {}.", Ne, M, N);
    msg("  norme motif = {}", norme_motif);
  }


  cfloat lc = 0.0f, lc0 = 0.0f;
  float alc = 0, alc0 = 0;
  entier delais_corr = 0;

  void step(const Veccf &x, Vecf &y)
  {
    soit &config = Configurable<DetecteurConfig>::lis_config();

    soit n = x.rows();

    mon[0].commence_op();
    soit en = filtre_energie->step(abs2(x));
    mon[0].fin_op();

    // Ligne à retard sur le signal d'énergie, pour synchroniser avec corr
    si(config.mode == DetecteurConfig::MODE_OLA)
      en = lar1->step(en);

    // corr = résultat du produit de corrélation entre le motif et le signal d'entrée
    mon[1].commence_op();
    Veccf corr;
    // Retard = Ne échantillon avec un OLA, M échantillons avec un filtre
    ola->step(x, corr);
    mon[1].fin_op();


    // Pb : rien ne garantit que corr.rows() == n.
    assertion_msg(corr.rows() == n, "Sortie OLA (corr) devrait faire {} échantillons, mais {}.", n, corr.rows());
    assertion(en.rows() == n);

    mon[2].commence_op();

    soit ratio = sqrt(1.0f * N) / sqrt(1.0f * M);

    //si(config.mode == DetecteurConfig::MODE_RIF)
      //ratio /= sqrt(1.0f * M);


    // Supprime les valeurs trop faibles
    soit seuil = 1e-12f;

    //corr = (seuil < abs2(corr)).select(corr, 0.0f).eval();
    soit n2 = corr.dim();
    pour(auto i = 0; i < n2; i++)
      si(abs(corr(i)) <= sqrt(seuil))
        corr(i) = 0;

    y = ratio * sqrt(abs2(corr) / (en + 1e-20f));


    // TODO: Pb à régler plus proprement (premier buffer nul)
    //si((itr == 0) && (config.mode == DetecteurConfig::MODE_OLA))
      //y.head(Ne).setZero();

    /*ArrayXf y1(n);
    y1(0) = y(0);     // TODO
    y1(n-1) = y(n-1); // TODO
    pour(auto i = 1; i + 1 < n; i++)
    {
      y1(i) = y(i-1) + y(i) + y(i+1);
    }*/

    // Pré-calcul pour l'érosion : sur chaque segment de largeur M, on
    // ne considère que le plus grand
    soit y2 = Vecf::zeros(n);
    pour(auto i = 0; i < n; i += M)
    {
      soit [mx, idx] = y.segment(i, min(M, n - i)).max();
      y2(idx+i) = mx;
    }

    vector<entier> lst = trouve(y2 > config.seuil), lst2;

    si(pic_final_a_traiter)
    {
      pic_final_a_traiter = non;
      lst2.push_back(-1);
    }

    // Attention : ci-dessous, coût quadratique en fonction du nombre de détections
    // (pb si seuil faible)
    // Il faut faire une sorte d'érosion
    pour(auto idx: lst)
    {
      //msg("Cand : {} ({})", idx, y(idx));
      soit ok = oui;
      pour(auto idx2: lst)
      {
        si((y(idx2) > y(idx)) && (abs(idx - idx2) < M))
        {
          ok = non;
          break;
        }
      }
      si(ok)
        lst2.push_back(idx);
    }

    si(config.debug_actif)
    {
      Figures f;
      //f.subplot().plot(x1, "", "Buffer d'avant");
      f.subplot().plot(x, "", "x");
      f.subplot().plot(sqrt(en), "", "Energie (sqrt)");
      f.subplot().plot(abs(corr), "", "Corrélation");

      si((lst2.size() > 0) && (lst2[0] > 10) && (lst2[0] + 10 < n))
        f.subplot().plot(y.segment(lst2[0]-10, 20), "", "Corrélation (zoom)");


      soit s = f.subplot();
      s.plot(y, "", "Corrélation norm.");
      pour(auto idx: lst2)
      {
        si(idx != -1)
          s.plot(idx, y(idx), "ro");
      }
      f.subplot().plot(y2, "", "Corrélation norm (érosion).");
      f.afficher(sformat("DBG Corr, itr {}", itr));
    }

    pour(auto idx: lst2)
    {
      assertion(idx <= (entier) (n-1));
      assertion(idx >= -1);

      //msg("Detection motif ok, score = %.2f, seuil = %.2f (energie = %f, row = %d).", score_max, config.seuil, en(maxRow), maxRow);
      Detection det;

      // Position fractionnaire
      float δ = 0;
      cfloat c0, c1, c2;
      float ac0, ac1, ac2;

      si(idx == -1)
      {
        det                = pic_final;
        det.position      -= dernier_n;
        det.position_prec -= dernier_n;
      }
      sinon
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

      si(idx == -1)
      {
        ac0 = alc0;
        c0  = lc0;
        ac1 = alc;
        c1  = lc;
        ac2 = y(0);
        c2  = corr(0);
        si(ac1 < ac2)
        {
          VERB(msg("Dernier échan : pas le pic.");)
          goto suite;
        }
      }
      sinon si(idx == 0)
      {
        ac0 = alc;
        c0  = lc;
        ac1 = y(idx);
        c1 = corr(idx);
        ac2 = y(1);
        c2  = corr(1);
        // si le dernier échantillons du buffer précédant était le pic.
        si(ac1 < ac0)
        {
          VERB(msg("Premier échan : pas le pic.");)
          goto suite;
        }
      }
      sinon si(idx == n-1)
      {
        // Ceci sera traité au début de l'itération suivante,
        // car on a besoin d'un échantillon de plus pour faire l'interpolation quadratique
        VERB(msg("idx = n-1 -> pic final");)
        pic_final_a_traiter = oui;
        pic_final = det;
        goto suite;
      }
      sinon
      {
        ac0 = y(idx-1);
        c0  = corr(idx-1);
        ac1 = y(idx);
        c1 = corr(idx);
        ac2 = y(idx+1);
        c2  = corr(idx+1);
      }

      /*si(config.mode == DetecteurConfig::MODE_RIF)
      {
        c0 /= sqrt(1.0f * M);
        c1 /= sqrt(1.0f * M);
        c2 /= sqrt(1.0f * M);
      }*/

      // Interpolation quadratique pour améliorer la précision
      δ = qint_loc(ac0, ac1, ac2);
      si((δ < -0.5) || (δ > 0.5))
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

      soit recu_theo = config.motif * std::polar(det.gain, det.θ);

      // TODO: un peu arbitraire
      //si(abs(δ) > 0.01)
      //{
      recu_theo = délais(recu_theo, δ);
      //}


      // On a besoin ici d'une sorte de délais variable
      // Donc, corr = [....... idx ...... ]
      //       x    = [.idx-D ........... ]

      // Délais maximum de la ligne à retard : dépends de n !

      Veccf recu_vraiment(M);

      entier id = idx - delais_corr;

      // Note : on suppose delais_corr >= M (pas besoin d'échantillon dans le futur)

      // Nb échantillons dispo dans le x en cours
      soit nspl_in_current = M;
      si(id >= 0)
        nspl_in_current = M;
      sinon si((id < 0) && (id > -M))
        nspl_in_current = M + id;
      sinon
        nspl_in_current = 0;
      // Nb échantillons à prendre dans la ligne à retard
      soit nspl_avant = M - nspl_in_current;

      //msg("nspl_in_curr={}, nspl_avant={}, M={}, lar.K={}, id={}", nspl_in_current, nspl_avant, M, lar.K, id);

      recu_vraiment.tail(nspl_in_current) = x.segment(max(id,0), nspl_in_current);
      si(nspl_avant < M)
        recu_vraiment.head(nspl_avant)      = lar.derniers(nspl_avant);
      sinon
        recu_vraiment.head(nspl_avant)      = lar.mem.segment(id + lar.K, nspl_avant); // Erreur ici

      soit bruit      = recu_vraiment - recu_theo;
      soit var_bruit  = abs2(bruit.segment(1, M-2)).moyenne();
      soit var_signal = pow(det.gain * norme_motif, 2.0f) / M;
      det.σ_noise = sqrt(var_bruit);
      det.SNR_dB  = pow2db(var_signal / var_bruit);

      si(config.debug_actif)
      {
        Figures f;
        f.subplot().plot(config.motif,  "",   "motif");
        f.subplot().plot(recu_theo,     "",   "recu_theo");
        f.subplot().plot(recu_vraiment, "",   "recu_vraiment");
        f.subplot().plot(bruit,         "r-", "bruit");
        f.afficher(sformat("DBG Det, itr {}", itr));
        msg("Détection : {}", det);
      }

      si(config.gere_detection)
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
  os << sformat("Détection : score={:.3f}, pos={} ({:.3f}), gain={:.5e}, θ={:.1f}°, σ={:.2e}, SNR={:.1f} dB.",
      det.score, det.position, det.position_prec, det.gain, rad2deg(det.θ), det.σ_noise, det.SNR_dB);
  retourne os;
}

sptr<Detecteur> détecteur_création(const DetecteurConfig &config)
{
  retourne make_shared<DetecteurImpl>(config);
}


}

