#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include "tsd/fourier/cqt.hpp"

#include <cmath>



using namespace tsd;
using namespace tsd::fourier;
using namespace tsd::vue;

namespace tsd::tf::cqt {

using namespace std;

struct CQTNoyaux
{
  // Initialization of a CQT kernel structure.
  //
  // Calling Sequence
  // cqtk = cqtk_init(fs, fmin, fmax, gamm, [Q, kprec]);
  //
  // Parameters
  // fs: Sample frequency, in Hz
  // fmin: Minimum frequency of interest (must be > 0)
  // fmax: Maximum frequency of interest (must be > fmin)
  // gamm: Desired ratio between two successive freqency bins
  // Q: Quality factor (default value: 34)
  // kprec: Kernel precision: between 0 (lowest precision) and 1 (highest precision). Default value is 0.99.
  // cqtk: Resulting CQT kernel structure
  //
  // Description
  // This function compute the CQT kernel (in frequency domain) given the specified time / frequency resolution. These kernels can then been used thereafter to compute the time / frequency analysis of one or several signals.
  //
  // The parameter <varname>kprec</varname> (kernel precision) can be used to speed up the following computations (in <link linkend="cqtk_process">cqtk_process</link>) by approximating more or less the frequency domain kernels.
  //
  // <refsection><title>Example</title></refsection>
  // <para>In this example, we create a frequency interval from A1 (55 Hz) to A8 (7.04 kHz), with a progression of one tone between two consecutive frequencies.</para>
  // <programlisting>
  //   // A1 is 55 Hz
  //   // A8 os 7040 Hz
  //   // A tone is a twelth of octave
  //   cqtk = cqtk_init(fs = 8e3, fmin = 55, fmax = 7.04e3, gamm = 2^(1/12));
  //   // Display the frequency bins
  //   clf(); plot2d3(cqtk.freqs);
  // </programlisting>
  // <mediaobject><imageobject><imagedata fileref="ex_cqtk_init.png" format="PNG"/></imageobject><caption><para>Frequency bins (frequency value vs bin number)</para></caption></mediaobject>
  // <para>As can be seen on the figure above, the frequency bins pour CQT analysis are geometrically distribued (constant ratio between two consecutive frequencies).</para>
  //
  // See also
  //  cqtk_plot
  //  cqtk_save
  //  cqtk_load
  //
  // Authors
  //  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/scriptscilab/cqt">http://www.tsdconseil.fr/log/scriptscilab/cqt</ulink>

  string id;
  float fs, fmin, fmax, γ;
  float Q;
  entier nfreqs, N;
  float précision_noyau;
  Vecf freqs;
  Veci ktime, kimin, ksize;
  Tabcf noyaux;

  void configure(const CQTConfig &config)//float fs, float fmin, float fmax, float gamm, float Q = 34, float kernel_precision = 0.99)
  {
    this->id      = "CQT kernel structure";
    this->fs      = config.fs;
    this->fmin    = config.fmin;
    this->fmax    = config.fmax;
    this->γ       = config.γ;
    this->Q       = config.Q;

    // fmin * γ^nfreqs = fmax
    // => log(fmax/fmin) = nfreqs * log(γ)
    nfreqs = (entier) ceil(log(fmax / fmin) / log(γ));

    this->précision_noyau = 0.99;

    si((précision_noyau <= 0) || (précision_noyau >= 1))
      echec("CQT: kernel_precision must be in ]0,1[ (typical : 0.98). Argument : {}.", précision_noyau);

    si(γ <= 1)
      echec("CQT: γ (ratio of 2 successive frequencies) must be > 1 (γ = {}).", γ);

    si(Q <= 0)
      echec("CQT: Quality factor (Q) must be > 0 (Q = {}).", Q);

    si(fs <= 0)
      echec("CQT: fs must be > 0 (fs = {}).", fs);

    si((fmin <= 0) || (fmin > fs / 2))
      echec("CQT: fmin must be strictly > 0 Hz and <= fs/2 (fmin = {}).", fmin);

    // (1) Compute FFT minimal size
    entier ideal_N = (entier) ceil((fs * Q) / fmin);
    N = 1;
    tantque(N < ideal_N)
      N *= 2;

    fmax = fmin * pow(γ, (float) nfreqs);
    msg("CQT: fs = {:.0f} Hz, fmin = {:.2f} Hz, fmax = {:.2f} Hz, γ = {:.3f}, nfreqs = {}, Q = {:.2f}.",
        fs, fmin, fmax, γ, nfreqs, Q);
    msg("Ideal FFT window size = {}, practical = {}.", ideal_N, N);
    //printf("Freq. interval: %.1f Hz - %.1f Hz.\n", fmin, fmin * gamm^nfreqs);


    freqs = Vecf::zeros(nfreqs);
    ksize = kimin = ktime = Veci::zeros(nfreqs);

    freqs(0) = fmin;
    pour(auto i = 1; i < nfreqs; i++)
      freqs(i) = freqs(i-1) * γ;

    //entier lp = 0;
    // (2) Initialize kernels
    msg("initialisation des noyaux...");
    pour(auto i = nfreqs-1; i >= 0; i--)
      init_noyau(i);
    msg("fait.");
  }

  // Taille max pour les noyaux : dÃ©terminÃ© par la plus haute frÃ©quence
  void init_noyau(entier i)
  {
    soit frequence = freqs(i);
    soit kernel_time = (entier) ceil((fs * Q) / frequence);
    si(kernel_time > N-1)
    {
      msg_avert("Attention, i = {}, ktime = {}, N = {}.", i, kernel_time, N);
      kernel_time = N-1;
    }
    si((kernel_time % 2) == 0)
      kernel_time++;

    ktime(i) = kernel_time;

    // Compute spatial response, centered
    soit ol = cos(2*π*(frequence/fs)*linspace(1,N,N));

    // Fenêtre de Hamming
    soit fen = tsd::filtrage::fenetre("hm", kernel_time, non);
    //float α = 0.5;
    //ArrayXf wnd = alpha - (1.0 - alpha) * cos((2 * π * linspace(1,kernel_time,kernel_time))/kernel_time);
    soit tmp = Vecf::zeros(N);

    soit k2 = kernel_time / 2;

    tmp.segment(N/2-k2, kernel_time) = ol.segment(N/2-k2,kernel_time) * fen;

    // FFT du noyau
    soit X = fft(tmp);

    // Localisation du noyau
    soit Xa = abs2(X);
    soit total = Xa.somme();

    soit jmax = Xa.head(N/2).index_max();

    // Recherche j tel que
    // sum(Xa,j = jmax-j:jmax+j) > total * kernel_precision
    // <=> Sp(jmax+j) - Sp(jmax-j) > total * kernel_precision
    // Appriximation:
    // <~> Sp(jmax+j) - Sp(jmax) > total * kernel_precision / 2
    // j = find(Sp(jmax:$) - Sp(jmax) > total * kernel_precision / 2, 1);

    soit Sp = cumsum(Xa);
    float pv;
    si(jmax == 1)
      pv = 0;
    sinon
      pv = Sp(jmax-1);

    entier j2;
    entier lc = trouve_premier(Sp.segment(jmax,N/2-jmax) > (Xa.segment(jmax,N/2-jmax)).somme() * précision_noyau + pv);
    si(lc == -1)
      j2 = N/2;
    sinon
      j2 = (jmax-1) + lc;
    entier j1 = trouve_premier(Sp.head(jmax) > (Xa.head(jmax)).somme() * (1-précision_noyau));
    si(j1 == -1)
      j1 = 1;


    // TODO
    si(j2 - j1 < 4)
      j2 = j1 + 4;

    entier kernel_index_min  = j1;//jmax - j;

    //si(kernel_index_min < 1) then
      //  kernel_index_min = 1;
    //end;
    kimin(i) = kernel_index_min;

    entier kernel_size = j2 - j1 + 1;//2 * j + 1;
    //    si(kernel_index_min + kernel_size >= N + 1) then
    //        kernel_size = N - kernel_index_min;
    //    end;
    ksize(i) = kernel_size;

    //    si(kernel_size <= 1)
    //        printf("ksize = 1 ?\n");
    //        pause;
    //    end;

    si(i == nfreqs - 1)
    {
      msg("Allocation des noyaux : {} * {} = {:.1f} ko...",
          kernel_size, nfreqs, kernel_size * nfreqs * 8.0 / 1024);
      noyaux = Tabcf::zeros(kernel_size, nfreqs);
    }

    // Store kernel
    noyaux.col(i).head(kernel_size) = X.segment(j1,j2-j1+1) / sqrt(total);
  }


  void plot()
  {
  // Plotting of the CQT kernels
  //
  // This function plots the CQT kernels, in the frequency domain.
  // For each frequency bin, a CQT kernel is stored as the Fourier transform of the corresponding CQT kernel in time domain (in time domain, a kernel is a windowed (pour compact support) sinusoidal wave).
  // Example
  // cqtk = cqtk_init(fs = 8e3, fmin = 100, fmax = 4e3, gamm = 2^(1/12));
  // cqtk_plot(cqtk);
  // As can be seen in the figure above, the CQT kernels are wider (in frequency) pour higher frequencies, which is consistent with the "Constant Q" objective.

    soit mv4 = 256;//1024;

    soit d = 1.0f;
    si(N / 2 > mv4)
      d = floor(N / (2 * mv4));


    soit np = (entier) floor(N / (2*d));
    Tabf nymap = Tabf::zeros(np, nfreqs);

    msg("plot noyaux cqt : N={}, d={}, N/(2d)={}, nfreqs={}", N, d, np, nfreqs);

    pour(auto i = 0; i < nfreqs; i++)
    {
      soit tmp1 = abs(noyaux.col(i).head(ksize(i)));
      soit tmp  = rééchan(tmp1,1.0f/d);
      soit ns = (entier) floor(ksize(i)/d),
           nm = (entier) floor(kimin(i)/d);
      si(nm < 1)
        nm = 1;

      si((ns > 0) && (tmp.rows() > 0))
      nymap.col(i).segment(nm,/*ns*/tmp.rows()) = tmp;
    }

    Figure f;
    f.plot_img({{0, 1}, Dimf{fs/2, nfreqs}}, nymap);
    f.afficher("Noyaux CQT");
  }
};


struct CQT::Impl
{
  CQTNoyaux cqtk;

  entier abs_position = 0, nb_data_attendu = 0;
  Veci kposition, kaposition;
  Vecf cirbuffer;
  sptr<tsd::fourier::FFTPlan> plan_fft;

  struct Sortie
  {
    entier id_freq;
    // Attention, temps = nombre de symboles
    float temps, valeur;
  };

  vector<Sortie> sortie;

  void configure(const CQTConfig &config)
  {
    cqtk.configure(config);
    // Zero-padding en début d'analyse par N/2 zéros
    abs_position   = -cqtk.N / 2;
    nb_data_attendu = cqtk.N / 2;

    kposition   = (cqtk.N / 2) * Veci::ones(cqtk.nfreqs);
    kaposition  = Veci::zeros(cqtk.nfreqs);
    cirbuffer   = Vecf::zeros(cqtk.N);

    plan_fft = tfrplan_création(cqtk.N, oui);
  }

  void step(const Vecf &x)
  {
    soit lon      = x.rows();
    soit décalage = 1;

    tantque(lon >= nb_data_attendu)
    {
      soit x1 = x.segment(décalage, nb_data_attendu);
      lon      -= nb_data_attendu;
      décalage += nb_data_attendu;
      step1(x1);
    }
  }




  void step1(const Vecf &x)
  {
    soit N = cqtk.N;

    tsd_assert(x.rows() == nb_data_attendu);

    // (1) Data rotation
    //      to remove      to keep
    //  |----------------|---------|
    //        n_waited      n_keep
    cirbuffer.head(N - nb_data_attendu)   = cirbuffer.tail(N - nb_data_attendu);
    cirbuffer.tail(nb_data_attendu) = x;

    // (2) Change to frequency domain
    //ArrayXcf X = fft(cirbuffer);

    // TODO : utiliser plan rfft
    soit X = plan_fft->step(cirbuffer);

    soit minimal_position = N;
    // (3) Kernels correlations
    pour(auto i = 0; i < cqtk.nfreqs; i++)
    {
      soit ksize    = cqtk.ksize(i),
           kmin     = cqtk.kimin(i),
           ktime    = cqtk.ktime(i),
           kpos     = kposition(i),
           kapos    = kaposition(i);

      soit kerncx = cqtk.noyaux.col(i).head(ksize) * (X.segment(kmin, ksize)).conjugate();
      soit dt = (entier) floor(ktime / 2);

      tantque(kpos + ktime <= N)
      {
        // Time shift, in samples, from the middle of the frame
        soit décalage = (-kpos + 0.5f*ktime) / N - 0.5f;

        // OL
        soit ol0 = std::polar(1.0f, 2 * π_f * décalage);
        // Mixer pour + freqs
        soit s0 = std::polar(1.0f, 2 * π_f * kmin * décalage);

        // Trop consommateur
        //float cor = abs((tsd::polar(shift*posindex) * kerncx).sum());
        cdouble correl = 0;
        pour(auto j = 0; j < ksize; j++)
        {
          /* shift kernel in time domain */
          correl += kerncx(j)  * s0;
          /* Update mixer */
          s0 *= ol0;
        }

        sortie.push_back({i, (float) kapos, (float) abs(correl)});

        // Shift kernel (with overlap 1/2)
        kpos  += dt;
        kapos += dt;
      } // end tantque pos
      si(kpos < minimal_position)
        minimal_position = kpos;
      kposition(i)  = kpos;
      kaposition(i) = kapos;
    } // end pour i=1::nfreqs

    // Nb data waited = minimal_position
    nb_data_attendu = minimal_position;
    kposition      -= minimal_position;
  }

  tuple<Vecf, Vecf, Tabf> interpolation(float ofs)
  {
    // Nombre total d'échantillons (irréguliers)
    soit ns = (entier) sortie.size();
    // Nombre d'échantillons d'entrée / échantillon de sortie
    soit oprd = cqtk.fs / ofs;

    msg("interpolation : ofs = {} Hz, nptr = {}.", ofs, sortie.size());

    soit lvalue = Vecf::zeros(cqtk.nfreqs);
    soit lindex = lvalue.clone();

    // Index d'écriture, pour chaque fréquence
    soit décalages = Veci::zeros(cqtk.nfreqs);

    // Durée totale de l'analyse (en secondes)
    soit duree = sortie[ns-1].temps / cqtk.fs;

    // Nombre d'échantillons à sortie (par fréquence)
    soit onmax = (entier) ceil(duree * ofs);

    msg("CQT interpolation : duration = {} samples = {} seconds.", sortie[ns-1].temps, duree);
    msg("CQT interpolation : creation of a matrix of size {} * {}.", onmax, cqtk.nfreqs);

    soit A = Tabf::zeros(onmax, cqtk.nfreqs);
    // Echantillons classés par indices de fréquence croissants
    pour(auto i = 0; i < ns; i++)
    {
      soit &spl = sortie[i];
      soit fr = spl.id_freq;
      si(fr >= cqtk.nfreqs)
        echec("ofreq invalid.");

      // Temps échantillon précédent (pour cette fréquence), et échantillon courant
      // (en nombre d'échantillons d'entrée)
      soit t0 = lindex(fr), t1 = spl.temps;
      // Valeur précédente et courante
      soit y0 = lvalue(fr), y1 = abs(spl.valeur);
      // Compléte la matrice de sortie
      tantque(décalages(fr) * oprd <= t1)
      {
        soit t = décalages(fr) * oprd;
        float val;
        si(t1 == t0)
            val = y1;
        sinon
            val = y0 + ((t - t0) / (t1 - t0)) * (y1 - y0);
        A(décalages(fr),fr) = val;
        décalages(fr)++;
      }
      lindex(fr) = t1;
      lvalue(fr) = y1;
    }
     soit t = linspace(0,A.rows()/ofs,A.rows());
     // [t,f,A]
     retourne {t, cqtk.freqs, A};
  }

};

CQT::CQT()
{
  impl = make_shared<Impl>();
}

void CQT::configure(const CQTConfig &config)
{
  impl->configure(config);
}

void CQT::step(const Vecf &x)
{
  impl->step(x);
}

void CQT::affiche_noyaux()
{
  impl->cqtk.plot();
}

tuple<Vecf, Vecf, Tabf> CQT::interpolation(float ofs)
{
  retourne impl->interpolation(ofs);
}



}

namespace tsd::tf {

  using namespace tsd::tf::cqt;

  Tabf periodogramme_cqt(const Vecf &x, float fs, float fmin, float fmax, float γ, float ofs)
  {
    CQTConfig config;
    CQT cqt;

    config.fs   = fs;
    config.fmin = fmin;
    config.fmax = fmax;
    config.γ    = γ;

    cqt.configure(config);
    //cqt.affiche_noyaux();
    cqt.step(x);
    soit [t, f, A] = cqt.interpolation(ofs);
    retourne A;
  }


}
