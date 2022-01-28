#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier/cqt.hpp"

using namespace tsd;
using namespace tsd::fourier;
using namespace tsd::vue;

namespace tsd::tf::cqt {

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
  // <para>As can be seen on the figure above, the frequency bins for CQT analysis are geometrically distribued (constant ratio between two consecutive frequencies).</para>
  //
  // See also
  //  cqtk_plot
  //  cqtk_save
  //  cqtk_load
  //
  // Authors
  //  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/scriptscilab/cqt">http://www.tsdconseil.fr/log/scriptscilab/cqt</ulink>

  std::string id;
  float fs, fmin, fmax, γ;
  float Q;
  int nfreqs, N;
  float kernel_precision;
  ArrayXf freqs;
  ArrayXi ktime, kimin, ksize;
  ArrayXXcf noyaux;

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
    nfreqs = (int) ceil(log(fmax / fmin) / log(γ));

    this->kernel_precision = 0.99;

    if((kernel_precision <= 0) || (kernel_precision >= 1))
      echec("CQT: kernel_precision must be in ]0,1[ (typical : 0.98). Argument : {}.", kernel_precision);

    if(γ <= 1)
      echec("CQT: γ (ratio of 2 successive frequencies) must be > 1 (γ = {}).", γ);

    if(Q <= 0)
      echec("CQT: Quality factor (Q) must be > 0 (Q = {}).", Q);

    if(fs <= 0)
      echec("CQT: fs must be > 0 (fs = {}).", fs);

    if((fmin <= 0) || (fmin > fs / 2))
      echec("CQT: fmin must be strictly > 0 Hz and <= fs/2 (fmin = {}).", fmin);

    // (1) Compute FFT minimal size
    int ideal_N = (int) ceil((fs * Q) / fmin);
    N = 1;
    while(N < ideal_N)
      N *= 2;

    fmax = fmin * std::pow(γ, (float) nfreqs);
    msg("CQT: fs = {:.0f} Hz, fmin = {:.2f} Hz, fmax = {:.2f} Hz, γ = {:.3f}, nfreqs = {}, Q = {:.2f}.",
        fs, fmin, fmax, γ, nfreqs, Q);
    msg("Ideal FFT window size = {}, practical = {}.", ideal_N, N);
    //printf("Freq. interval: %.1f Hz - %.1f Hz.\n", fmin, fmin * gamm^nfreqs);


    freqs = ArrayXf::Zero(nfreqs);
    ksize = kimin = ktime = ArrayXi::Zero(nfreqs);

    freqs(0) = fmin;
    for(auto i = 1; i < nfreqs; i++)
      freqs(i) = freqs(i-1) * γ;

    //int lp = 0;
    // (2) Initialize kernels
    msg("initialisation des noyaux...");
    for(auto i = nfreqs-1; i >= 0; i--)
    {
      //int percent_done = (int) (100 * (nfreqs-i+1) / nfreqs);
      kernel_init(i);
      /*if(int(percent_done / 20) <> int(lp / 20)) then
        lp = percent_done;
        printf("CQT kernel initialization %d/%d (%.d %% done)...\n", i, nfreqs, percent_done);
      end;*/
    }
    msg("fait.");
  }

  // Taille max pour les noyaux : dÃ©terminÃ© par la plus haute frÃ©quence
  void kernel_init(int i)
  {
    float frequency = freqs(i);
    int kernel_time = (int) ceil((fs * Q) / frequency);
    if(kernel_time > N-1)
    {
      msg_avert("Attention, i = {}, ktime = {}, N = {}.", i, kernel_time, N);
      kernel_time = N-1;
    }
    if((kernel_time % 2) == 0)
      kernel_time++;

    ktime(i) = kernel_time;

    // Compute spatial response, centered
    ArrayXf ol = (2*π*(frequency/fs)*linspace(1,N,N)).cos();

    // Fenêtre de Hamming
    ArrayXf wnd = tsd::filtrage::fenetre("hm", kernel_time, false);
    //float α = 0.5;
    //ArrayXf wnd = alpha - (1.0 - alpha) * cos((2 * π * linspace(1,kernel_time,kernel_time))/kernel_time);
    ArrayXf tmp = ArrayXf::Zero(N);

    int k2 = kernel_time / 2;

    tmp.segment(N/2-k2, kernel_time) = ol.segment(N/2-k2,kernel_time) * wnd;
    //tmp((N/2 - int(kernel_time/2)):(N/2 + int(kernel_time/2))) = ol((N/2 - int(kernel_time)/2):(N/2 + int(kernel_time/2))) .* wnd;

    // FFT du noyau
    ArrayXcf X = fft(tmp);

    // Localisation du noyau
    ArrayXf Xa = X.abs2();
    float total = Xa.sum();

    int jmax;
    Xa.head(N/2).maxCoeff(&jmax);
    // float emax = Xa.head(N/2).maxCoeff(&jmax);

    // Recherche j tel que
    // sum(Xa,j = jmax-j:jmax+j) > total * kernel_precision
    // <=> Sp(jmax+j) - Sp(jmax-j) > total * kernel_precision
    // Appriximation:
    // <~> Sp(jmax+j) - Sp(jmax) > total * kernel_precision / 2
    // j = find(Sp(jmax:$) - Sp(jmax) > total * kernel_precision / 2, 1);

    ArrayXf Sp = cumsum(Xa);
    float pv;
    if(jmax == 1)
      pv = 0;
    else
      pv = Sp(jmax-1);

    int j2;
    int lc = trouve_premier(Sp.segment(jmax,N/2-jmax) > (Xa.segment(jmax,N/2-jmax)).sum() * kernel_precision + pv);
    if(lc == -1)
      j2 = N/2;
    else
      j2 = (jmax-1) + lc;
    int j1 = trouve_premier(Sp.head(jmax) > (Xa.head(jmax)).sum() * (1-kernel_precision));
    if(j1 == -1)
      j1 = 1;


    //    if(i == 30)
    //        pause;
    //    end;
    //

    //      m1 = Xa(jmax:2*jmax-1);
    //      //m1 = Xa(jmax:N/2);
    //      m2 = Xa(jmax:-1:1);
    //      Sp = cumsum(m1+m2);
    //      j = find(Sp - Xa(jmax) > total * kernel_precision / 2, 1);

    // TODO
    if(j2 - j1 < 4)
      j2 = j1 + 4;

    int kernel_index_min  = j1;//jmax - j;

    //if(kernel_index_min < 1) then
      //  kernel_index_min = 1;
    //end;
    kimin(i) = kernel_index_min;

    int kernel_size = j2 - j1 + 1;//2 * j + 1;
    //    if(kernel_index_min + kernel_size >= N + 1) then
    //        kernel_size = N - kernel_index_min;
    //    end;
    ksize(i) = kernel_size;

    //    if(kernel_size <= 1)
    //        printf("ksize = 1 ?\n");
    //        pause;
    //    end;

    if(i == nfreqs - 1)
    {
      msg("Allocation des noyaux : {} * {} = {:.1f} ko...", kernel_size, nfreqs, kernel_size * nfreqs * 8.0 / 1024);
      noyaux = ArrayXXcf::Zero(kernel_size, nfreqs);
    }

    // Store kernel
    //cqtk.noyaux(1:kernel_size,i) = X(kernel_index_min:kernel_index_min+kernel_size-1) ./ sqrt(total);
    noyaux.col(i).head(kernel_size) = X.segment(j1,j2-j1+1) / sqrt(total);
  }


  void plot()
  {
  // Plotting of the CQT kernels
  //
  // Calling Sequence
  // cqtk_plot(cqtk)
  //
  // Parameters
  // cqtk: CQT kernel structure, as computed by <link linkend="cqtk_init">cqtk_init</link>
  //
  // Description
  // This function plots the CQT kernels, in the frequency domain.
  // For each frequency bin, a CQT kernel is stored as the Fourier transform of the corresponding CQT kernel in time domain (in time domain, a kernel is a windowed (for compact support) sinusoidal wave).
  // <refsection><title>Example</title></refsection>
  // <programlisting>
  // cqtk = cqtk_init(fs = 8e3, fmin = 100, fmax = 4e3, gamm = 2^(1/12));
  // cqtk_plot(cqtk);
  // </programlisting>
  // <mediaobject><imageobject><imagedata fileref="ex_cqtk_plot.png" format="PNG"/></imageobject><caption><para>CQT kernels (frequency domain)</para></caption></mediaobject>
  // As can be seen in the figure above, the CQT kernels are wider (in frequency) for higher frequencies, which is consistent with the "Constant Q" objective.
  //
  // See also
  //  cqtk_init
  //
  // Authors
  //  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/scriptscilab/cqt">http://www.tsdconseil.fr/log/scriptscilab/cqt</ulink>



    int mv4 = 256;//1024;

    float d = 1;
    if(N / 2 > mv4)
      d = floor(N / (2*mv4));


    int np = (int) floor(N / (2*d));
    ArrayXXf nymap = ArrayXXf::Zero(np, nfreqs);

    msg("plot noyaux cqt : N={}, d={}, N/(2d)={}, nfreqs={}", N, d, np, nfreqs);

    for(auto i = 0; i < nfreqs; i++)
    {
      ArrayXf tmp1 = (noyaux.col(i).head(ksize(i))).abs();
      ArrayXf tmp = reechan(tmp1,1.0f/d);
      int ns = (int) floor(ksize(i)/d);
      int nm = (int) floor(kimin(i)/d);
      if(nm < 1)
        nm = 1;

      if((ns > 0) && (tmp.rows() > 0))

      //if(ns != tmp.rows())
      //  msg_erreur("ns = {}, tmp:{}, d:{}, ksize:{}.", ns, tmp.rows(), d, ksize(i));
      nymap.col(i).segment(nm,/*ns*/tmp.rows()) = tmp;
    }

    Figure f;
    f.plot_img(0, fs/2, 1, nfreqs, nymap);
    f.afficher("Noyaux CQT");
    //mySgrayplot(linspace(0,cqtk.fs/2,floor(cqtk.N / (2*d))), linspace(0, cqtk.nfreqs, cqtk.nfreqs), abs(nymap));
    //xtitle("CQT kernels", "Frequency (Hz)", "Kernel number")
  }
};


struct CQT::Impl
{
  CQTNoyaux cqtk;

  int abs_position = 0, nb_data_attendu = 0;
  ArrayXi kposition, kaposition;
  ArrayXf cirbuffer;
  sptr<tsd::fourier::FFTPlanAbstrait> plan_fft;

  struct Sortie
  {
    int id_freq;
    // Attention, temps = nombre de symboles
    float temps, valeur;
  };

  std::vector<Sortie> sortie;

  void configure(const CQTConfig &config)
  {
    cqtk.configure(config);
    // Zero-padding en début d'analyse par N/2 zéros
    abs_position   = -cqtk.N / 2;
    nb_data_attendu = cqtk.N / 2;

    kposition   = (cqtk.N / 2) * ArrayXi::Ones(cqtk.nfreqs);
    kaposition  =  ArrayXi::Zero(cqtk.nfreqs);
    cirbuffer   = ArrayXf::Zero(cqtk.N);

    plan_fft = creation_fft_plan(cqtk.N, true);
  }

  void step(const ArrayXf &x)
  {
    int len = x.rows();
    int offset = 1;
    //float lp = 0;

    while(len >= nb_data_attendu)
    {
      //percent_done = int(100 * offset / length(x));
      //if(int(percent_done / 20) <> int(lp / 20)) then
        //lp = percent_done;
        //printf("CQT iteration: %.d %% done...\n", percent_done);
      //end;
      ArrayXf x1 = x.segment(offset, nb_data_attendu);
      len    -= nb_data_attendu;
      offset += nb_data_attendu;
      step1(x1);
    }
  }




  void step1(const ArrayXf &x)
  {
    int N = cqtk.N;

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
    ArrayXcf x2 = cirbuffer;
    ArrayXcf X;
    plan_fft->step(x2, X);

    int minimal_position = N;


    // (3) Kernels correlations
    for(auto i = 0; i < cqtk.nfreqs; i++)
    {
      int ksize    = cqtk.ksize(i);
      int kmin     = cqtk.kimin(i);
      int ktime    = cqtk.ktime(i);
      int kpos     = kposition(i);
      int kapos    = kaposition(i);

      ArrayXcf kerncx = cqtk.noyaux.col(i).head(ksize) * (X.segment(kmin, ksize)).conjugate();
      int dt = (int) floor(ktime / 2);

      while(kpos + ktime <= N)
      {
        // Time shift, in samples, from the middle of the frame
        float shift = (-kpos + 0.5*ktime) / N - 0.5f;

        // OL
        cfloat ol0 = std::polar(1.0f, 2 * π_f * shift);
        // Mixer for + freqs
        cfloat s0 = std::polar(1.0f, 2 * π_f * kmin * shift);

        // Trop consommateur
        //float cor = std::abs((tsd::polar(shift*posindex) * kerncx).sum());
        cdouble correl = 0;
        for(auto j = 0; j < ksize; j++)
        {
          /* shift kernel in time domain */
          correl += kerncx(j)  * s0;
          /* Update mixer */
          s0 *= ol0;
        }

        sortie.push_back({i, (float) kapos, (float) std::abs(correl)});

        // Shift kernel (with overlap 1/2)
        kpos  += dt;
        kapos += dt;
      } // end while pos
      if(kpos < minimal_position)
        minimal_position = kpos;
      kposition(i)  = kpos;
      kaposition(i) = kapos;
    } // end for i=1::nfreqs

    // Nb data waited = minimal_position
    nb_data_attendu = minimal_position;
    kposition     -= minimal_position;
  }

  std::tuple<ArrayXf, ArrayXf, ArrayXXf> interpolation(float ofs)
  {
    // Nombre total d'échantillons (irréguliers)
    int ns = sortie.size();
    // Nombre d'échantillons d'entrée / échantillon de sortie
    float oprd = cqtk.fs / ofs;

    msg("interpolation : ofs = {} Hz, nptr = {}.", ofs, sortie.size());

    ArrayXf lvalue = ArrayXf::Zero(cqtk.nfreqs);
    ArrayXf lindex = lvalue;

    // Index d'écriture, pour chaque fréquence
    ArrayXi offset = ArrayXi::Zero(cqtk.nfreqs);

    // Durée totale de l'analyse (en secondes)
    float duree = sortie[ns-1].temps / cqtk.fs;

    // Nombre d'échantillons à sortie (par fréquence)
    int onmax = (int) ceil(duree * ofs);

    msg("CQT interpolation : duration = {} samples = {} seconds.", sortie[ns-1].temps, duree);
    msg("CQT interpolation : creation of a matrix of size {} * {}.", onmax, cqtk.nfreqs);

    ArrayXXf A = ArrayXXf::Zero(onmax, cqtk.nfreqs);
    // Echantillons classés par indices de fréquence croissants
    for(auto i = 0; i < ns; i++)
    {
      auto &spl = sortie[i];
      int fr = spl.id_freq;
      if(fr >= cqtk.nfreqs)
        echec("ofreq invalid.");

      // Temps échantillon précédent (pour cette fréquence), et échantillon courant
      // (en nombre d'échantillons d'entrée)
      float t0 = lindex(fr), t1 = spl.temps;
      // Valeur précédente et courante
      float y0 = lvalue(fr), y1 = std::abs(spl.valeur);
      // Compléte la matrice de sortie
      while(offset(fr) * oprd <= t1)
      {
        float t = offset(fr) * oprd;
        float val;
        if(t1 == t0)
            val = y1;
        else
            val = y0 + ((t - t0) / (t1 - t0)) * (y1 - y0);
        A(offset(fr),fr) = val;
        offset(fr)++;
      }
      lindex(fr) = t1;
      lvalue(fr) = y1;
    }
     ArrayXf t = linspace(0,A.rows()/ofs,A.rows());
     // [t,f,A]
     return {t, cqtk.freqs, A};
  }

};

CQT::CQT()
{
  impl = std::make_shared<Impl>();
}

void CQT::configure(const CQTConfig &config)
{
  impl->configure(config);
}

void CQT::step(const ArrayXf &x)
{
  impl->step(x);
}

void CQT::affiche_noyaux()
{
  impl->cqtk.plot();
}

std::tuple<ArrayXf, ArrayXf, ArrayXXf> CQT::interpolation(float ofs)
{
  return impl->interpolation(ofs);
}



}

namespace tsd::tf {

  using namespace tsd::tf::cqt;

  ArrayXXf periodogramme_cqt(const ArrayXf &x, float fs, float fmin, float fmax, float γ, float ofs)
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
    auto [t, f, A] = cqt.interpolation(ofs);
    return A;
  }


}
