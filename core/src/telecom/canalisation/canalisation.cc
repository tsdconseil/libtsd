#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"

using namespace tsd::filtrage;

namespace tsd::telecom {








// Creation of a polyphase decimation filter structure
//
// Calling Sequence
// H = polyphase_filter(h,m)
//
//  Parameters
//  h: FIR filter coefficients (1d vector)
//  m: number of polyphase branches e.g. decimation ratio
//  H: resulting polyphase matrix form
//
// Description
// Creation of the polyphase matrix. The filter is decomposed into m
// phases sub-filters (m rows of the H matrix). Each sub-filter can
// process a decimated version of the input signal, with a different phase.
/*ArrayXXf filtre_polyphase(IArrayXf h, unsigned int m)
{
  ArrayXf h2 = h;
  soit ntaps = h.rows();
  soit reste = ntaps % m;
  si(reste != 0)
  {
    //printf("Ntaps not multiple of nchn.\nPadding with %d zeros: ntaps %d", m - reste, ntaps);
    h2 = vconcat(h2, ArrayXf::Zero(m-reste));
    ntaps = h2.rows();
  }
  tsd_assert((ntaps % m) == 0);
  // en ligne : les filtres
  // en colonne : les canaux
  retourne Eigen::Map<ArrayXXf>(h2.data(), m, ntaps/m);
  //H = flipdim(H,1);
}*/





// Channelization: frequency multiplexing of m input signals into a single signal (but with bandwidth multiplied by m)
//
// Calling Sequence
// y = channelize(X)
//
// Parameters
// X: input matrix, size [n x m], with n: n samples / channel, and m: number of channels. Each column of X is a different signal. The number of columns (m) is the number of channels.
// y: output vector, size [nm x 1]
//
// Description
// Merge m different signals into a single vector, by frequency multiplexing.
// The signals are shifted at the following normalized frequencies: 0, 1/m, 2/m, ..., (m-1)/m.
// Note: this could be done more effectively using a modulated filter bank (Harris method, reciprocal algorithm of the unchannelize function).
// <programlisting>
// fs = 1e3;
// mod = mod_init("bpsk", fs=1e3,fi=0,fsymb=50);
// [mod,x] = mod_process(mod,prbs(1000));
// nchn = 8;
// // In this example, just duplicate the same channel 8 times
// X = repmat(x,1,nchn);
// y = channelize(X); // output sample rate is 8 times higher
// clf(); plot_psd(y,8*fs,'b');
// </programlisting>
// <imageobject><imagedata fileref="ex_channelize.png" format="PNG"/></imageobject>
Veccf canalisation(const Tabcf X, const Vecf &h)
{

  // A FAIRE : remplacer par un filtre au fil de l'eau

  soit n = X.rows(); // Nb éch / canal
  soit m = X.cols(); // Nb canaux

//    y = zeros(n*m,1);
//
//    pour i = 1:m
//        nu = (i-1) / m; // normalized frequency
//        xu = intdec(real(X(:,i)), m) + %i * intdec(imag(X(:,i)), m);
//        xm = xu .* exp(2*%pi*%i*nu*(0:n*m-1)');
//        y = y + xm;
//    end
//


  // IDFT suivant les colonnes (dim n°1)
  Tabcf Y(n, m);
  pour(auto k = 0; k < m; k++)
    Y.col(k) = tsd::fourier::ifft(X.col(k));

    //soit Y = tsd::fourier::ifft(X);

  // Y : chaque ligne = une fréquence donnée

  // (3) polyphase partition of the filter
  soit H = forme_polyphase(h, m); // Transférer ceci dans le module de filtrage
  // H : m lignes, len(h)/m colonnes

    // (2) Sortie du filtre polyphase
  soit n2 = n;
  soit n3 = n2 + (length(h)-1)/m;
  soit XF = Tabcf::zeros(m, n3);
  pour(auto i = 0; i < m; i++)
    XF.row(i) = tsd::filtrage::convol<cfloat,float>(H.row(i),Y.row(i));

  // (1) polyphase partition of the signal
  // On suppose qu'on a un bloc de nchn échantillons
  // A partir du signal d'entrée de n échantillons,
  // on forme nchn signaux de n / nchn échantillons
  retourne iforme_polyphase<cfloat>(XF);
  // Chaque ligne de X peut être traitée séparément
}


// Apply a FIR filter and decimate the output, using an efficient polyphase structure
//
// Calling Sequence
// y = polyphase_decimation(x,h,R)
//
//  Parameters
//  x: input signal (1d vector)
//  h: FIR filter coefficients (1d vector)
//  R: number of polyphase branches e.g. decimation ratio
//  y: filtered and decimated output
//
//  See also
//   unchannelize
//   channelize
//
// Bibliography
//   F.J. HARRIS, Digital Receivers and Transmitters Using Polyphase Filter Banks pour Wireless Communications, 2003


// TODO : à supprimer et remplacer par filtre_rif_decim
/*ArrayXf decimation_polyphase(IArrayXf x, IArrayXf h, unsigned int R)
{
  ArrayXXf H = forme_polyphase(h, R);
  soit ntaps = length(h);
  soit X = forme_polyphase(x, R);
  soit n2 = X.cols();
  soit n3 = ceil(n2 + ((float)ntaps)/R - 1);
  ArrayXXf XF = ArrayXXf::Zero(R, n3);
  pour(auto i = 0u; i < R; i++)
    XF.row(i) = tsd::filtrage::convol<float,float>(H.row(i), X.row(i));
  ArrayXf y = ArrayXf::Zero(n3);
  pour(auto i = 0; i < n3; i++)
      y(i) = XF.col(i).sum();
  retourne y;
}*/


}
