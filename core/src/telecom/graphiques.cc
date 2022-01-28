#include "tsd/telecom.hpp"
#include "tsd/fourier.hpp"
#include "tsd/figure.hpp"

using namespace tsd::vue;


namespace tsd::telecom {



void plot_eye(Figure &f, const ArrayXf &x, float T)
{


# if 0
  auto trig2 = [](const ArrayXf &x, unsigned int T) -> ArrayXf
  {
    // Reference signal with pulse @ t = 0
    int n = x.rows();
    ArrayXf A = ArrayXf::Zero(T);
    A(0) = 1;
    ArrayXf ref = repmat(A, floor(n / T), 1);

    // Signal with pulses at the transitions
    ArrayXf y = diff(x).square();

    float score;
    int delais = tsd::fourier::estimation_delais_entier(ref, y, score);
    delais = 1 + (delais % T);
    ArrayXf tt = delais + linspace(T/2,n,n/T);    //T/2:T:length(x);
    return tt;
  };

  int n = length(x);
  ArrayXf tt = trig2(x,T);
  //tt(find(tt < 1)) = [];
  ArrayXf t = ArrayXf::Zero(n);
  ArrayXf nt = ArrayXf::Zero(2*T*length(tt));
  ArrayXf nx = nt;

  for(auto i = 0u; i < tt.rows(); i++)
  {
    ti = tt(i):(tt(i+1)-1);
    t(ti) = ti-tt(i);

    itrv = 1+(i-1)*2*T:i*2*T;

    nt(itrv) = (0:2*T-1);
    if(tt(i)+2*T-1 <= n)
      nx(itrv) = x(tt(i):(tt(i)+2*T-1));
    else
    {
      nt(1+(i-1)*2*T:$) = [];
      x(1+(i-1)*2*T:$) = [];
      break;
    }
  }

  //plot(nt, nx, 'b.');

  // Accumulation
  auto sx = (int) (2*T), sy = (int) floor(sx * 2 / 3);
  auto img = ArrayXXf::Zero(sx,sy);
  nt = 1 + uint32_t(floor(nt*(sx-1)/max(nt)));
  y  = 1 + uint32_t(floor((sy-1) .* (nx - min(nx)) ./ (max(nx) - min(nx))));

  n = length(nt);
  // maintenant on sait que sur des intervalles kT...k+1 T,
  // nt est une simple rampe.
  for(auto i = 0u; i < n; i++)
      // nt(i) = x, ny(i) = y
      img(nt(i),y(i)) = img(nt(i),y(i)) + 1;

  // xset("colormap",graycolormap(512));
  // Sgrayplot(1:sx,1:sy,-img+max(img));
  // xtitle("Eye diagram");
# endif
}


}
