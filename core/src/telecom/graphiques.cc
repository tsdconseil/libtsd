#include "tsd/telecom.hpp"
#include "tsd/tsd-all.hpp"

namespace tsd::telecom {

void plot_eye(Figure &f, const Vecf &x, float T)
{
  soit trig2 = [](const Vecf &x, entier T) -> Vecf
  {
    // Reference signal with pulse @ t = 0
    soit n  = x.rows();
    soit A   = sigimp(T);
    soit ref = repmat(A, floor(n / T), 1);

    // Signal with pulses at the transitions
    soit y = square(diff(x));

    float score;
    entier delais = tsd::fourier::estimation_delais_entier(ref, y, score);
    delais = 1 + (delais % T);
    retourne  delais + linspace(T/2,n,n/T);    //T/2:T:length(x);
  };

  soit n = x.rows();
  soit tt = trig2(x,T);
  //tt(find(tt < 1)) = [];
  soit t  = Vecf::zeros(n),
       nt = Vecf::zeros(2*T*tt.rows()),
       nx = nt;

  pour(auto i = 0u; i < tt.rows(); i++)
  {
    ti = tt(i):(tt(i+1)-1);
    t(ti) = ti-tt(i);

    itrv = 1+(i-1)*2*T:i*2*T;

    nt(itrv) = (0:2*T-1);
    si(tt(i)+2*T-1 <= n)
      nx(itrv) = x(tt(i):(tt(i)+2*T-1));
    sinon
    {
      nt(1+(i-1)*2*T:$) = [];
      x(1+(i-1)*2*T:$) = [];
      break;
    }
  }

  //plot(nt, nx, 'b.');

  // Accumulation
  soit sx = (entier) (2*T), sy = (entier) floor(sx * 2 / 3);
  soit img = Tabf::zeros(sx,sy);
  nt = 1 + floor(nt*(sx-1)/max(nt));
  y  = 1 + floor((sy-1) .* (nx - min(nx)) ./ (max(nx) - min(nx)));

  n = length(nt);
  // maintenant on sait que sur des intervalles kT...k+1 T,
  // nt est une simple rampe.
  pour(auto i = 0u; i < n; i++)
      // nt(i) = x, ny(i) = y
      img(nt(i),y(i)) = img(nt(i),y(i)) + 1;

  // xset("colormap",graycolormap(512));
  // Sgrayplot(1:sx,1:sy,-img+max(img));
  // xtitle("Eye diagram");
}


}
