#include "tsd/tsd-all.hpp"


namespace tsd::fourier
{

  // D'après wikipédia
  Vecf fwht(const Vecf &x)
  {
    soit n = x.rows(), h = 1;
    Vecf y = x;

    tantque(h < n)
    {
      pour(auto i = 0; i < n; i += h * 2)
      {
        pour(auto j = i; j < i + h; j++)
        {
          soit x1 = y(j);
          soit y1 = y(j + h);
          y(j)     = x1 + y1;
          y(j + h) = x1 - y1;
        }
      }
      y /= 2;
      h *= 2;
    }

    retourne y;
  }

  Tabf hadamard_matrice(entier n)
  {
    si(n == 1)
     retourne Tabf::ones(1,1);

    soit M = Tabf::zeros(n, n);
    soit Mp = hadamard_matrice(n/2);

    M.block(0, n/2, 0, n/2)     = Mp;
    M.block(n/2, n/2, 0, n/2)   = Mp;
    M.block(0, n/2, n/2, n/2)   = Mp;
    M.block(n/2, n/2, n/2, n/2) = -Mp;

    retourne M;
  }

}


