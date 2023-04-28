# libtsd
Libtsd provides functions in C++ for digital signal processing, with the following design goals:
  - Ease of use: introducing vector-type classes (Vecf, Veci, ...) for easy vector computing (Matlab / Numpy like),
  - A clear distinction between API (header files) and implementation (source files). 
  - A multilanguage API (english using dsp/dsp.hpp, french using tsd/tsd.hpp)

# Simple example (english API)

```cpp
#include "dsp/dsp-all.hpp"

int main()
{
  {
    // Example 1: simple plot
    Figure f;
    let x = linspace(0, 8 * π, 100);
    f.plot(x, cos(x));
    f.save("./cosinus.png");
  }
  {
    // Example 2: filter design
    let h = design_fir_wnd(31, "lp", 0.25);
    
    Figure f;
    f.plot(h, "|o", "impulse response");
    f.save("./lp-filter.png");
    
    let n = 500;
    let x = sigcos(0.01, n) + 0.1 * randn(n);
    let y = filter(h, x);
    
    f.clear();
    f.plot(x, "b-", "Noisy signal");
    f.plot(y, "r-", "Filtered signal");
    f.save("./filtering.png");
  }
  return 0;
}

```

# Simple example (french API)

```cpp
#include "tsd/tsd-all.hpp"

int main()
{
  {
    // Exemple 1 : tracé simple
    Figure f;
    soit x = linspace(0, 8 * π, 100);
    f.plot(x, cos(x));
    f.enregistrer("./cosinus.png");
  }
  {
    // Exemple 2 : conception de filtre
    soit h = design_rif_fen(31, "pb", 0.25);
    
    Figure f;
    f.plot(h, "|o", "réponse impulsionnelle");
    f.enregistrer("./filtre.png");
    
    soit n = 500;
    soit x = sigcos(0.01, n) + 0.1 * randn(n);
    soit y = filtrer(h, x);
    
    f.clear();
    f.plot(x, "b-", "Signal bruité");
    f.plot(y, "r-", "Signal filtré");
    f.enregistrer("./filtrage.png");
  }
  return 0;
}
```


# Resources

 - English documentation: https://tsdconseil.github.io/libtsd/en
 - Documentation en français : https://tsdconseil.github.io/libtsd/fr
 - Graphical demonstrator: (will soon be provided)

# License
Libtsd is available under two distinct licensing options:
 - <b>Free Open Source License:</b> The GNU GPL version of libtsd is made available via Open Source licensing. You are free to redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. This license does not include technical support.
 - <b>Commercial License:</b> Non-exclusive commercial licensing rights can be granted. This license include technical support (by email).

## What type of license is right for me?

If you are building an Open Source software, then the GPL license is best suited for your purpose.
Otherwise, e.g. if you are building a proprietary software, then you should purchase a commercial license.
Another reason to purchase a license is that a one year technical support is included with it.

For more information about contracting a commercial license, please  e-mail contact.AT.tsdconseil.fr

## Contributions
As this project is dual licensed (GPL V3 for open source projects, or commercial license for proprietary project), I need that any contributor accept a copyright assignment (see https://github.com/tsdconseil/libtsd/wiki/Copyright-assignment), otherwise I could not be able to grant proprietary licenses.  
