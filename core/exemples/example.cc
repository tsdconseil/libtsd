#include <dsp/dsp-all.hpp>

int main()
{
  Figure f;

  let x = linspace(0, 8*π, 100);

  f.plot(x, cos(x));

  f.save("./cosinus.png");

  return 0;
}
