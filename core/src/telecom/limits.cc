#include "tsd/telecom.hpp"



namespace tsd::telecom {

float capacite_canal_awgn(float snr, float B)
{
  return B * log2(1 + snr);
}

}







