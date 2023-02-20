#include "tsd/telecom.hpp"



namespace tsd::telecom {

float capacite_canal_awgn(float snr, float B)
{
  retourne B * log2(1 + snr);
}

}







