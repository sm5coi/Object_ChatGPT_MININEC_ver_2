#include "frequencyInput.hpp"



void frequencyInput::set_fq(double fq)
{
    if (Fmhz == 0.0) Fmhz = 299.8; // BASIC default

    double lambda = 299.8 / Fmhz;  // meters
    st_.SRM = 0.0001 * lambda;

    double k = 2.0 * M_PI / lambda;
    st_.W  = k;
    st_.W2 = (k * k) / 2.0;

}
