#ifndef MININEC3_STATE_HPP
#define MININEC3_STATE_HPP

#include <vector>
#include <array>

struct Mininec3State
{
    int N = 0;     // number of pulses
    int G = 1;     // image loop flag (usually 1)

    int WpulseScratchObs = 0;
    int WpulseScratchSrc = 0;

    int curI = 0;   // nuvarande observationspuls
    int curJ = 0;   // nuvarande källpuls

    double fq = 0.0;            // Frequency in Hz
    double lambda = 0.0;        // 299.8/fq (MHz)
    double Mw = 0.0;            // 1140 REM -- M = 1 / (4 * PI * OMEGA * EPSILON)
    double W = 0.0;             // k = (2*pi)/lambda
    double W2 = 0.0;            // (k^2)/2
    double SRM = 0.0;           // 0.0001*lambda
    double c0 = 299792458.0;    // speed of light

    std::vector<std::array<int,2>> C;

    std::vector<double> S, A, CA, CB, CG;
    std::vector<int> Wpulse;

    std::vector<std::vector<double>> ZR, ZI;

    double T1=0.0, T2=0.0;
    double T5=0.0, T6=0.0, T7=0.0;

    std::vector<double> X, Y, Z;

    std::vector<std::vector<int>> wireNodes;
    std::vector<int> wireBaseP;

    std::vector<double> BX, BY, BZ;

    std::vector<double> E, L, M;

    std::vector<std::array<int,2>> wirePulseRange;

    std::vector<int> J1;
    std::vector<std::array<int,2>> J2;

    int Gmode = 1;
    int NW = 0;

    //std::array<double, 14> Qbasic{};

    std::array<double, 14> Qbasic{
        0.288675135, 0.5,
        0.430568156, 0.173927423,
        0.169990522, 0.326072577,
        0.480144928, 0.050614268,
        0.398333239, 0.111190517,
        0.262766205, 0.156853323,
        0.091717321, 0.181341892
    };
};

#endif // MININEC3_STATE_HPP
