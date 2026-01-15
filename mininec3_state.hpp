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

    int curI = 0;
    int curJ = 0;

    double W = 0.0;
    double W2 = 0.0;
    double SRM = 0.0;

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

    std::array<double, 14> Qbasic{};
};

#endif // MININEC3_STATE_HPP
