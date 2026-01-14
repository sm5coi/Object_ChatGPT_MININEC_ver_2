#ifndef MININEC3_COMPAT_HPP
#define MININEC3_COMPAT_HPP

#include <vector>
#include <array>
#include <cmath>
#include "vec3.hpp"

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

    // C%(I,1..2) in BASIC (signed segment IDs)
    std::vector<std::array<int,2>> C;  // size N

    // Per-segment arrays (index = segment id)
    std::vector<double> S, A, CA, CB, CG;

    // Wire index per pulse (BASIC W%(I))
    std::vector<int> Wpulse; // size N

    // Impedance matrix parts
    std::vector<std::vector<double>> ZR, ZI;

    // Scratch (matches BASIC globals)
    double T1=0.0, T2=0.0;
    double T5=0.0, T6=0.0, T7=0.0;

    // Node coordinates (BASIC: X(), Y(), Z())
    std::vector<double> X, Y, Z;

    // Wire node chains:
    // wireNodes[w][i] = nodeId at integer position i along wire w
    std::vector<std::vector<int>> wireNodes;

    // Base offset in "P-space" for each wire.
    // BASIC: P = 2*W%(I) + I - 1  => 2*W%(I) acts like a base.
    std::vector<int> wireBaseP;

    // BASIC-style breakpoint coordinates X(),Y(),Z()
    std::vector<double> BX, BY, BZ;

    // Wire endpoints storage (BASIC E(),L(),M())
    std::vector<double> E, L, M;

    // For each wire: pulse range N(I,1), N(I,2) (BASIC N(,))
    std::vector<std::array<int,2>> wirePulseRange;

    // Connections (BASIC J2(I,1..2) and J1(I))
    std::vector<int> J1;
    std::vector<std::array<int,2>> J2;

    // MININEC ground flag (we use G=1 now)
    int Gmode = 1;
    int NW = 0; // number of wires

    std::array<double, 14> Qbasic{};
};

class Mininec3Compat
{
public:
    explicit Mininec3Compat(Mininec3State& st) : st_(st) {}

    void buildZ();

    void psiGaussCoreCompat(double P1, double P2, double P3, int P4, int Kimg);

private:
    Mininec3State& st_;

    void computeObsDeltaR(int I);

    int computeF8flag(int I, int J, int I1, int I2, int J1, int J2) const;

    void psiVector(double P1, double P2, double P3, int P4, int Kimg); // GOSUB 102
    void psiScalar(double P1, double P2, double P3, int P4, int Kimg); // GOSUB 87


    void computeZijForImage(int I, int J, int K,
                            int I1, int I2, int J1, int J2,
                            double F4, double F5,
                            double F6, double F7,
                            int F8);

    void applySymmetryAndToeplitzCopy(int I, int J, int F8, int J2);

    Vec3 pointFromP(double P, int wireIndex, int Kimg) const;

    Vec3 pointFromGlobalP(double Pglobal, int wireIndex, int Kimg) const;

    Vec3 pointFromP_Basic(double P, int Kimg) const;

    void psiCore102_28(const Vec3& X1,
                                       double P2, double P3,
                                       int P4, int Kimg,
                                       int FVS);

    void kernel28(double T, int Kimg, int P4,
                                  const Vec3& X2, const Vec3& V,
                                  double A2, double I6,
                                  double& T3, double& T4) const;


};


#endif // MININEC3_COMPAT_HPP
