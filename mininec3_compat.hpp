#ifndef MININEC3_COMPAT_HPP
#define MININEC3_COMPAT_HPP

#include "mininec3_state.hpp"
#include "vec3.hpp"
#include <cmath>

class Mininec3Compat
{
public:
    explicit Mininec3Compat(Mininec3State& st) : st_(st) {}

    void buildZ();

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
