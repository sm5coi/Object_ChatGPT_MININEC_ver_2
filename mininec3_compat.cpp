#include "mininec3_compat.hpp"
#include "mininec3_state.hpp"
#include <iostream>

static int hits = 0;

void Mininec3Compat::computeObsDeltaR(int I)
{
    int c1 = st_.C[I][0];
    int c2 = st_.C[I][1];

    int I1 = std::abs(c1);
    int I2 = std::abs(c2);

    double F4 = ((c1 >= 0) ? +1.0 : -1.0) * st_.S[I1];
    double F5 = ((c2 >= 0) ? +1.0 : -1.0) * st_.S[I2];

    st_.T5 = F4 * st_.CA[I1] + F5 * st_.CA[I2];
    st_.T6 = F4 * st_.CB[I1] + F5 * st_.CB[I2];
    st_.T7 = F4 * st_.CG[I1] + F5 * st_.CG[I2];

    if (c1 == -c2)
        st_.T7 = st_.S[I1] * (st_.CG[I1] + st_.CG[I2]);
}

// BASIC 235–245
int Mininec3Compat::computeF8flag(int I, int J, int I1, int I2, int J1, int J2) const
{
    int F8 = 0;

    if (I1 != I2) return 0;

    constexpr double eps = 1e-12;

    if (std::abs(st_.CA[I1] + st_.CB[I1]) > eps)
    {
        if (st_.C[I][0] != st_.C[I][1]) return 0;
    }

    if (J1 != J2) return 0;

    if (std::abs(st_.CA[J1] + st_.CB[J1]) > eps)
    {
        if (st_.C[J][0] != st_.C[J][1]) return 0;
    }

    if (I1 == J1) F8 = 1;
    if (I == J)   F8 = 2;

    return F8;
}


void Mininec3Compat::kernel28(double T, int Kimg, int P4,
                              const Vec3& X2, const Vec3& V,
                              double A2, double I6,
                              double& T3, double& T4) const
{
    // BASIC 28–35: choose direction depending on K (+1 or -1)
    Vec3 X3;
    if (Kimg >= 0)
    {
        // X3 = X2 + T*(V - X2)
        X3 = X2 + (V - X2) * T;
    }
    else
    {
        // X3 = V + T*(X2 - V)
        X3 = V + (X2 - V) * T;
    }

    double D3 = X3.x*X3.x + X3.y*X3.y + X3.z*X3.z;

    double D = 0.0;

    // BASIC 38–40
    if (st_.A[P4] <= st_.SRM)
    {
        D = std::sqrt(D3);
    }
    else
    {
        D = std::sqrt(D3 + A2);
    }

    // double Dmin = (st_.SRM > 0.0) ? st_.SRM : 1e-6;
    // if (D < Dmin) D = Dmin;

    // BASIC 42–48 exact kernel branch (elliptic integral approx)
    // if (I6 != 0.0)
    // {
    //     // B = D3 / (D3 + 4*A2)
    //     double B = D3 / (D3 + 4.0 * A2);

    //     // These coefficients must exist in state: C0..C9 and P (pi)
    //     double W0 = st_.C0 + B * (st_.C1 + B * (st_.C2 + B * (st_.C3 + B * st_.C4)));
    //     double W1 = st_.C5 + B * (st_.C6 + B * (st_.C7 + B * (st_.C8 + B * st_.C9)));

    //     double V0 = (W0 - W1 * std::log(B)) * std::sqrt(1.0 - B);

    //     T3 += (V0 + std::log(D3 / (64.0 * A2)) / 2.0) / M_PI / st_.A[P4] - 1.0 / D;
    // }

    // BASIC 49–52
    double B1 = D * st_.W;
    T3 += std::cos(B1) / D;
    T4 -= std::sin(B1) / D;
}


void Mininec3Compat::psiCore102_28(const Vec3& X1,
                                   double P2, double P3,
                                   int P4, int Kimg,
                                   int FVS)
{
    // ------------------------------------------------------------
    // BASIC 112–134: build X2 and V relative to X1
    // ------------------------------------------------------------

    // IMPORTANT: BASIC uses:
    //   X1,Y1,Z1 = observation point (NOT multiplied by K)
    //   Z(P2), Z(P3) are multiplied by K before subtracting Z1
    //
    // We already pass X1 as un-imaged point from psiVector/psiScalar.
    // Now we fetch Ru/Rv WITHOUT image and apply K only to z.
    Vec3 Ru = pointFromP_Basic(P2, +1); // no image here
    Vec3 Rv = pointFromP_Basic(P3, +1); // no image here

    Ru.z *= (double)Kimg;
    Rv.z *= (double)Kimg;

    Vec3 X2 = Ru - X1;
    Vec3 V  = Rv - X1;

    // ------------------------------------------------------------
    // BASIC 135–139: D0 and D3 magnitudes
    // ------------------------------------------------------------
    double D0 = std::sqrt(X2.x*X2.x + X2.y*X2.y + X2.z*X2.z);
    double D3 = std::sqrt(V.x*V.x  + V.y*V.y  + V.z*V.z);

    // ------------------------------------------------------------
    // BASIC 141: A2 = A(P4)^2
    // ------------------------------------------------------------
    double A2 = st_.A[P4] * st_.A[P4];

    // ------------------------------------------------------------
    // BASIC 143: S4 = (P3 - P2)*S(P4)
    // ------------------------------------------------------------
    double S4 = (P3 - P2) * st_.S[P4];

    // ------------------------------------------------------------
    // BASIC 146–151 init
    // ------------------------------------------------------------
    double T1sum = 0.0;
    double T2sum = 0.0;

    double I6 = 0.0;   // BASIC I6! (exact kernel add-on)
    double F2 = 1.0;   // BASIC F2

    int L = 7;         // BASIC default quadrature order
    double Tcrit = (D0 + D3) / st_.S[P4];  // BASIC 151

    // ------------------------------------------------------------
    // BASIC 153–164: EXACT KERNEL CRITERIA
    // ------------------------------------------------------------
    // BASIC also checks "C$ = N" to disable exact kernel.
    // We assume enabled for now.
    bool exactKernelEnabled = true;

    bool usedExactKernel = false;

    if (Tcrit <= 1.1 && exactKernelEnabled)
    {
        // BASIC checks junction connectivity:
        //   J2(W%(I),1/2) compared to J2(W%(J),1/2)
        //
        // We use the scratch wires set in buildZ():
        //   st_.WpulseScratchObs = W%(I)
        //   st_.WpulseScratchSrc = W%(J)
        int wI = st_.WpulseScratchObs;
        int wJ = st_.WpulseScratchSrc;

        int a1 = st_.J2[wI][0];
        int a2 = st_.J2[wI][1];
        int b1 = st_.J2[wJ][0];
        int b2 = st_.J2[wJ][1];

        bool shareJunction =
            (a1 == b1) || (a1 == b2) || (a2 == b1) || (a2 == b2);

        if (shareJunction)
        {
            // BASIC 160–161: if small radius -> jump to thin-wire self special
            if (st_.A[P4] <= st_.SRM)
            {
                // BASIC:
                //   IF FVS = 1 THEN 91 ELSE GOTO 106
                // i.e. scalar uses 2*log, vector uses log.
                if (FVS == 1)
                {
                    // BASIC 91–93
                    st_.T1 = 2.0 * std::log(st_.S[P4] / st_.A[P4]);
                    st_.T2 = -st_.W * st_.S[P4];
                    return;
                }
                else
                {
                    // BASIC 106–108
                    st_.T1 = std::log(st_.S[P4] / st_.A[P4]);
                    st_.T2 = -st_.W * st_.S[P4] / 2.0;
                    return;
                }
            }

            // BASIC 162
            F2 = 2.0 * (P3 - P2);

            // BASIC 163:
            // I6! = (1 - LOG(S4 / F2 / 8 / A(P4))) / P / A(P4)
            // where P = pi
            double denom = (F2 * 8.0 * st_.A[P4]);
            if (denom <= 0.0) denom = 1e-30;

            double arg = S4 / denom;
            if (arg <= 0.0) arg = 1e-30;

            I6 = (1.0 - std::log(arg)) / (M_PI * st_.A[P4]);

            usedExactKernel = true;
        }
    }

    // ------------------------------------------------------------
    // BASIC 165–166: choose integration order (only if NOT exact kernel)
    // ------------------------------------------------------------
    if (!usedExactKernel)
    {
        if (Tcrit > 6.0)  L = 3;
        if (Tcrit > 10.0) L = 1;
    }

    // ------------------------------------------------------------
    // BASIC 167–178: Gaussian quadrature loop using packed Q()
    // ------------------------------------------------------------
    auto Q = [&](int idx1based) -> double {
        return st_.Qbasic[idx1based - 1];
    };

    int I5 = L + L; // BASIC 167: I5 = L + L

    // BASIC uses L as an index into Q and increments it inside the loop.
    // We'll mimic that exactly using idx = L.
    int idx = L;

    while (idx < I5)
    {
        double T3 = 0.0;
        double T4 = 0.0;

        // BASIC 170–173
        double T = (Q(idx) + 0.5) / F2;
        kernel28(T, Kimg, P4, X2, V, A2, I6, T3, T4);

        T = (0.5 - Q(idx)) / F2;
        kernel28(T, Kimg, P4, X2, V, A2, I6, T3, T4);

        // BASIC 174–176
        idx = idx + 1;
        T1sum += Q(idx) * T3;
        T2sum += Q(idx) * T4;

        // BASIC 177
        idx = idx + 1;
    }

    // ------------------------------------------------------------
    // BASIC 179–180: final results
    // ------------------------------------------------------------
    st_.T1 = S4 * (T1sum + I6);
    st_.T2 = S4 * T2sum;
}

static inline bool nearlyEqual(double a, double b, double eps = 1e-9)
{
    return std::abs(a - b) <= eps;
}

void Mininec3Compat::psiVector(double P1, double P2, double P3, int P4, int Kimg)
{
    // BASIC 102: FVS = 0
    int FVS = 0;

    if (P4 == 0)
    {
        st_.T1 = 0.0;
        st_.T2 = 0.0;
        return;
    }

    // -------------------------------
    // BASIC 103–108: thin-wire vector self special
    // IF K < 1 THEN 109
    // IF A(P4) >= SRM THEN 109
    // IF (I = J AND P3 = P2 + .5) THEN 106
    // -------------------------------
    if (Kimg >= 1)
    {
        if (st_.A[P4] < st_.SRM)
        {
            bool samePulse = (st_.curI == st_.curJ);

            // BASIC: P3 = P2 + 0.5
            if (samePulse && nearlyEqual(P3, P2 + 0.5))
            {
                st_.T1 = std::log(st_.S[P4] / st_.A[P4]);
                st_.T2 = -st_.W * st_.S[P4] / 2.0;

                // Debug counter (optional)
                st_.psiVectorHits++;

                std::cout << "psiVector self special HIT: curI=" << st_.curI
                          << " curJ=" << st_.curJ
                          << " P2=" << P2 << " P3=" << P3
                          << " P4=" << P4
                          << " A=" << st_.A[P4] << " SRM=" << st_.SRM
                          << "\n";

                return;
            }
        }
    }

    // -------------------------------
    // BASIC 109–112: X1 = X(P1),Y(P1),Z(P1) (NO K on Z here!)
    // -------------------------------
    Vec3 X1 = pointFromP_Basic(P1, +1); // IMPORTANT: K=+1 to avoid mirroring observation point

    // Jump to shared core (BASIC 113)
    psiCore102_28(X1, P2, P3, P4, Kimg, FVS);
}

void Mininec3Compat::psiScalar(double P1, double P2, double P3, int P4, int Kimg)
{
    // BASIC 87: FVS = 1
    int FVS = 1;

    if (P4 == 0)
    {
        st_.T1 = 0.0;
        st_.T2 = 0.0;
        return;
    }

    // -------------------------------
    // BASIC 88–93: thin-wire scalar self special
    // IF K < 1 THEN 94
    // IF A(P4) > SRM THEN 94
    // IF (P3 = P2 + 1 AND P1 = (P2 + P3)/2) THEN 91
    // -------------------------------
    if (Kimg >= 1)
    {
        if (st_.A[P4] <= st_.SRM)
        {
            bool cond1 = nearlyEqual(P3, P2 + 1.0);
            bool cond2 = nearlyEqual(P1, 0.5 * (P2 + P3));

            if (cond1 && cond2)
            {
                st_.T1 = 2.0 * std::log(st_.S[P4] / st_.A[P4]);
                st_.T2 = -st_.W * st_.S[P4];
                return;
            }
        }
    }

    // -------------------------------
    // BASIC 94–98: midpoint between INT(P1) and INT(P1)+1
    // X1 = (X(I4)+X(I5))/2  (NO K on X1.z here!)
    // -------------------------------
    int I4 = (int)std::floor(P1 + 1e-12);
    int I5 = I4 + 1;

    Vec3 X1;
    X1.x = 0.5 * (st_.BX[I4] + st_.BX[I5]);
    X1.y = 0.5 * (st_.BY[I4] + st_.BY[I5]);
    X1.z = 0.5 * (st_.BZ[I4] + st_.BZ[I5]); // IMPORTANT: no K here

    // Jump into shared core (BASIC 113)
    psiCore102_28(X1, P2, P3, P4, Kimg, FVS);
}


static inline double CA0(const Mininec3State& st, int w) { return (w > 0) ? st.CA[w] : 0.0; }
static inline double CB0(const Mininec3State& st, int w) { return (w > 0) ? st.CB[w] : 0.0; }
static inline double CG0(const Mininec3State& st, int w) { return (w > 0) ? st.CG[w] : 0.0; }

static inline double Ssafe(const Mininec3State& st, int w)
{
    // w==0 means "no half-pulse" (BASIC C%=0). Must never divide by S(0).
    if (w <= 0) return 1.0; // denominator dummy, caller should skip contribution anyway
    double s = st.S[w];
    return (s > 1e-12) ? s : 1.0;
}

void Mininec3Compat::computeZijForImage(int I, int J, int K,
                                        int /*I1*/, int /*I2*/, int J1, int J2,
                                        double F4, double F5,
                                        double F6, double F7,
                                        int F8)
{
    // (valfritt) debug
    // if (J1==0 || J2==0)
    //     std::cout << "DEBUG: J1="<<J1<<" J2="<<J2<<" for I="<<I<<" J="<<J<<"\n";

    std::cout << "DEBUG S: S[J1]=" << st_.S[J1]
              << " S[J2]=" << st_.S[J2]
              << " Ssafe(J1)=" << Ssafe(st_, J1)
              << " Ssafe(J2)=" << Ssafe(st_, J2)
              << "\n";


    int pI = I + 1;
    int pJ = J + 1;

    double P1 = 2.0 * st_.Wpulse[I] + pI - 1;
    double P2 = 2.0 * st_.Wpulse[J] + pJ - 1;
    double P3 = P2 + 0.5;

    int P4 = J2;

    psiVector(P1, P2, P3, P4, K);
    double U1 = F5 * st_.T1;
    double U2 = F5 * st_.T2;

    P3 = P2;
    P2 = P2 - 0.5;
    P4 = J1;

    if (F8 < 2)
        psiVector(P1, P2, P3, P4, K);

    double V1 = F4 * st_.T1;
    double V2 = F4 * st_.T2;

    // --- vector potential contribution ---
    double X3 = U1 * CA0(st_, J2) + V1 * CA0(st_, J1);
    double Y3 = U1 * CB0(st_, J2) + V1 * CB0(st_, J1);
    double Z3 = (F7 * U1 * CG0(st_, J2) + F6 * V1 * CG0(st_, J1)) * (double)K;

    double D1 = st_.W2 * (X3 * st_.T5 + Y3 * st_.T6 + Z3 * st_.T7);

    X3 = U2 * CA0(st_, J2) + V2 * CA0(st_, J1);
    Y3 = U2 * CB0(st_, J2) + V2 * CB0(st_, J1);
    Z3 = (F7 * U2 * CG0(st_, J2) + F6 * V2 * CG0(st_, J1)) * (double)K;

    double D2 = st_.W2 * (X3 * st_.T5 + Y3 * st_.T6 + Z3 * st_.T7);

    // --- scalar part ---
    P1 = P1 + 0.5;
    if (F8 == 2) P1 = P1 - 1.0;

    P2 = P3;
    P3 = P3 + 1.0;
    P4 = J2;

    double U5 = 0.0, U6 = 0.0;

    if (F8 == 1)
    {
        U5 = F5 * U1 + st_.T1;
        U6 = F5 * U2 + st_.T2;
    }
    else
    {
        psiScalar(P1, P2, P3, P4, K);
        U5 = st_.T1;
        U6 = st_.T2;

        if (F8 >= 2)
        {
            U1 = (2.0 * st_.T1 - 4.0 * U1 * F5) / Ssafe(st_, J1);
            U2 = (2.0 * st_.T2 - 4.0 * U2 * F5) / Ssafe(st_, J1);

            st_.ZR[I][J] += (double)K * (D1 + U1);
            st_.ZI[I][J] += (double)K * (D2 + U2);
            return;
        }
    }

    P1 = P1 - 1.0;
    psiScalar(P1, P2, P3, P4, K);

    // U1 = (st_.T1 - U5) / Ssafe(st_, J2);
    // U2 = (st_.T2 - U6) / Ssafe(st_, J2);

    if (J2 != 0)
    {
        U1 = (st_.T1 - U5) / st_.S[J2];
        U2 = (st_.T2 - U6) / st_.S[J2];
    }
    else
    {
        U1 = 0.0;
        U2 = 0.0;
    }

    P1 = P1 + 1.0;
    P3 = P2;
    P2 = P2 - 1.0;
    P4 = J1;

    psiScalar(P1, P2, P3, P4, K);
    double U3 = st_.T1;
    double U4 = st_.T2;

    if (F8 >= 1)
    {
        st_.T1 = U5;
        st_.T2 = U6;
    }
    else
    {
        P1 = P1 - 1.0;
        psiScalar(P1, P2, P3, P4, K);
    }

    // U1 = U1 + (U3 - st_.T1) / Ssafe(st_, J1);
    // U2 = U2 + (U4 - st_.T2) / Ssafe(st_, J1);

    if (J1 != 0)
    {
        U1 = U1 + (U3 - st_.T1) / st_.S[J1];
        U2 = U2 + (U4 - st_.T2) / st_.S[J1];
    }

    st_.ZR[I][J] += (double)K * (D1 + U1);
    st_.ZI[I][J] += (double)K * (D2 + U2);
}

void Mininec3Compat::applySymmetryAndToeplitzCopy(int I, int J, int F8, int J2)
{
    // BASIC 317
    if (J < I) return;

    // BASIC 318
    if (F8 == 0) return;

    // BASIC 319–320 reciprocity
    st_.ZR[J][I] = st_.ZR[I][J];
    st_.ZI[J][I] = st_.ZI[I][J];

    // BASIC 322–323
    int P1 = J + 1;
    if (P1 >= st_.N) return;

    // BASIC 324
    if (st_.C[P1][0] != st_.C[P1][1]) return;

    // BASIC 325–327
    if (st_.C[P1][1] != st_.C[J][1])
    {
        if (st_.C[P1][1] != -st_.C[J][1]) return;
        if (std::abs(st_.CA[J2] + st_.CB[J2]) > 1e-12) return;
    }

    // BASIC 328–329
    int P2 = I + 1;
    if (P2 >= st_.N) return;

    // BASIC 330–331
    st_.ZR[P2][P1] = st_.ZR[I][J];
    st_.ZI[P2][P1] = st_.ZI[I][J];
}

void Mininec3Compat::buildZ()
{
    for (int I = 0; I < st_.N; ++I)
    {
        computeObsDeltaR(I);

        int cI1 = st_.C[I][0];
        int cI2 = st_.C[I][1];
        int I1 = std::abs(cI1);
        int I2 = std::abs(cI2);

        for (int J = 0; J < st_.N; ++J)
        {
            st_.curI = I + 1;
            st_.curJ = J + 1;

            int cJ1 = st_.C[J][0];
            int cJ2 = st_.C[J][1];

            int J1 = std::abs(cJ1);
            int J2 = std::abs(cJ2);

            double F4 = (cJ1 >= 0) ? +1.0 : -1.0;
            double F5 = (cJ2 >= 0) ? +1.0 : -1.0;

            double F6 = 1.0;
            double F7 = 1.0;

            // BASIC: image loop FOR K = 1 TO G STEP -2
            for (int K : {+1, -1})
            {
                // BASIC 231–234
                if (cJ1 == -cJ2)
                {
                    if (K < 0) break;
                    F6 = F4;
                    F7 = F5;
                }

                // BASIC 235–245
                int F8 = 0;
                if (K >= 0)
                    F8 = computeF8flag(I, J, I1, I2, J1, J2);

                // BASIC 246: reuse
                if (st_.ZR[I][J] != 0.0)
                {
                    applySymmetryAndToeplitzCopy(I, J, F8, J2);
                    continue;
                }

                // BASIC scratch wires used by exact kernel logic
                st_.WpulseScratchObs = st_.Wpulse[I];
                st_.WpulseScratchSrc = st_.Wpulse[J];

                // BASIC 248–316
                computeZijForImage(I, J, K, I1, I2, J1, J2, F4, F5, F6, F7, F8);

                // BASIC 317–331
                applySymmetryAndToeplitzCopy(I, J, F8, J2);
            }
        }
    }
}


Vec3 Mininec3Compat::pointFromP_Basic(double P, int /*Kimg*/) const
{
    auto bp = [&](int idx) -> Vec3 {
        if (idx < 0 || idx >= (int)st_.BX.size())
            throw std::out_of_range("Breakpoint index out of range");
        return Vec3(st_.BX[idx], st_.BY[idx], st_.BZ[idx]);
    };

    double iFloor = std::floor(P + 1e-12);
    double frac   = P - iFloor;
    int i = (int)iFloor;

    if (std::abs(frac) < 1e-9)
        return bp(i);

    if (std::abs(frac - 0.5) < 1e-9)
        return (bp(i) + bp(i+1)) * 0.5;

    // fallback
    return bp(i);
}


