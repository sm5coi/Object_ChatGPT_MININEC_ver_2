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
    // BASIC 235-245 (only evaluated for K>0 in caller)

    // Default
    int F8 = 0;

    // IF I1 <> I2 THEN 246
    if (I1 != I2) return 0;

    // IF (CA(I1) + CB(I1)) = 0 THEN 241
    if ((st_.CA[I1] + st_.CB[I1]) != 0.0)
    {
        // IF C%(I,1) <> C%(I,2) THEN 246
        if (st_.C[I][0] != st_.C[I][1]) return 0;
    }

    // IF J1 <> J2 THEN 246
    if (J1 != J2) return 0;

    // IF (CA(J1) + CB(J1)) = 0 THEN 244
    if ((st_.CA[J1] + st_.CB[J1]) != 0.0)
    {
        // IF C%(J,1) <> C%(J,2) THEN 246
        if (st_.C[J][0] != st_.C[J][1]) return 0;
    }

    // IF I1 = J1 THEN F8 = 1
    if (I1 == J1) F8 = 1;

    // IF I = J THEN F8 = 2
    if (I == J) F8 = 2;

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
    // BASIC 112–134: build X2 and V relative to X1
    Vec3 Ru = pointFromP_Basic(P2, +1); // ospeglad punkt
    Vec3 Rv = pointFromP_Basic(P3, +1); // ospeglad punkt

    // BASIC: Z2 = K*Z(P2) - Z1, V3 = K*Z(P3) - Z1
    Ru.z *= (double)Kimg;
    Rv.z *= (double)Kimg;

    Vec3 X2 = Ru - X1;
    Vec3 V  = Rv - X1;


    // BASIC 135–139
    double D0 = std::sqrt(X2.x*X2.x + X2.y*X2.y + X2.z*X2.z);
    double D3 = std::sqrt(V.x*V.x  + V.y*V.y  + V.z*V.z);

    // BASIC 141
    double A2 = st_.A[P4] * st_.A[P4];

    // BASIC 143–151
    double S4 = (P3 - P2) * st_.S[P4];
    double F2 = 1.0;
    int L = 7;

    double Tcrit = (D0 + D3) / st_.S[P4];

    // BASIC 153–166: choose integration order
    // (NOTE: BASIC uses L=7 initially, then maybe L=3 or L=1)
    if (Tcrit > 6.0)  L = 3;
    if (Tcrit > 10.0) L = 1;

    // BASIC 148,163
    double I6 = 0.0;

    // BASIC 153–164: exact-kernel branch (we can add later)
    // It requires:
    // - C$ != "N"
    // - connectivity compare using J2(W%(I),1/2)
    // - and radius/SRM tests + FVS handling
    //
    // For now: skip (I6 stays 0)

    // BASIC 167–180
    double T1sum = 0.0;
    double T2sum = 0.0;

    // BASIC uses packed Q() array.
    // We'll store it 0-based in st_.Qbasic[0..13].
    auto Q = [&](int idx1based) -> double {
        return st_.Qbasic[idx1based - 1];
    };

    int I5 = L + L; // 2*L
    int idx = L;    // BASIC L is reused as index into Q()

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

        // BASIC 177–178
        idx = idx + 1;
    }

    // BASIC 179–180
    st_.T1 = S4 * (T1sum + I6);
    st_.T2 = S4 * T2sum;

    (void)FVS; // används först när vi lägger in exact-kernel

    if (Kimg < 0) {
        std::cout << "K=-1: X1.z=" << X1.z
                  << " Ru.z(beforeK?)=" << (Ru.z/(double)Kimg)
                  << " Ru.z(after)=" << Ru.z << "\n";
    }
}


void Mininec3Compat::psiVector(double P1, double P2, double P3, int P4, int Kimg)
{
    int FVS = 0;

    if (P4 == 0)
    {
        st_.T1 = 0.0;
        st_.T2 = 0.0;
        return;
    }


    // BASIC 103–108 special case
    if (Kimg >= 1)
    {
        if (st_.A[P4] < st_.SRM)
        {
            auto isHalfStep = [&](double a, double b) {
                return std::abs((a - b) - 0.5) < 1e-6;
            };

            if (st_.curI == st_.curJ && isHalfStep(P3, P2))
            {
                ++hits;
                std::cout << "hits = " << hits << std::endl;

                st_.T1 = std::log(st_.S[P4] / st_.A[P4]);
                st_.T2 = -st_.W * st_.S[P4] / 2.0;
                return;
            }
        }
    }

    // BASIC 109–111: X1 = X(P1) (integer point)
    Vec3 X1 = pointFromP_Basic(P1, +1);

    // jump to core (BASIC 113)
    psiCore102_28(X1, P2, P3, P4, Kimg, FVS);
}


void Mininec3Compat::psiScalar(double P1, double P2, double P3, int P4, int Kimg)
{
    // BASIC 87
    int FVS = 1;
    (void)FVS; // används i 102 för exact-kernel val, vi tar in senare

    if (P4 == 0)
    {
        st_.T1 = 0.0;
        st_.T2 = 0.0;
        return;
    }

    // BASIC 88–93: thin-wire scalar self special
    if (Kimg >= 1)
    {
        if (st_.A[P4] <= st_.SRM)
        {
            // (P3 = P2 + 1 AND P1 = (P2 + P3)/2)
            if (std::abs(P3 - (P2 + 1.0)) < 1e-9 &&
                std::abs(P1 - (P2 + P3) * 0.5) < 1e-9)
            {
                st_.T1 = 2.0 * std::log(st_.S[P4] / st_.A[P4]);
                st_.T2 = -st_.W * st_.S[P4];
                return;
            }
        }
    }

    // BASIC 94–98: P1 is half-integer -> midpoint between INT(P1) and INT(P1)+1
    // In scalar calls, P1 is like N+0.5, so INT(P1)=N.
    int I4 = (int)std::floor(P1);
    int I5 = I4 + 1;

    Vec3 X1;
    X1.x = 0.5 * (st_.BX[I4] + st_.BX[I5]);
    X1.y = 0.5 * (st_.BY[I4] + st_.BY[I5]);
    X1.z = 0.5 * (st_.BZ[I4] + st_.BZ[I5]);

    // and then jump into the shared core (BASIC GOTO 113)
    psiCore102_28(X1, P2, P3, P4, Kimg, /*FVS=*/1);
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

    U1 = (st_.T1 - U5) / Ssafe(st_, J2);
    U2 = (st_.T2 - U6) / Ssafe(st_, J2);

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

    U1 = U1 + (U3 - st_.T1) / Ssafe(st_, J1);
    U2 = U2 + (U4 - st_.T2) / Ssafe(st_, J1);

    st_.ZR[I][J] += (double)K * (D1 + U1);
    st_.ZI[I][J] += (double)K * (D2 + U2);
}


// void Mininec3Compat::applySymmetryAndToeplitzCopy(int I, int J, int F8, int J2)
// {
//     if (J < I) return;

//     if (J >= I) {
//         st_.ZR[J][I] = st_.ZR[I][J];
//         st_.ZI[J][I] = st_.ZI[I][J];
//     }
//     if (F8 == 0) return;

//     // if (F8 == 0) return;

//     // st_.ZR[J][I] = st_.ZR[I][J];
//     // st_.ZI[J][I] = st_.ZI[I][J];

//     int P1 = J + 1;
//     if (P1 >= st_.N) return;

//     if (st_.C[P1][0] != st_.C[P1][1]) return;

//     if (st_.C[P1][1] != st_.C[J][1])
//     {
//         if (st_.C[P1][1] != -st_.C[J][1]) return;
//         // if ((st_.CA[J2] + st_.CB[J2]) != 0) return;
//         if (std::abs(st_.CA[J2] + st_.CB[J2]) > 1e-12) return;
//     }

//     int P2 = I + 1;
//     if (P2 >= st_.N) return;

//     st_.ZR[P2][P1] = st_.ZR[I][J];
//     st_.ZI[P2][P1] = st_.ZI[I][J];
// }

void Mininec3Compat::applySymmetryAndToeplitzCopy(int I, int J, int F8, int J2)
{
    if (J < I) return;

    // alltid symmetri
    st_.ZR[J][I] = st_.ZR[I][J];
    st_.ZI[J][I] = st_.ZI[I][J];

    // DEBUG: stäng av Toeplitz-copy
    return;
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

            // Beräkna F8 en gång (för hela (I,J))
            int F8 = computeF8flag(I, J, I1, I2, J1, J2);

            // Nollställ innan vi summerar K
            st_.ZR[I][J] = 0.0;
            st_.ZI[I][J] = 0.0;

            st_.WpulseScratchObs = st_.Wpulse[I];
            st_.WpulseScratchSrc = st_.Wpulse[J];

            if (cJ1 == -cJ2)
            {
                // special: bara K=+1
                F6 = F4;
                F7 = F5;
                computeZijForImage(I, J, +1, I1, I2, J1, J2, F4, F5, F6, F7, F8);
            }
            else
            {
                // normal: summera K=+1 och K=-1
                computeZijForImage(I, J, +1, I1, I2, J1, J2, F4, F5, F6, F7, F8);
                computeZijForImage(I, J, -1, I1, I2, J1, J2, F4, F5, F6, F7, F8);
            }

            // Kopiera först när Z[I][J] är färdig
            applySymmetryAndToeplitzCopy(I, J, F8, J2);
        }
    }

    std::cout << " BuildZ:  No of hits = " << hits << std::endl;

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


