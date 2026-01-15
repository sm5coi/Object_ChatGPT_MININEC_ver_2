#include <iostream>
#include <iomanip>
#include <cmath>
#include "mininec3_compat.hpp"
#include "mininec3_loader.hpp"

// st.Qbasic = {
//     0.288675135, 0.5,
//     0.430568156, 0.173927423,
//     0.169990522, 0.326072577,
//     0.480144928, 0.050614268,
//     0.398333239, 0.111190517,
//     0.262766205, 0.156853323,
//     0.091717321, 0.181341892
// };


int main()
{
    Mininec3State st;

    if (!Mininec3Loader::loadFromFile("geom.txt", st))
    {
        std::cerr << "Failed to load geom.txt\n";
        return 1;
    }

    std::cout << "\nBreakpoints:\n";
    for (int i = 0; i < (int)st.BX.size(); ++i)
    {
        std::cout << "  BP[" << i << "] = "
                  << st.BX[i] << ", "
                  << st.BY[i] << ", "
                  << st.BZ[i] << "\n";
    }

    // Physics setup
    double freq = 10e6;
    double c0 = 299792458.0;
    st.W = 2.0 * M_PI * freq / c0;

    // W2 is not correct MININEC scaling yet; keep 1 for now
    st.W2 = 1.0;

    std::cout << "Loaded geometry:\n";
    std::cout << "  NW(wires) = " << st.NW << "\n";
    std::cout << "  pulsesN   = " << st.N << "\n";
    std::cout << "  BX size   = " << st.BX.size() << "\n";
    std::cout << "  SRM       = " << st.SRM << "\n";

    Mininec3Compat solver(st);

    std::cout << "C matrix:\n";
    for (int p = 0; p < st.N; ++p)
    {
        std::cout << "pulse " << (p+1)
        << " W=" << st.Wpulse[p]
        << " C1=" << st.C[p][0]
        << " C2=" << st.C[p][1]
        << "\n";
    }
    std::cout << "S(wire):\n";
    for (int w = 1; w <= st.NW; ++w)
    {
        std::cout << "wire " << w << " S=" << st.S[w]
                  << " CA=" << st.CA[w]
                  << " CB=" << st.CB[w]
                  << " CG=" << st.CG[w]
                  << "\n";
    }

    solver.buildZ();

    std::cout << "\nZR:\n";
    for (int i = 0; i < st.N; ++i)
    {
        for (int j = 0; j < st.N; ++j)
            std::cout << std::setw(14) << st.ZR[i][j] << " ";
        std::cout << "\n";
    }

    std::cout << "\nZI:\n";
    for (int i = 0; i < st.N; ++i)
    {
        for (int j = 0; j < st.N; ++j)
            std::cout << std::setw(14) << st.ZI[i][j] << " ";
        std::cout << "\n";
    }

    return 0;
}


/*
#include "mininec3_compat.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

static void setupTiny3SegDipole(Mininec3State& st)
{
    // --- Geometry: 1 wire, 4 nodes along z, 3 segments ---
    // node 0: z=-1
    // node 1: z=-0.333...
    // node 2: z=+0.333...
    // node 3: z=+1
    st.X = {0,0,0,0};
    st.Y = {0,0,0,0};
    st.Z = {-1.0, -1.0/3.0, +1.0/3.0, +1.0};

    st.wireNodes.resize(1);
    st.wireNodes[0] = {0,1,2,3};
    st.wireBaseP = {0};  // wire 0 starts at P=0


    // --- Segments: 3 ---
    int nSeg = 3;
    st.S.assign(nSeg, 2.0/3.0); // each segment length
    st.A.assign(nSeg, 0.001);   // radius dummy

    // All segments are along +z direction
    st.CA.assign(nSeg, 0.0);
    st.CB.assign(nSeg, 0.0);
    st.CG.assign(nSeg, 1.0);

    // --- Pulses: for 3 segments -> 2 pulses ---
    st.N = 2;
    st.G = 1; // no ground image for now

    st.C.resize(st.N);

    // Pulse 0 uses seg 0 and seg 1
    st.C[0][0] = 0;
    st.C[0][1] = 1;

    // Pulse 1 uses seg 1 and seg 2
    st.C[1][0] = 1;
    st.C[1][1] = 2;

    // All pulses on wire 0
    st.Wpulse.assign(st.N, 0);

    // --- Physics ---
    double freq = 10e6;
    double c0 = 299792458.0;
    st.W = 2.0 * M_PI * freq / c0; // k

    // W2 scaling in BASIC is not exactly k*mu/4pi etc.
    // For now set 1 so we see numbers.
    st.W2 = 1.0;
    st.SRM = 0.002;

    // --- Matrices ---
    st.ZR.assign(st.N, std::vector<double>(st.N, 0.0));
    st.ZI.assign(st.N, std::vector<double>(st.N, 0.0));
}

int main()
{
    Mininec3State st;
    setupTiny3SegDipole(st);

    Mininec3Compat solver(st);
    solver.buildZ();

    std::cout << "ZR:\n";
    for (int i = 0; i < st.N; ++i)
    {
        for (int j = 0; j < st.N; ++j)
            std::cout << std::setw(14) << st.ZR[i][j] << " ";
        std::cout << "\n";
    }

    std::cout << "\nZI:\n";
    for (int i = 0; i < st.N; ++i)
    {
        for (int j = 0; j < st.N; ++j)
            std::cout << std::setw(14) << st.ZI[i][j] << " ";
        std::cout << "\n";
    }

    return 0;
}

*/
