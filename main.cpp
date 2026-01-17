#include <iostream>
#include <iomanip>
#include <cmath>
#include "mininec3_compat.hpp"
#include "mininec3_loader.hpp"
#include "mininec3_state.hpp"
#include "frequencyInput.hpp"

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

    frequencyInput input1(299.792458, "MHz", st);

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
