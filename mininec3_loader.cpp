#include "mininec3_loader.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

/*
static double dist3(double x1,double y1,double z1,double x2,double y2,double z2)
{
    double dx=x2-x1, dy=y2-y1, dz=z2-z1;
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}
*/

bool Mininec3Loader::loadFromFile(const std::string& filename, Mininec3State& st)
{
    std::vector<WireDef> wires;
    if (!parseFile(filename, wires))
        return false;

    buildMininecConnectivity(wires, st);
    return true;
}

bool Mininec3Loader::parseFile(const std::string& filename, std::vector<WireDef>& wires)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "Cannot open file: " << filename << "\n";
        return false;
    }

    wires.clear();

    std::string line;
    int lineNo = 0;
    while (std::getline(in, line))
    {
        lineNo++;

        // strip comments
        auto hashPos = line.find('#');
        if (hashPos != std::string::npos)
            line = line.substr(0, hashPos);

        // trim
        auto isSpace = [](unsigned char c){ return std::isspace(c); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [&](char c){ return !isSpace(c); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [&](char c){ return !isSpace(c); }).base(), line.end());

        if (line.empty())
            continue;

        std::istringstream iss(line);
        WireDef w;
        if (!(iss >> w.wireId >> w.nSeg >> w.x1 >> w.y1 >> w.z1 >> w.x2 >> w.y2 >> w.z2 >> w.radius))
        {
            std::cerr << "Parse error line " << lineNo << ": " << line << "\n";
            return false;
        }

        if (w.nSeg < 1)
        {
            std::cerr << "Invalid nSeg at line " << lineNo << "\n";
            return false;
        }

        wires.push_back(w);
    }

    // Ensure wires are sorted by wireId and contiguous (0..W-1)
    std::sort(wires.begin(), wires.end(), [](const WireDef& a, const WireDef& b){
        return a.wireId < b.wireId;
    });

    for (int i = 0; i < (int)wires.size(); ++i)
    {
        if (wires[i].wireId != i)
        {
            std::cerr << "Wire IDs must be contiguous starting at 0. Missing wireId=" << i << "\n";
            return false;
        }
    }

    return true;
}

/*
void Mininec3Loader::buildStateFromWires(const std::vector<WireDef>& wires, Mininec3State& st)
{
    // Reset state
    st = Mininec3State();

    int nWires = (int)wires.size();
    st.wireNodes.resize(nWires);
    st.wireBaseP.resize(nWires);

    // We'll build global nodes and global segments
    int nodeCounter = 0;
    int segCounter  = 0;

    // We need per-segment arrays sized later
    // But we also need mapping from (wire, localSeg) -> globalSegId
    std::vector<std::vector<int>> wireSegId(nWires);

    // ---- Build nodes + wireNodes + wireBaseP ----
    for (int w = 0; w < nWires; ++w)
    {
        const auto& wd = wires[w];

        int nNodes = wd.nSeg + 1;

        // BaseP is global node offset in our "global P space"
        st.wireBaseP[w] = nodeCounter;

        st.wireNodes[w].resize(nNodes);

        for (int i = 0; i < nNodes; ++i)
        {
            double t = (wd.nSeg == 0) ? 0.0 : (double)i / (double)wd.nSeg;

            double x = wd.x1 + (wd.x2 - wd.x1) * t;
            double y = wd.y1 + (wd.y2 - wd.y1) * t;
            double z = wd.z1 + (wd.z2 - wd.z1) * t;

            int nodeId = nodeCounter++;
            st.wireNodes[w][i] = nodeId;

            st.X.push_back(x);
            st.Y.push_back(y);
            st.Z.push_back(z);
        }

        wireSegId[w].resize(wd.nSeg);
        for (int s = 0; s < wd.nSeg; ++s)
            wireSegId[w][s] = segCounter++;
    }

    int totalSeg = segCounter;

    // Allocate per-segment arrays
    st.S.assign(totalSeg, 0.0);
    st.A.assign(totalSeg, 0.0);
    st.CA.assign(totalSeg, 0.0);
    st.CB.assign(totalSeg, 0.0);
    st.CG.assign(totalSeg, 0.0);

    // ---- Fill per-segment arrays ----
    for (int w = 0; w < nWires; ++w)
    {
        const auto& wd = wires[w];

        for (int s = 0; s < wd.nSeg; ++s)
        {
            int n0 = st.wireNodes[w][s];
            int n1 = st.wireNodes[w][s + 1];

            double x0 = st.X[n0], y0 = st.Y[n0], z0 = st.Z[n0];
            double x1 = st.X[n1], y1 = st.Y[n1], z1 = st.Z[n1];

            double len = dist3(x0,y0,z0, x1,y1,z1);
            if (len < 1e-12) len = 1e-12;

            double dx = (x1 - x0) / len;
            double dy = (y1 - y0) / len;
            double dz = (z1 - z0) / len;

            int segId = wireSegId[w][s];

            st.S[segId]  = len;
            st.A[segId]  = wd.radius;
            st.CA[segId] = dx;
            st.CB[segId] = dy;
            st.CG[segId] = dz;
        }
    }

    // ---- Build pulses C% and W% ----
    // For each wire with Ns segments => Np = Ns - 1 pulses
    // Pulse p uses segments (s, s+1)
    std::vector<std::array<int,2>> C;
    std::vector<int> Wpulse;

    for (int w = 0; w < nWires; ++w)
    {
        int Ns = wires[w].nSeg;
        int Np = std::max(0, Ns - 1);

        for (int p = 0; p < Np; ++p)
        {
            int seg1 = wireSegId[w][p];
            int seg2 = wireSegId[w][p + 1];

            C.push_back({seg1, seg2});
            Wpulse.push_back(w);
        }
    }

    st.C = std::move(C);
    st.Wpulse = std::move(Wpulse);
    st.N = (int)st.C.size();

    // Allocate Z matrices
    st.ZR.assign(st.N, std::vector<double>(st.N, 0.0));
    st.ZI.assign(st.N, std::vector<double>(st.N, 0.0));

    // SRM is typically 2*radius in MININEC
    // Here we just set it to 2*min(radius) as a first approximation
    double minR = 1e9;
    for (const auto& wd : wires)
        minR = std::min(minR, wd.radius);
    if (minR > 0 && minR < 1e8)
        st.SRM = 2.0 * minR;
    else
        st.SRM = 1e-6;
}
*/

static int sgnInt(int v) { return (v >= 0) ? +1 : -1; }

void Mininec3Loader::buildMininecConnectivity(const std::vector<WireDef>& wires, Mininec3State& st)
{
    // st = Mininec3State();

    st.NW = (int)wires.size();
    st.Gmode = 1; // user requested G=1

    int NW = st.NW;

    // BASIC arrays are 1-based for wires. We'll keep wire index 1..NW in logic,
    // but store in vectors with size NW+1.
    st.E.assign(2*NW + 1, 0.0);
    st.L.assign(2*NW + 1, 0.0);
    st.M.assign(2*NW + 1, 0.0);

    st.J1.assign(NW + 1, 0);
    st.J2.assign(NW + 1, {0,0});

    st.wirePulseRange.assign(NW + 1, {0,0});

    // Per-wire CA/CB/CG/S/A like BASIC uses index = wireId
    st.CA.assign(NW + 1, 0.0);
    st.CB.assign(NW + 1, 0.0);
    st.CG.assign(NW + 1, 0.0);
    st.S.assign(NW + 1, 0.0);
    st.A.assign(NW + 1, 0.0);

    // Pulse arrays: we don't know N yet, build incrementally like BASIC
    st.C.clear();
    st.Wpulse.clear();

    // Breakpoints BX/BY/BZ: BASIC uses indices around:
    // I1 = N1 + 2*(I-1)
    // I6 = N  + 2*I
    // So allocate "big enough" as we go (push_back won't match indices),
    // we will resize to max index used.
    st.BX.assign(1, 0.0);
    st.BY.assign(1, 0.0);
    st.BZ.assign(1, 0.0);

    auto ensureBP = [&](int idx){
        if (idx >= (int)st.BX.size())
        {
            st.BX.resize(idx + 1, 0.0);
            st.BY.resize(idx + 1, 0.0);
            st.BZ.resize(idx + 1, 0.0);
        }
    };

    int N = 0; // BASIC pulse count

    for (int I = 1; I <= NW; ++I)
    {
        const auto& wd = wires[I-1];

        int S1 = wd.nSeg;
        double X1 = wd.x1, Y1 = wd.y1, Z1 = wd.z1;
        double X2 = wd.x2, Y2 = wd.y2, Z2 = wd.z2;

        st.A[I] = wd.radius;

        // --- BASIC 1299-1304: store endpoints
        st.E[I]      = X1;
        st.L[I]      = Y1;
        st.M[I]      = Z1;
        st.E[I + NW] = X2;
        st.L[I + NW] = Y2;
        st.M[I + NW] = Z2;

        // --- BASIC 1305-1310 init connections
        int I1 = 0;
        int I2 = 0;
        st.J1[I] = 0;
        st.J2[I][0] = -I;
        st.J2[I][1] = -I;

        // --- BASIC 1311: IF G = 1 THEN 1323
        // We are G=1 so we skip ground connection checks entirely.

        // --- BASIC 1323..1357: check connections to previous wires
        if (I != 1)
        {
            // END1 comparisons
            for (int J = 1; J <= I-1; ++J)
            {
                // END1 to END1
                if (X1 == st.E[J] && Y1 == st.L[J] && Z1 == st.M[J])
                {
                    I1 = -J;
                    st.J2[I][0] = J;
                    if (st.J2[J][0] == -J) st.J2[J][0] = J;
                    goto end1_done;
                }

                // END1 to END2
                if (X1 == st.E[J+NW] && Y1 == st.L[J+NW] && Z1 == st.M[J+NW])
                {
                    I1 = J;
                    st.J2[I][0] = J;
                    if (st.J2[J][1] == -J) st.J2[J][1] = J;
                    goto end1_done;
                }
            }
        }
    end1_done:

        // END2 comparisons
        if (I != 1)
        {
            for (int J = 1; J <= I-1; ++J)
            {
                // END2 to END2
                if (X2 == st.E[J+NW] && Y2 == st.L[J+NW] && Z2 == st.M[J+NW])
                {
                    I2 = -J;
                    st.J2[I][1] = J;
                    if (st.J2[J][1] == -J) st.J2[J][1] = J;
                    goto end2_done;
                }

                // END2 to END1
                if (X2 == st.E[J] && Y2 == st.L[J] && Z2 == st.M[J])
                {
                    I2 = J;
                    st.J2[I][1] = J;
                    if (st.J2[J][0] == -J) st.J2[J][0] = J;
                    goto end2_done;
                }
            }
        }
    end2_done:

        // --- BASIC 1190-1197 direction cosines and segment length
        double X3 = X2 - X1;
        double Y3 = Y2 - Y1;
        double Z3 = Z2 - Z1;
        double D = std::sqrt(X3*X3 + Y3*Y3 + Z3*Z3);
        if (D < 1e-12) D = 1e-12;

        st.CA[I] = X3 / D;
        st.CB[I] = Y3 / D;
        st.CG[I] = Z3 / D;
        st.S[I]  = D / (double)S1;

        // --- BASIC 1198-1208: compute pulse range N1..N
        int N1 = N + 1;
        st.wirePulseRange[I][0] = N1;

        // BASIC 1201..1207 handle special case S1=1 and I1/I2=0
        // For now we keep the same behavior.
        if (S1 == 1 && I1 == 0) st.wirePulseRange[I][0] = 0;

        N = N1 + S1;
        if (I1 == 0) N = N - 1;
        if (I2 == 0) N = N - 1;

        st.wirePulseRange[I][1] = N;
        if (S1 == 1 && I2 == 0) st.wirePulseRange[I][1] = 0;

        // If N < N1 => no pulses (single segment, no connections) special case
        if (N < N1)
        {
            // BASIC 1247..1254 single segment 0 pulse case
            int bp1 = N1 + 2*(I-1);
            ensureBP(bp1);
            st.BX[bp1] = X1; st.BY[bp1] = Y1; st.BZ[bp1] = Z1;

            int bp2 = bp1 + 1;
            ensureBP(bp2);
            st.BX[bp2] = X2; st.BY[bp2] = Y2; st.BZ[bp2] = Z2;

            continue;
        }

        // --- BASIC 1209-1213: fill pulses J=N1..N with default C=(I,I), W=I
        // We'll append to st.C and st.Wpulse but remember BASIC uses 1-based pulse indices.
        // We'll store pulses 0-based but keep the values inside C% as signed wire indices like BASIC.
        int pulsesOnWire = N - N1 + 1;

        int pulseBaseIndex0 = (int)st.C.size(); // 0-based index in vector

        for (int j = 0; j < pulsesOnWire; ++j)
        {
            st.C.push_back({I, I});   // default
            st.Wpulse.push_back(I);   // BASIC wire number
        }

        // --- BASIC 1214-1215: modify first/last pulse connection flags
        // First pulse: C%(N1,1)=I1
        // Last  pulse: C%(N,2)=I2
        st.C[pulseBaseIndex0 + 0][0] = I1;
        st.C[pulseBaseIndex0 + (pulsesOnWire - 1)][1] = I2;

        // --- BASIC 1216-1245: compute breakpoint coordinates
        int bpStart = N1 + 2*(I-1);
        int I3 = bpStart;

        ensureBP(bpStart);
        st.BX[bpStart] = X1;
        st.BY[bpStart] = Y1;
        st.BZ[bpStart] = Z1;

        // If first pulse has connection, extend backward along connected wire direction
        if (st.C[pulseBaseIndex0 + 0][0] != 0)
        {
            int I2abs = std::abs(st.C[pulseBaseIndex0 + 0][0]);
            double F3 = (double)sgnInt(st.C[pulseBaseIndex0 + 0][0]) * st.S[I2abs];

            st.BX[bpStart] -= F3 * st.CA[I2abs];
            st.BY[bpStart] -= F3 * st.CB[I2abs];

            if (st.C[pulseBaseIndex0 + 0][0] == -I) F3 = -F3;
            st.BZ[bpStart] -= F3 * st.CG[I2abs];

            I3 = I3 + 1;
        }

        int bpEnd = N + 2*I;
        ensureBP(bpEnd);

        for (int bp = bpStart + 1; bp <= bpEnd; ++bp)
        {
            int J = bp - I3;
            st.BX[bp] = X1 + (double)J * X3 / (double)S1;
            st.BY[bp] = Y1 + (double)J * Y3 / (double)S1;
            st.BZ[bp] = Z1 + (double)J * Z3 / (double)S1;
        }

        // If last pulse has connection, extend forward along connected wire direction
        if (st.C[pulseBaseIndex0 + (pulsesOnWire - 1)][1] != 0)
        {
            int I2abs = std::abs(st.C[pulseBaseIndex0 + (pulsesOnWire - 1)][1]);
            double F3 = (double)sgnInt(st.C[pulseBaseIndex0 + (pulsesOnWire - 1)][1]) * st.S[I2abs];

            int I3b = bpEnd - 1;
            st.BX[bpEnd] = st.BX[I3b] + F3 * st.CA[I2abs];
            st.BY[bpEnd] = st.BY[I3b] + F3 * st.CB[I2abs];

            if (I == -st.C[pulseBaseIndex0 + (pulsesOnWire - 1)][1]) F3 = -F3;
            st.BZ[bpEnd] = st.BZ[I3b] + F3 * st.CG[I2abs];
        }
    }

    // Final pulse count
    st.N = (int)st.C.size();

    // Allocate Z matrices
    st.ZR.assign(st.N, std::vector<double>(st.N, 0.0));
    st.ZI.assign(st.N, std::vector<double>(st.N, 0.0));

    // ---- SRM (MININEC effective radius floor) ----
    double minR = 1e9;
    for (int i = 1; i <= st.NW; ++i)
        minR = std::min(minR, st.A[i]);

    if (minR > 0.0 && minR < 1e8)
        st.SRM = 2.0 * minR;
    else
        st.SRM = 1e-6;

}
