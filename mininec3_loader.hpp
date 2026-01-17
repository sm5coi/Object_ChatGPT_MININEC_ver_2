#ifndef MININEC3_LOADER_HPP
#define MININEC3_LOADER_HPP

#include "mininec3_state.hpp"

#include <string>

struct WireDef
{
    int wireId = 0;
    int nSeg = 0;
    double x1=0, y1=0, z1=0;
    double x2=0, y2=0, z2=0;
    double radius=0.001;
};

class Mininec3Loader
{
public:
    // Reads simple format:
    // wireId nSegments x1 y1 z1 x2 y2 z2 radius
    static bool loadFromFile(const std::string& filename, Mininec3State& st);

private:
    static bool parseFile(const std::string& filename, std::vector<WireDef>& wires);
    //static void buildStateFromWires(const std::vector<WireDef>& wires, Mininec3State& st);

    static void buildMininecConnectivity(const std::vector<WireDef>& wires, Mininec3State& st);
};

#endif // MININEC3_LOADER_HPP
