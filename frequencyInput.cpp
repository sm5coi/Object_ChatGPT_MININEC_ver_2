#include "frequencyInput.hpp"

// Konstruktorimplementation
frequencyInput::frequencyInput(double freq, const std::string& unitType, Mininec3State& st)
    : frequency(0), unit(unitType) // Initiera datafält
{
    // Beräkna frekvensen beroende på enheten med hjälp av en privat funktion
    frequency = convertFrequencyToHz(freq, unitType);

    st.fq = frequency;
    st.lambda = st.c0/st.fq;
    // 1140 REM ----- 1 / (4 * PI * OMEGA * EPSILON)
    st.Mw = 4.77783352*st.lambda;
    st.W = 2.0 * M_PI * st.fq/st.c0;
    st.W2 = (st.W * st.W)/2.0;

    st.SRM = 0.0001*st.lambda;
}

// Privat metod: Konvertera frekvens till Hertz
double frequencyInput::convertFrequencyToHz(double freq, const std::string& unitType) {
    if (unitType == "Hz") {
        return freq; // Ingen konvertering behövs
    } else if (unitType == "kHz") {
        return freq * 1e3; // kHz till Hz
    } else if (unitType == "MHz") {
        return freq * 1e6; // MHz till Hz
    } else if (unitType == "GHz") {
        return freq * 1e9; // MHz till Hz
    } else {
        std::cerr << "Ogiltig enhet: " << unitType << ". Använder Hz som standard.\n";
        return freq; // Standard till Hz
    }
}

/*
1131 REM ********** FREQUENCY INPUT **********
1132 REM ----- SET FLAG
1133 PRINT
1134 INPUT "FREQUENCY (MHZ)"; F
1135 IF F = 0 THEN F = 299.8
1136 IF O$ > "C" THEN PRINT #3, " ": PRINT #3, "FREQUENCY (MHZ):"; F
1137 W = 299.8 / F
1138 REM -----VIRTUAL DIPOLE LENGTH FOR NEAR FIELD CALCULATION
1139 S0 = .001 * W
1140 REM ----- 1 / (4 * PI * OMEGA * EPSILON)
1141 M = 4.77783352# * W
1142 REM ----- SET SMALL RADIUS MODIFICATION CONDITION
1143 SRM = .0001 * W
1144 PRINT #3, "    WAVE LENGTH = "; W; " METERS"
1145 REM ----- 2 PI / WAVELENGTH
1146 W = 2 * P / W
1147 W2 = W * W / 2
1148 FLG = 0
1149 RETURN
*/
