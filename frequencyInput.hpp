#ifndef FREQUENCYINPUT_H
#define FREQUENCYINPUT_H

#include <iostream>
#include <string>
#include <cmath>
#include "mininec3_state.hpp"

class frequencyInput {

public:
    // Konstruktor
    frequencyInput(double freq, const std::string& unitType, Mininec3State& st);
      //  : frequency(freq), unit(unitType);


    // Getter-funktioner (valfritt)
    double getFrequency() const {
        return frequency;
    }

    std::string getUnit() const {
        return unit;
    }

    // Funktion för att visa information
    void display() const {
        std::cout << "Frekvens: " << frequency << " " << unit << std::endl;
    }

private:
    double frequency; // Frekvens (i Hertz)
    std::string unit; // Enhet (exempelvis "Hz", "kHz")

    // Privat metod som kan användas för interna beräkningar
    double convertFrequencyToHz(double freq, const std::string& unitType);

};

#endif // FREQUENCYINPUT_H
