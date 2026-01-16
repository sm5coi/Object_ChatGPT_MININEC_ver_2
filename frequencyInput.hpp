#ifndef FREQUENCYINPUT_H
#define FREQUENCYINPUT_H

#include <iostream>
#include <string>
#include "mininec3_state.hpp"

class frequencyInput {

public:
    // Konstruktor
    frequencyInput(double freq, const std::string& unitType)
        : frequency(freq), unit(unitType)
    {
        std::cout << "FrequencyInput-objekt skapat med frekvens: "
                  << frequency << " " << unit << std::endl;

        if (unit == "Hz") {
            std::cout << "You chose START\n";
            multip_factor = 1.0;

        } else if (unit == "kHz") {
            std::cout << "You chose STOP\n";
        } else if (unit == "Mhz") {
            std::cout << "You chose PAUSE\n";
        } else if (unit == "Ghz") {
            std::cout << "You chose PAUSE\n";


        } else {
            std::cout << "Unknown command\n";
        }



    }

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

    void set_fq(double fq);

private:
    double frequency; // Frekvens (i Hertz)
    std::string unit; // Enhet (t.ex., "Hz", "kHz")
    double multip_factor; // Faktor för att erhålla frekvens i Hz

};

#endif // FREQUENCYINPUT_H
