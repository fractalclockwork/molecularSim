// File:  simMC.cpp
// -----------------
// This file contains the implementation of the Simulation class,
// which defines all methods and data for a Molecular Simulation,
// using approaches such as Molecular Dynamics, Monte Carlo etc.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "simMC.h"

#include <iostream>
#include <toml++/toml.h>

using namespace std;

// Method:  readSimParameters
// Usage:   readSimParameters();
// -----------------------------
// Reads in parameters for NVE ensemble

void Simulation::readSimParameters() {
    int simulation;

    // Path to the TOML configuration file
    std::string configPath = std::string(DATA_PATH) + "/config_mc.toml";

    cout << "Attempting to parse: " << configPath << endl;

    try {
        // Parse the TOML file
        auto config = toml::parse_file(configPath);

        // Read the simulation choice
        simulation = config["simulation"]["choice"].value_or(0);
        if (simulation <= 0) {
            std::cerr << "Invalid simulation choice in configuration!" << endl;
            return;
        }

        // Process the selected simulation type
        switch (simulation) {
        case 2: // Monte Carlo
            mc = new MonteCarlo(); // Instantiate Monte Carlo object
            mc->run();             // Run Monte Carlo simulation
            break;

        // Placeholder for additional simulation types
        // case 3: ...
        // case 4: ...

        default:
            cout << "Invalid simulation selected! Aborting." << endl;
            break;
        }
    } catch (const toml::parse_error& err) {
        // Handle TOML parsing errors
        std::cerr << "Error parsing TOML file: " << err << std::endl;
    } catch (const std::exception& e) {
        // Handle other unexpected errors
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }

    return;
}