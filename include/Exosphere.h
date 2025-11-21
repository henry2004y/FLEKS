#ifndef _EXOSPHERE_H_
#define _EXOSPHERE_H_

#include <string>
#include <vector>
#include <cmath>

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>

#include "ReadParam.h"
#include "Grid.h"
#include "Constants.h"

class Exosphere {
private:
    // Configuration parameters
    std::string neutral_profile = "exponential";
    double n0 = 0.0; // Number density at reference radius [m^-3]
    double H0 = 0.0; // Scale height [m]
    double T0 = 0.0; // Temperature [K]
    double k0 = 0.0; // Power-law exponent
    double exobase_radius = 0.0; // Altitude below which neutral density is zero
    double total_production_rate = 0.0; // Total photoionization rate [#/s]

    // Neutral density field (stored on cell centers)
    amrex::Vector<amrex::MultiFab> neutralDensity;

    // Pointer to the Grid object to access geometry and hierarchy
    Grid* grid = nullptr;

public:
    Exosphere(Grid* gridIn) : grid(gridIn) {}
    ~Exosphere() = default;

    // Read parameters from input file
    void read_param(const std::string& command, ReadParam& param);

    // Initialize the neutral density profile
    void initialize();

    // Get neutral density at a specific level and cell index
    amrex::Real get_density(const int iLev, const amrex::MFIter& mfi, const amrex::IntVect& ijk) const;

    // Get temperature
    double get_temperature() const { return T0; }

private:
    // Helper function to calculate profile value at a given position
    double calculate_profile(const amrex::RealVect& pos) const;
};

#endif
