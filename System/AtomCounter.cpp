// AtomCounter.cpp
#include "AtomCounter.h"
#include <iostream>
#include <cmath>

const float AtomCounter::AVOGADRO = 6.022e23f;
const float AtomCounter::WATER_MOLAR_MASS = 18.015f; // g/mol

const std::map<std::string, float> AtomCounter::water_models = {
    {"SPCE", 0.998f}, // g/cm³
    {"TIP3P", 0.982f} // g/cm³
};

AtomCounter::AtomCounter(float lx, float ly, float lz,
                         int sodium, int chlorine, const std::string &model,
                         int gx, int gy, int gz)
    : cell_volume(lx * ly * lz),
      added_sodium(sodium), added_chlorine(chlorine), water_model(model),
      grid_x(gx), grid_y(gy), grid_z(gz) {}

float AtomCounter::calculateWaterMolecules() const
{
    float volume_cm3 = 1000.0f * cell_volume * 1e-24f;            // convert Å³ to cm³
    float water_mass = water_models.at(water_model) * volume_cm3; // mass of water in g
    return (water_mass / WATER_MOLAR_MASS) * AVOGADRO;
}

std::map<std::string, float> AtomCounter::calculateAtomCounts() const
{
    std::string used_model = water_model;
    if (water_models.find(water_model) == water_models.end())
    {
        std::cout << "Invalid water model. Using SPC/E as default." << std::endl;
        used_model = "SPCE";
    }

    float water_molecules = calculateWaterMolecules();
    int total_grid_points = grid_x * grid_y * grid_z;

    std::map<std::string, float> atom_counts;
    atom_counts["O"] = water_molecules / total_grid_points;
    atom_counts["H"] = 2.0f * water_molecules / total_grid_points;
    atom_counts["Na"] = static_cast<float>(added_sodium) / total_grid_points;
    atom_counts["Cl"] = static_cast<float>(added_chlorine) / total_grid_points;

    return atom_counts;
}