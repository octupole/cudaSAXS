#include "Atoms.h"

Atoms::Atoms(string lab, vector<vector<float>> xc)
{
    label = lab;
    x = std::vector<std::vector<float>>(xc.size(), std::vector<float>(DIM, 0.0f));

    auto oc = Cell::getOC();

    for (int o = 0; o < xc.size(); o++)
    {
        x[o][XX] = oc[XX][XX] * xc[o][XX] + oc[XX][YY] * xc[o][YY] + oc[XX][ZZ] * xc[o][ZZ];
        x[o][YY] = oc[YY][XX] * xc[o][XX] + oc[YY][YY] * xc[o][YY] + oc[YY][ZZ] * xc[o][ZZ];
        x[o][ZZ] = oc[ZZ][XX] * xc[o][XX] + oc[ZZ][YY] * xc[o][YY] + oc[ZZ][ZZ] * xc[o][ZZ];
    }
}
Atoms::~Atoms()
{
}