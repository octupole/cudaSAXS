#include "Atoms.h"

Atoms::Atoms(string lab, vector<vector<float>> xc)
{
    label = lab;
    x = std::vector<std::vector<float>>(xc.size(), std::vector<float>(DIM, 0.0f));

    auto oc = Cell::getOC();

    for (int o = 0; o < xc.size(); o++)
    {
        x[XX][o] = oc[XX][XX] * xc[o][XX] + oc[XX][YY] * xc[o][YY] + oc[XX][ZZ] * xc[o][ZZ];
        x[YY][o] = oc[YY][XX] * xc[o][XX] + oc[YY][YY] * xc[o][YY] + oc[YY][ZZ] * xc[o][ZZ];
        x[ZZ][o] = oc[ZZ][XX] * xc[o][XX] + oc[ZZ][YY] * xc[o][YY] + oc[ZZ][ZZ] * xc[o][ZZ];
    }
}

Atoms::~Atoms()
{
}