#ifndef CUDASAXSMANAGER_H
#define CUDASAXSMANAGER_H
#include <iostream>
#include "Array.h"
#include "Atoms.h"
#include <cufft.h> // Include cuFFT header for cufftComplex

#pragma once

class CudaSaxsManager
{
public:
    CudaSaxsManager(Atoms &);
    ~CudaSaxsManager();

private:
    Array::Array3<cufftComplex> *Rho_p;
};

#endif