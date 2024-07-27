#ifndef PARTICLES_H
#define PARTICLES_H
#include <cuda_runtime.h>
#pragma once

class Particle
{
public:
    float x, y, z;

    __host__ __device__ Particle(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
    ~Particle();

private:
};

#endif