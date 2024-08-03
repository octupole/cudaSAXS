#ifndef OPTIONS_H
#define OPTIONS_H
#include <string>
#include <vector>
#pragma once

class Options
{
public:
    static std::string tpr_file, xtc_file;
    static int nx, ny, nz;
    static int nnx, nny, nnz;
    static float sigma, Dq;
    static int order;

private:
    Options() {};
    ~Options() {};
};

#endif