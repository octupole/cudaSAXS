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
    static float sigma, Dq;
    static std::vector<int> myN;
    static std::vector<int> myNS;

private:
    Options() {};
    ~Options() {};
};

#endif