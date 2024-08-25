#ifndef OPTIONS_H
#define OPTIONS_H
#include <string>
#include <vector>
#include <map>
#pragma once
enum class padding
{
    avg,
    given
};

class Options
{
public:
    static std::string tpr_file, xtc_file;
    static int nx, ny, nz;
    static int nnx, nny, nnz;
    static float sigma, Dq, Qcut;
    static int order;
    static std::string Wmodel;
    static int Sodium, Chlorine;
    static padding myPadding;
    static std::map<std::string, float> myWmodel;
    static std::string outFile;

private:
    Options() {};
    ~Options() {};
};

#endif