#ifndef PYTHON3_H
#define PYTHON3_H
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <vector>
#include <map>
#include "Atoms.h"
#include "Cell.h"

#pragma once
namespace py = pybind11;

class Python3
{
public:
    Python3(const std::string &tpr, const std::string &traj) : tpr_file(tpr), xtc_file(traj) {};
    void get_atoms(int);

    ~Python3();

private:
    std::string tpr_file, xtc_file;
    std::vector<Atoms> atoms;
    std::vector<std::vector<float>> box_dimensions;
};

#endif