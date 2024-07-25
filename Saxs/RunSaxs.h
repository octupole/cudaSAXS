#ifndef RUNSAXS_H
#define RUNSAXS_H
#include "Python3.h"
#include <string>

#pragma once

class RunSaxs
{
public:
    RunSaxs(std::string, std::string);

    void Run(int, int, int);
    ~RunSaxs();

private:
    Python3 *script{nullptr};
};

#endif