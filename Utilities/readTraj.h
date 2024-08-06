#ifndef READTRAJ_H
#define READTRAJ_H
#include "FstreamC.h"
#include "Units.h"
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <xdrfile_seek.h>

#pragma once

class readTraj
{
public:
    readTraj(string, int);
    void readFrame(int);
    void readaStep(int);
    vector<vector<float>> &getCoord() { return coords; }
    vector<vector<float>> &getCO() { return co; }
    vector<vector<float>> &getOC();

    ~readTraj();

private:
    XDRFILE *fin{nullptr};
    XDR *xdr{nullptr};

    FILE *fp{nullptr};
    rvec *x0{nullptr};
    matrix box;
    void setCoord(rvec *);
    vector<vector<float>> getInverseMatrix(vector<vector<float>> &);
    FstreamC *mytraj{nullptr};
    vector<vector<float>> coords;
    vector<vector<float>> co{3, vector<float>(3)};
    vector<vector<float>> oc{3, vector<float>(3)};
    int natoms;
    int step_c;
    float time_c, prec_c;
    int dstep{0};
    static bool firsttime;
};
#endif