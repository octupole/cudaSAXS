#ifndef ATOMS_H
#define ATOMS_H
#include <vector>
#include <string>
#include "Cell.h"
#include "Ftypedefs.h"
using namespace std;

class Atoms
{
public:
    Atoms(string, vector<vector<float>>);
    ~Atoms();

private:
    string label;
    vector<vector<float>> x;
};

#endif