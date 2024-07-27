#ifndef ATOMS_H
#define ATOMS_H
#include <vector>
#include <string>
#include <map>
#include "Cell.h"
#include "Ftypedefs.h"
using namespace std;

class Atoms
{
public:
    Atoms(string, vector<vector<float>>);
    string getLabel() { return label; }
    std::vector<float> operator[](int index) { return x[index]; };
    ~Atoms();

private:
    string label;
    vector<vector<float>> x;
};

#endif