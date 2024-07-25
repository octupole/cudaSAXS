#include "RunSaxs.h"

RunSaxs::RunSaxs(std::string tpr_file, std::string xtc_file)
{
    script = new Python3(tpr_file, xtc_file);
}

void RunSaxs::Run(int beg, int end, int dt)
{
    script->get_atoms(100);
};
RunSaxs::~RunSaxs() {};