#include "RunSaxs.h"
#include <iostream>
#include <vector>
#include <cstdlib> // for std::getenv
#include <CLI/CLI.hpp>
#include "Options.h"

int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string tpr_file;
    std::string xtc_file;
    std::vector<int> myN;
    std::vector<int> myNS;
    float sigma{2.5f};
    float Dq{0.05f};
    // Add required option for topology file
    app.add_option("-s,--topology", tpr_file, "Topology file (.tpr)")->required();

    // Add optional option for trajectory file
    app.add_option("-x,--trajectory", xtc_file, "Trajectory file (.xtc)")->required();

    // Add optional option grid size
    app.add_option("-g,--grid ", myN, "Grid length")->required();

    // Add optional option grid size
    app.add_option("--gridS ", myNS, "Scaled Grid length");

    // Add optional option grid scale
    app.add_option("--Scale ", sigma, "Multiplicative factor to obtained the largest scaled grid");
    // Add optional option grid scale
    app.add_option("--bin,--Dq ", Dq, "Binsize of the histogram");

    // Add help option
    app.set_help_flag("-h,--help", "Print usage");

    CLI11_PARSE(app, argc, argv);
    int nx{0}, ny{0}, nz{0};
    if (myN.size() == 1)
    {
        nx = myN[0];
        ny = myN[0];
        nz = myN[0];
    }
    else if (myN.size() == 3)
    {
        nx = myN[0];
        ny = myN[1];
        nz = myN[2];
    }
    else
    {
        std::cout << "Invalid number of arguments for grid size. Please provide either 1 or 3 integers." << std::endl;
        return 1;
    }
    std::cout << "Grid size: " << nx << "x" << ny << "x" << nz << std::endl;
    std::cout << "Topology file: " << tpr_file << std::endl;
    if (!xtc_file.empty())
    {
        std::cout << "Trajectory file: " << xtc_file << std::endl;
    }
    Options::tpr_file = tpr_file;
    Options::xtc_file = xtc_file;
    Options::nx = nx;
    Options::ny = ny;
    Options::nz = nz;
    Options::sigma = sigma;
    // Get the input string from the command-line arguments
    RunSaxs saxs(tpr_file, xtc_file);
    saxs.Run(100, 200, 10);
    // run_python_script_with_args(tpr_file, xtc_file);

    // Print the result from the Python script
    std::cout << "Program completed " << std::endl;

    return 0;
}
