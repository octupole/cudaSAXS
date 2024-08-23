#include "RunSaxs.h"
#include <iostream>
#include <vector>
#include <cstdlib> // for std::getenv
#include <CLI/CLI.hpp>
#include "Options.h"
#include "Ftypedefs.h"

/**
 * @brief Main entry point of the cudaSAXS application.
 *
 * This function parses command-line arguments using the CLI11 library and sets up the necessary parameters for running the SAXS simulation. It then creates a `RunSaxs` object and calls its `Run()` method to execute the simulation.
 *
 * The command-line arguments include:
 * - `-s,--topology`: Required option for the topology file (.tpr)
 * - `-x,--trajectory`: Required option for the trajectory file (.xtc)
 * - `-g,--grid`: Required option for the grid length (can be a single value or a vector of 3 values)
 * - `-b,--begin`: Required option for the initial frame
 * - `-e,--end`: Required option for the last frame
 * - `--dt`: Optional option for the frame read interval
 * - `--order`: Optional option for the BSpline order
 * - `--gridS`: Optional option for the scaled grid length (can be a single value or a vector of 3 values)
 * - `--Scale`: Optional option for the multiplicative factor to obtain the largest scaled grid
 * - `--bin,--Dq`: Optional option for the bin size of the histogram
 * - `-h,--help`: Prints the usage information
 *
 * After parsing the command-line arguments, the function sets the appropriate values in the `Options` struct and then creates a `RunSaxs` object to run the SAXS simulation. Finally, it prints a message indicating that the program has completed.
 */
int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string tpr_file;
    std::string xtc_file;
    std::vector<int> myN;
    std::vector<int> myNS;
    int start, end, dt{1};

    // Add required option for topology file
    app.add_option("-s,--topology", tpr_file, "Topology file (.tpr)")->required();

    // Add optional option for trajectory file
    app.add_option("-x,--trajectory", xtc_file, "Trajectory file (.xtc)")->required();

    // Add optional option grid size
    app.add_option("-g,--grid ", myN, "Grid length")->required();

    // Add optional option grid size
    app.add_option("-b,--begin ", start, "Initial frame")->required();

    // Add optional option grid size
    app.add_option("-e,--end ", end, "Last Frame")->required();

    // Add optional option grid size
    app.add_option("--dt ", dt, "Read two frames every dt");

    // Add optional option grid size
    app.add_option("--order ", Options::order, "BSpline order");

    // Add optional option grid size
    app.add_option("--gridS ", myNS, "Scaled Grid length");

    // Add optional option grid scale
    app.add_option("--Scale ", Options::sigma, "Multiplicative factor to obtained the largest scaled grid");

    // Add optional option grid scale
    app.add_option("--bin,--Dq ", Options::Dq, "Binsize of the histogram");

    // Add optional cutoff on reciprocal space
    app.add_option("-q,--qcut ", Options::Qcut, "Cutoff on reciprocal space");

    // Add optional cutoff on reciprocal space
    app.add_option("--water", Options::Wmodel, "Model to use for the weighting function");

    // Add optional cutoff on reciprocal space
    app.add_option("--na ", Options::Sodium, "Sodium atoms");

    // Add optional cutoff on reciprocal space
    app.add_option("--cl ", Options::Chlorine, "Chlorine atoms");

    // Add help option
    app.set_help_flag("-h,--help", "Print usage");

    CLI11_PARSE(app, argc, argv);
    int nx{0}, ny{0}, nz{0};
    if (myN.size() == 1)
    {
        nx = myN[XX];
        ny = myN[XX];
        nz = myN[XX];
    }
    else if (myN.size() == 3)
    {
        nx = myN[XX];
        ny = myN[YY];
        nz = myN[ZZ];
    }
    else
    {
        std::cout << "Invalid number of arguments for grid size. Please provide either 1 or 3 integers." << std::endl;
        return 1;
    }
    int nnx{0}, nny{0}, nnz{0};
    if (myNS.size() == 1)
    {
        nx = myNS[XX];
        ny = myNS[XX];
        nz = myNS[XX];
    }
    else if (myN.size() == 3)
    {
        nx = myNS[XX];
        ny = myNS[YY];
        nz = myNS[ZZ];
    }
    else if (myNS.size() != 0)
    {
        std::cout << "Invalid number of arguments for grid size. Please provide either 1 or 3 integers." << std::endl;
        return 1;
    }
    Options::tpr_file = tpr_file;
    Options::xtc_file = xtc_file;
    Options::nx = nx;
    Options::ny = ny;
    Options::nz = nz;
    Options::nnx = nnx;
    Options::nny = nny;
    Options::nnz = nnz;
    if (Options::Wmodel.size())
    {
        Options::myPadding = padding::given;
    }
    // Get the input string from the command-line arguments
    RunSaxs saxs(tpr_file, xtc_file);
    saxs.Run(start, end, dt);
    // run_python_script_with_args(tpr_file, xtc_file);

    // Print the result from the Python script
    std::cout << "Program completed " << std::endl;

    return 0;
}
