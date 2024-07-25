#include "RunSaxs.h"
#include <iostream>
#include <vector>
#include <cstdlib> // for std::getenv
#include <CLI/CLI.hpp>

int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string tpr_file;
    std::string xtc_file;

    // Add required option for topology file
    app.add_option("-s,--topology", tpr_file, "Topology file (.tpr)")->required();

    // Add optional option for trajectory file
    app.add_option("-x,--trajectory", xtc_file, "Trajectory file (.xtc)")->required();

    // Add help option
    app.set_help_flag("-h,--help", "Print usage");

    CLI11_PARSE(app, argc, argv);

    std::cout << "Topology file: " << tpr_file << std::endl;
    if (!xtc_file.empty())
    {
        std::cout << "Trajectory file: " << xtc_file << std::endl;
    }

    // Get the input string from the command-line arguments
    RunSaxs saxs(tpr_file, xtc_file);
    saxs.Run(100, 200, 10);
    // run_python_script_with_args(tpr_file, xtc_file);

    // Print the result from the Python script
    std::cout << "Program completed " << std::endl;

    return 0;
}
