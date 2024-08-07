#include "RunSaxs.h"
#include <fstream>
#include <fmt/core.h>
#include <fmt/format.h>

/// Creates a vector of integers with a specified start, end, and step.
///
/// This function calculates the size of the vector based on the given start, end, and step values.
/// It then creates a vector and fills it with sequential integers starting from 0, and transforms
/// the values to match the desired sequence.
///
/// @param start The starting value for the sequence.
/// @param end The ending value for the sequence.
/// @param step The step size between values.
/// @return A vector of integers representing the desired sequence.
std::vector<int> RunSaxs::createVector(int start, int end, int step)
{
    // Calculate the size of the vector
    int size = (end - start) / step + 1;

    // Create a vector to hold the values
    std::vector<int> result(size);

    // Fill the vector with sequential integers starting from 0
    std::iota(result.begin(), result.end(), 0);

    // Transform the values to match the desired sequence
    std::transform(result.begin(), result.end(), result.begin(),
                   [start, step](int x)
                   { return start + x * step; });

    return result;
}
/// Runs the SAXS (Small-Angle X-ray Scattering) analysis on a range of frames.
///
/// This function creates a vector of frame indices to process, and then iterates over each frame.
/// For each frame, it retrieves the centered coordinates and box dimensions, calculates the
/// transformation matrices, and runs the SAXS kernel on the coordinates. The elapsed time for
/// the entire process is measured and printed.
///
/// @param beg The starting frame index.
/// @param end The ending frame index.
/// @param dt The step size between frames.
void RunSaxs::Run(int beg, int end, int dt)
{
    auto NToA = [](std::vector<std::vector<float>> &vec)
    {
        for (auto &row : vec)
        {
            for (auto &col : row)
                col = col * 10.0f;
        };
    };
    auto args = createVector(beg, end, dt);
    py::scoped_interpreter guard{}; // Start the interpreter and keep it alive for this scope
    std::string result;
    std::map<std::string, std::vector<std::vector<float>>> coord_map;
    std::map<std::string, std::vector<int>> index_map;
    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")(PY_SOURCE_DIR);
    // Set up the arguments to pass to the Python script

    // Load the Python script
    py::module_ script = py::module_::import("topology");

    // Create a StructureCentering instance (adjust arguments as needed)
    py::object trajStructure = script.attr("Topology");
    py::object analyzer = trajStructure(tpr_file, xtc_file);
    py::dict gather_dict = analyzer.attr("get_atom_index")();
    for (auto item : gather_dict)
    {
        std::string key = py::str(item.first);

        std::vector<int> value = item.second.cast<std::vector<int>>();
        index_map[key] = value;
    }

    auto start = std::chrono::high_resolution_clock::now();
    saxsKernel myKernel(Options::nx, Options::ny, Options::nz, Options::order);
    myKernel.setnpx(8);
    myKernel.scaledCell();
    for (auto frame : args)
    {

        try
        {
            analyzer.attr("read_frame")(frame);
            float simTime = analyzer.attr("get_time")().cast<py::float_>();
            py::object coords_obj = analyzer.attr("get_coordinates")();
            auto coords = analyzer.attr("get_coordinates")().cast<std::vector<std::vector<float>>>();
            NToA(coords);
            auto box_dimensions = analyzer.attr("get_box")().cast<std::vector<std::vector<float>>>();

            Cell::calculateMatrices(box_dimensions);
            auto co = Cell::getCO();
            auto oc = Cell::getOC();
            myKernel.runPKernel(frame, simTime, coords, index_map, oc);
        }

        catch (const py::error_already_set &e)
        {
            std::cerr << "Python error: " << e.what() << std::endl;
        }
    }
    auto myhisto = myKernel.getSaxs();
    std::ofstream myfile;
    myfile.open("saxs.dat");
    for (auto data : myhisto)
    {
        myfile << data[0] << " " << data[1] << std::endl;
    }
    std::cout << "Done " << args.size() << " Steps " << std::endl;
    myfile.close();
    auto end0 = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start);
    auto cudaTime = myKernel.getCudaTime();
    auto totalTime = duration_ms.count() / (float)args.size();
    auto readTime = totalTime - cudaTime;

    std::string banner = fmt::format(
        "\n=========================================================\n"
        "=                                                       =\n"
        "=                   cudaSAXS Timing                     =\n"
        "=                                                       =\n"
        "=           CUDA Time:     {:<10.2f} ms/per step       =\n"
        "=           Read Time:     {:<10.2f} ms/per step       =\n"
        "=           Total Time:    {:<10.2f} ms/per step       =\n"
        "=                                                       =\n"
        "=========================================================\n\n",
        cudaTime, readTime, totalTime);

    fmt::print("{}", banner);
};
RunSaxs::~RunSaxs() {};