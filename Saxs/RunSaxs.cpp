#include "RunSaxs.h"

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
void RunSaxs::Run(int beg, int end, int dt)
{
    std::cout << "Running SAXS" << Options::nx << ::endl;
    auto args = createVector(beg, end, dt);
    py::scoped_interpreter guard{}; // Start the interpreter and keep it alive for this scope
    std::string result;
    std::map<std::string, std::vector<std::vector<float>>> coord_map;
    std::map<std::string, std::vector<int>> index_map;
    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")(PY_SOURCE_DIR);
    // Set up the arguments to pass to the Python script

    // Load the Python script
    py::module_ script = py::module_::import("trajectories");

    // Create a StructureCentering instance (adjust arguments as needed)
    py::object trajStructure = script.attr("TrajectoryStructures");
    py::object analyzer = trajStructure(tpr_file, xtc_file);
    py::dict gather_dict = analyzer.attr("get_atom_index")();
    for (auto item : gather_dict)
    {
        std::string key = py::str(item.first);

        std::vector<int> value = item.second.cast<std::vector<int>>();
        index_map[key] = value;
    }

    auto start = std::chrono::high_resolution_clock::now();

    for (auto frame : args)
    {

        try
        {
            std::cout << "Frame: " << frame << std::endl;
            py::object dims = analyzer.attr("get_dimensions_frame")(frame);
            auto box_dimensions = dims.cast<std::vector<float>>();
            Cell::calculateMatrices(box_dimensions);
            auto &co = Cell::getCO();
            auto &oc = Cell::getOC();

            py::object coords_obj = analyzer.attr("get_centered_frame")(frame);
            auto coords = coords_obj.cast<std::vector<std::vector<float>>>();
            size_t O{0};
            for (const auto &pair : index_map)
            {
                std::vector<std::vector<float>> vec0;
                std::string key = pair.first;
                std::vector<int> value = pair.second;
                std::transform(value.begin(), value.end(), std::back_inserter(vec0), [&coords](int i)
                               { return coords[i]; });
                Atoms myAtoms(key, vec0);
                atoms.push_back(myAtoms);
            }
            atoms.clear();
        }

        catch (const py::error_already_set &e)
        {
            std::cerr << "Python error: " << e.what() << std::endl;
        }
    }
    std::cout << "Done " << args.size() << std::endl;
    auto end0 = std::chrono::high_resolution_clock::now();

    auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end0 - start);
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start);
    auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(end0 - start);

    std::cout << "Elapsed time: " << duration_ns.count() / (float)args.size() << " nanoseconds\n";
    std::cout << "Elapsed time: " << duration_us.count() / (float)args.size() << " microseconds\n";
    std::cout << "Elapsed time: " << duration_ms.count() / (float)args.size() << " milliseconds\n";
    std::cout << "Elapsed time: " << duration_sec.count() / (float)args.size() << " seconds\n";
};
RunSaxs::~RunSaxs() {};