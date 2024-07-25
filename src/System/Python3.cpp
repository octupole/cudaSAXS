#include "Python3.h"
#include "Atoms.h"

void Python3::get_atoms(int frame)
{
    py::scoped_interpreter guard{}; // Start the interpreter and keep it alive for this scope
    std::string result;
    std::map<std::string, std::vector<std::vector<float>>> coord_map;

    try
    {
        py::module_ sys = py::module_::import("sys");
        sys.attr("path").attr("append")(PY_SOURCE_DIR);
        // Set up the arguments to pass to the Python script

        py::list argv;
        argv.append("count_molecules.py");
        sys.attr("argv") = argv;

        // Load the Python script
        py::module_ script = py::module_::import("trajectories");

        // Create a StructureCentering instance (adjust arguments as needed)
        py::object trajStructure = script.attr("TrajectoryStructures");
        py::object analyzer = trajStructure(tpr_file, xtc_file);

        std::vector<int> frames{100, 1000, 1001, 1002};
        for (auto frame : frames)
        {
            std::cout << "Frame: " << frame << std::endl;
            py::object dims = analyzer.attr("get_dimensions_frame")(frame);
            auto box_dimensions = dims.cast<std::vector<float>>();
            Cell::calculateMatrices(box_dimensions);
            auto &co = Cell::getCO();
            auto &oc = Cell::getOC();

            py::object coords_obj = analyzer.attr("get_centered_frame")(frame);
            py::dict coords_dict = coords_obj.cast<py::dict>();
            std::cout << "Coords: " << coords_dict.size() << std::endl;
            for (auto item : coords_dict)
            {
                std::string key = py::str(item.first);

                py::array_t<float> value = item.second.cast<py::array_t<float>>();
                // Check if array is 2D and has correct dimensions
                std::vector<std::vector<float>> vec0 = value.cast<std::vector<std::vector<float>>>();
                coord_map[key] = vec0;
            }
        }
    }

    catch (const py::error_already_set &e)
    {
        std::cerr << "Python error: " << e.what() << std::endl;
    }
}

Python3::~Python3()
{
}