#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "filamentFields.h"

namespace py = pybind11;

PYBIND11_MODULE(filamentFields, m) {
    // std::cout << "Filament Processing module initialized" << std::endl;
    py::class_<filamentFields>(m, "filamentFields")
        .def(py::init<const std::vector<Eigen::MatrixXd>&>())
        .def(py::init<const std::vector<Eigen::MatrixXd>&, const Eigen::MatrixXd&>())
        .def("return_filament_nodes_list", &filamentFields::return_filament_nodes_list)
        .def("return_filament_edges_list", &filamentFields::return_filament_edges_list)
        .def("return_number_of_labels", &filamentFields::return_number_of_labels)
        .def("return_volume_fraction", &filamentFields::return_volume_fraction)
        .def("return_orientational_order_parameter", &filamentFields::return_orientational_order_parameter)
        .def("return_entanglement", &filamentFields::return_entanglement)
        .def("return_number_of_local_contacts", &filamentFields::return_number_of_local_contacts)
        .def("return_force_sum", &filamentFields::return_force_sum)
        .def("return_all_nodes", &filamentFields::return_all_nodes)
        .def("return_all_edges", &filamentFields::return_all_edges)
        .def("return_node_labels", &filamentFields::return_node_labels)
        .def("return_edge_labels", &filamentFields::return_edge_labels)
        .def("return_total_linking_matrix", &filamentFields::return_total_linking_matrix)
        .def("compute_total_linking_matrix", &filamentFields::compute_total_linking_matrix)
        .def("sample_edges_locally", &filamentFields::sample_edges_locally)
        .def("analyzeLocalVolume", &filamentFields::analyzeLocalVolume)
        .def("analyzeLocalVolumeFromPrecomputed", &filamentFields::analyzeLocalVolumeFromPrecomputed);
        
}