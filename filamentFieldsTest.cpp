#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "filamentFields.h"

int main() {
    std::vector<Eigen::MatrixXd> filament_nodes_list;

    // Uncomment and use these nodes for testing
    // Eigen::MatrixXd nodes1(2, 3);
    // nodes1 << -2, 0, 0,              
    //           2, 0, 0;
    // Eigen::MatrixXd nodes2(2, 3);
    // nodes2 << 0, -2, 0.0001,
    //           0, 2, 0.0001;
    // filament_nodes_list.push_back(nodes1);
    // filament_nodes_list.push_back(nodes2);

    int num_rods = 100;
    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    // generate random rods
    for (int i = 0; i < num_rods; i++) {
        Eigen::MatrixXd nodes(10, 3);
        for (int j = 0; j < 10; j++) {
            nodes(j, 0) = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
            nodes(j, 1) = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
            nodes(j, 2) = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        }
        filament_nodes_list.push_back(nodes);
    }

    filamentFields filament(filament_nodes_list);
    Eigen::Vector3d query_point(0.0, 0.0, 0.0);
    double R_omega = 0.1;
    double rod_radius = 0.1;
    filament.precompute(R_omega);
    
    // Uncomment these lines to use the filament analysis
    filament.analyzeLocalVolumeFromPrecomputed(query_point, R_omega, rod_radius);
    std::cout << "Number of labels: " << filament.return_number_of_labels() << std::endl;
    std::cout << "Volume fraction: " << filament.return_volume_fraction() << std::endl;
    std::cout << "Orientational order parameter: " << filament.return_orientational_order_parameter() << std::endl;
    std::cout << "Entanglement: " << filament.return_entanglement() << std::endl;
    
    return 0;
}
