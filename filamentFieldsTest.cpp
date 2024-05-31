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

    int num_rods = 1000;
    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    // generate random rods
    for (int i = 0; i < num_rods; i++) {
        Eigen::MatrixXd nodes(2, 3);
        nodes << 2.0 * (std::rand() % 1000) / 1000.0 - 1.0, 2.0 * (std::rand() % 1000) / 1000.0 - 1.0, 0.0001,
                 2.0 * (std::rand() % 1000) / 1000.0 - 1.0, 2.0 * (std::rand() % 1000) / 1000.0 - 1.0, 0.0001;
        filament_nodes_list.push_back(nodes);
    }

    filamentFields filament(filament_nodes_list);
    Eigen::Vector3d query_point(0.5, 0.5, 0.0);
    double R_omega = 10.;
    double rod_radius = 0.1;
    
    // Uncomment these lines to use the filament analysis
    filament.analyzeLocalVolume(query_point, R_omega, rod_radius);
    std::cout << "Number of labels: " << filament.return_number_of_labels() << std::endl;
    std::cout << "Volume fraction: " << filament.return_volume_fraction() << std::endl;
    std::cout << "Orientational order parameter: " << filament.return_orientational_order_parameter() << std::endl;
    std::cout << "Entanglement: " << filament.return_entanglement() << std::endl;
    
    return 0;
}
