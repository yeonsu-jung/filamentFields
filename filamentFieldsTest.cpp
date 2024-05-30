#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "filamentFields.h"


int main() {
    std::vector<Eigen::MatrixXd> filament_nodes_list;

    Eigen::MatrixXd nodes1(2, 3);
    nodes1 << -2, 0, 0,              
              2, 0, 0;
    Eigen::MatrixXd nodes2(2, 3);
    nodes2 << 0, -2, 0.0001,
              0, 2, 0.0001;
    filament_nodes_list.push_back(nodes1);
    filament_nodes_list.push_back(nodes2);

    filamentFields filament(filament_nodes_list);
    Eigen::Vector3d query_point(0.5, 0.5, 0.0);
    double R_omega = 10.0;
    double rod_radius = 0.1;
    
    // filament.analyzeLocalVolume(query_point, R_omega, rod_radius);
    // std::cout << "Number of labels: " << filament.return_number_of_labels() << std::endl;
    // std::cout << "Volume fraction: " << filament.return_volume_fraction() << std::endl;
    // std::cout << "Orientational order parameter: " << filament.return_orientational_order_parameter() << std::endl;
    // std::cout << "Entanglement: " << filament.return_entanglement() << std::endl;
    return 0;

}
