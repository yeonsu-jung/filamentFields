//
//  filamentFields_volume_fraction_only.cpp
//  
//
//  Created by Romy Aran on 6/12/24.
//
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "filamentFields.h"
#include <fstream>
#include <string>
#include <sstream>
#include <map>

std::vector<Eigen::MatrixXd> read_csv_to_eigen_matrices(const std::string& filename) {
    std::vector<Eigen::MatrixXd> eigen_matrices;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return eigen_matrices;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> values;

        // Read shape information
        std::getline(ss, value, ',');
        int rows = std::stoi(value);
        std::getline(ss, value, ',');
        int cols = std::stoi(value);

        // Read matrix data
        while (std::getline(ss, value, ',')) {
            values.push_back(std::stod(value));
        }

        // Ensure correct number of values
        if (values.size() != rows * cols) {
            std::cerr << "Mismatch between expected and actual number of matrix elements" << std::endl;
            continue;
        }

        // Convert to Eigen matrix
        Eigen::MatrixXd matrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix(i, j) = values[i * cols + j];
            }
        }

        eigen_matrices.push_back(matrix);
    }

    return eigen_matrices;
}

// Function to write Eigen matrices to CSV file
void write_eigen_matrices_to_csv(const std::vector<Eigen::MatrixXd>& matrices, const std::string& filename) {
    std::ofstream file(filename);

    for (const auto& matrix : matrices) {
        file << matrix.rows() << "," << matrix.cols();
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << "," << matrix(i, j);
            }
        }
        file << "\n";
    }
}


int main() {
    std::vector<Eigen::MatrixXd> eigen_matrices = read_csv_to_eigen_matrices("/Users/romyaran/Downloads/arrays.csv");

    // Check if the matrices are read correctly
    if (eigen_matrices.empty()) {
        std::cerr << "No matrices to process." << std::endl;
        return 1;
    }

    // Write the Eigen matrices to a new CSV file
    write_eigen_matrices_to_csv(eigen_matrices, "eigen_matrices.csv");

    filamentFields filament(eigen_matrices);

    double R_omega = 100;
    double rod_radius = 1;

    // Define the volume and interval
    int x_max = 2000;
    int y_max = 2000;
    int z_max = 100;
    int interval = 10;

    // Open a CSV file to write the entanglement results
    std::ofstream result_file("volume_fraction_results_z100.csv");
    result_file << "x,y,z,vf\n";

    // Iterate over the volume
    for (int x = 0; x <= x_max; x += interval) {
        for (int y = 0; y <= y_max; y += interval) {
            for (int z = z_max; z <= z_max; z += interval) {
                Eigen::Vector3d query_point(x, y, z);
                filament.analyze_local_volume(query_point, R_omega, rod_radius);
                double vf = filament.return_volume_fraction();

                // Write the result to the CSV file
                result_file << x << "," << y << "," << z << "," << vf << "\n";
            }
        }
    }

    result_file.close();
}

