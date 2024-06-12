//For finding the covariance between the entanglement field at different values as a function of distance (e.g. what is the typical covariance between a point and a point x units away from it?)

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
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <cmath> // For isnan
#include <unordered_map>

using namespace std;

// Helper function to check if a string is NaN
bool isNotNumber(const string& value) {
    try {
        stod(value);
    } catch (const invalid_argument& e) {
        return true;
    } catch (const out_of_range& e) {
        return true;
    }
    return false;
}

// Function to parse CSV file and ignore rows with NaN values
vector<vector<double>> parseCSV(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        vector<double> row;
        bool hasNaN = false;

        while (getline(ss, value, ',')) {
            if (isNotNumber(value)) {
                hasNaN = true;
                break;
            }
            row.push_back(stod(value));
        }

        if (!hasNaN && row.size() == 4) {
            data.push_back(row);
        }
    }

    return data;
}

// Function to calculate the mean entanglement
double calculateMeanEntanglement(const vector<vector<double>>& data) {
    double sum = 0;
    for (const auto& row : data) {
        sum += row[3];
    }
    return sum / data.size();
}

// Function to calculate A(_)
vector<double> calculateA(const vector<vector<double>>& data, double meanEntanglement) {
    vector<double> A(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        A[i] = data[i][3] - meanEntanglement;
    }
    return A;
}

// Function to calculate <A(*) * A(* + d)> for a given column index
vector<pair<int, double>> calculateCorrelationFunction(const vector<vector<double>>& data, const vector<double>& A, int columnIndex) {
    unordered_map<int, pair<double, int>> correlations;

    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = i + 1; j < data.size(); ++j) {
            int d = abs(static_cast<int>(data[j][columnIndex] - data[i][columnIndex]));
            correlations[d].first += A[i] * A[j];
            correlations[d].second += 1;
        }
    }

    vector<pair<int, double>> result;
    for (const auto& pair : correlations) {
        result.push_back({pair.first, pair.second.first / pair.second.second});
    }
    return result;
}

// Function to write the results to a CSV file
void writeCSV(const string& filename, const vector<pair<int, double>>& data, const string& col1, const string& col2) {
    ofstream file(filename);
    file << col1 << "," << col2 << endl;
    for (const auto& pair : data) {
        file << pair.first << "," << pair.second << endl;
    }
}

int main() {
    string filename = "/Users/romyaran/entanglement_results_z100.csv";
    vector<vector<double>> data = parseCSV(filename);

    if (data.empty()) {
        cout << "No valid data available." << endl;
        return 1;
    }

    // Calculate mean entanglement
    double meanEntanglement = calculateMeanEntanglement(data);

    // Calculate A(_)
    vector<double> A = calculateA(data, meanEntanglement);

    // Calculate <A(x) * A(x + a)>
    vector<pair<int, double>> correlationFunctionX = calculateCorrelationFunction(data, A, 0);
    writeCSV("correlation_x.csv", correlationFunctionX, "a", "<A(x)*A(x+a)>");

    // Calculate <A(y) * A(y + b)>
    vector<pair<int, double>> correlationFunctionY = calculateCorrelationFunction(data, A, 1);
    writeCSV("correlation_y.csv", correlationFunctionY, "b", "<A(y)*A(y+b)>");

    // Calculate <A(z) * A(z + c)>
    vector<pair<int, double>> correlationFunctionZ = calculateCorrelationFunction(data, A, 2);
    writeCSV("correlation_z.csv", correlationFunctionZ, "c", "<A(z)*A(z+c)>");

    cout << "Correlation functions calculated and saved to CSV files." << endl;

    return 0;
}


