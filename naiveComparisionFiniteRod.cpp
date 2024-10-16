#include <bits/stdc++.h>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>

#define double long double
using namespace std;
using namespace Eigen;

int main() {
        ofstream outFile("error.txt");
    // outFile << setprecision(40) << fixed;

    // Parameters
    double L = 1.0;               // Length of the wire (m)
    double a = 0.001;             // Radius of the wire (m)
    double V = 1.0;               // Potential (V)
    double ep_0 = 8.854e-12;      // Permittivity of free space (F/m)

    // Function to compute the Z matrix and solve for charge density
    cout << M_PI << endl;
    auto computeChargeDensity = [&](int N) {
        double dx = L / N;
        MatrixXd Z = MatrixXd::Zero(N, N);
        VectorXd b = 4 * M_PI * ep_0 * VectorXd::Ones(N);

        // Compute the Z matrix
        for (int m = 0; m < N; ++m) {
            for (int n = 0; n < N; ++n) {
                if (m == n) {
                    Z(m, n) = abs(log(abs((dx + sqrt(dx*dx + 4*a*a)) / (dx - sqrt(dx*dx + 4*a*a)))));
                } else {
                    Z(m, n) = abs( log(abs((m + 0.5) * dx - (n + 1)*dx) /abs((m + 0.5) * dx - (n) * dx)));
                    // Z(m, n) = 1 / abs((m-n) * dx); 
                }
            }
        }

        // Solve the system Z * a = b for charge densities a
        MatrixXd Zinv = Z.inverse();
        VectorXd A = Zinv * b;
        return A;
    };

    // Compute charge density for N=160 (reference solution)
    int N_ref = 320;
    VectorXd A_ref = computeChargeDensity(N_ref);

    // List of segment counts to test
    vector<int> segment_counts = {10, 20, 40, 80, 160, 320, 640,1280};

    // File to store errors

    outFile << "Segments Error1 Error2\n";

    // Compute error for different segment counts
    int referr = -1; // Reference error
    outFile << "[" << " " ;
    int cnt = 0;
    for(auto x : A_ref){
        outFile << x << ", "; 
        cnt++;
    }
    outFile << "] \n";
    // for (int j = 0; j < 8; j++) {
    //     int N = segment_counts[j];
    //     VectorXd A_N = computeChargeDensity(N);
    //     double error1 = 0.0;
    //     double error2 = 0.0;
    //     int commonPoints = 10;
    //     outFile << N << " "; 
    //     for (int i = 0; i < commonPoints; ++i) {
    //         error1 = abs(A_N[(N/10) * i]);
    //         outFile  << error1  << " "; 
    //     }
    //     outFile << endl; 
    //     // error1 /= commonPoints;
    //     // if(referr == -1)referr = error1; 

    //     // Write error to file
    //     cout << "Computed error for N = " << N << endl;
    // }

    // outFile.close();
    // cout << "Error data saved to error_data.txt" << endl;

    return 0;
}
