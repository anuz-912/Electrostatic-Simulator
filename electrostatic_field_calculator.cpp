#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

double k(double p, double rm) {
    double b = 1.12; // b = 1/e0 * 10^-11
    return p * b / rm; // assumed grid spacing h = 1.00
}

int main() {
    int row = 75, col = 75;
    double tp, btm, lf, rg;
    double p, rm;

    vector<vector<double>> T(row, vector<double>(col, 0.00));
    vector<vector<int>> is_fixed(row, vector<int>(col, 0));

    cout << "Please provide the boundary conditions in order of TOP BOTTOM LEFT RIGHT: ";
    cin >> tp >> btm >> lf >> rg;
    cout << endl;
    cout << "Please provide the charge density(magnitude of order 10^-11) and relative permittivity of the medium: ";
    cin >> p >> rm;
    cout << endl;

    double mew = k(p, rm);

    // Initialize boundary conditions
    for (int j = 0; j < col; j++) {
        T[0][j] = btm;           is_fixed[0][j] = 1;      // Bottom row (y = 0)
        T[row - 1][j] = tp;      is_fixed[row - 1][j] = 1; // Top row (y = row-1)
    }
    for (int i = 1; i < row - 1; i++) {
        T[i][0] = lf;            is_fixed[i][0] = 1;      // Left column (x = 0)
        T[i][col - 1] = rg;      is_fixed[i][col - 1] = 1; // Right column (x = col-1)
    }

    // Stage 8: Internal Geometry Selection
    cout << "Select internal electrode geometry:" << endl;
    cout << "1. Free Space (No internal electrodes)" << endl;
    cout << "2. Parallel Plate Capacitor" << endl;
    cout << "3. Circular Conductor" << endl;
    cout << "4. Central Rod Electrode" << endl;
    cout << "5. Dual Electrode " << endl;

    cout << "Choice: ";
    int choice;
    cin >> choice;
    cout << endl;

    switch (choice) {
        case 2: {
            double v_pos, v_neg;
            cout << "Enter top plate voltage: "; cin >> v_pos;
            cout << "Enter bottom plate voltage: "; cin >> v_neg;
            cout << endl;
            int y_top = (row/2 + 5), y_bot = (row/2 - 5);
            int x_start = (col/2 - 5), x_end = (col/2 + 5);
            for (int j = x_start; j <= x_end; j++) {
                T[y_top][j] = v_pos; is_fixed[y_top][j] = 1;
                T[y_bot][j] = v_neg; is_fixed[y_bot][j] = 1;
            }
            break;
        }
        case 3: {
            double v_center;
            int radius;
            cout << "Enter central conductor voltage: "; cin >> v_center;
            cout << "Enter radius in grid units: "; cin >> radius;
            cout << endl;
            int cx = col / 2, cy = row / 2;
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    if ((i - cy) * (i - cy) + (j - cx) * (j - cx) <= radius * radius) {
                        T[i][j] = v_center;
                        is_fixed[i][j] = 1;
                    }
                }
            }
            break;
        }
        case 4: {
            double v_rod;
            cout << "Enter central rod voltage: "; cin >> v_rod;
            cout << endl;
            int cy = row / 2;
            for (int j = 0; j < col; j++) {
                T[cy][j] = v_rod;
                is_fixed[cy][j] = 1;
            }
            break;
        }
        case 5: {
            double v_left_el, v_right_el;
            cout << "Enter left internal electrode voltage: "; cin >> v_left_el;
            cout << "Enter right internal electrode voltage: "; cin >> v_right_el;
            cout << endl;

            int cy = row / 2;
            int thickness = max(1, row / 20);
            int length = col / 4;

            int y_start = cy - thickness;
            int y_end = cy + thickness;

            for (int i = y_start; i <= y_end; i++) {
                for (int j = 0; j <= length; j++) {
                    T[i][j] = v_left_el;
                    is_fixed[i][j] = 1;
                }
            }

            for (int i = y_start; i <= y_end; i++) {
                for (int j = col - 1 - length; j < col; j++) {
                    T[i][j] = v_right_el;
                    is_fixed[i][j] = 1;
                }
            }
            break;
        }
        default:
            break;
    }

    double tolerance = 1e-5;
    double omega = 1.75; // SOR relaxation parameter for ~75x75 grid
    int n = 0;
    double max_change;

    // In-place Gauss-Seidel Iteration with SOR acceleration
    do {
        n++;
        max_change = 0.0;

        for (int i = 1; i < row - 1; i++) {
            for (int j = 1; j < col - 1; j++) {
                if (is_fixed[i][j]) continue;

                double T_old = T[i][j];
                double T_gs = (T[i - 1][j] + T[i + 1][j] + T[i][j - 1] + T[i][j + 1] + mew) / 4.0;
                T[i][j] = (1.0 - omega) * T_old + omega * T_gs;

                double change = fabs(T[i][j] - T_old);
                if (change > max_change) {
                    max_change = change;
                }
            }
        }
    } while (max_change > tolerance);

    cout << endl << "Potential saturated under " << n << " iterations (Gauss-Seidel + SOR)" << endl;
    cout << "Calculated potential matrix:" << endl;

    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f\t", T[i][j]);
        }
        cout << endl;
    }

    vector<vector<double>> Ex(row, vector<double>(col, 0.00));
    vector<vector<double>> Ey(row, vector<double>(col, 0.00));
    double h = 1.00;

    // 1. INTERIOR NODES (Central Differences)
    for (int i = 1; i < row - 1; i++) {
        for (int j = 1; j < col - 1; j++) {
            Ex[i][j] = -(T[i][j + 1] - T[i][j - 1]) / (2.0 * h);
            Ey[i][j] = -(T[i + 1][j] - T[i - 1][j]) / (2.0 * h);
        }
    }

    // 2. BOUNDARY EDGES
    for (int j = 1; j < col - 1; j++) {
        Ex[0][j] = -(T[0][j + 1] - T[0][j - 1]) / (2.0 * h);
        Ey[0][j] = -(T[1][j] - T[0][j]) / h;
    }
    for (int j = 1; j < col - 1; j++) {
        Ex[row - 1][j] = -(T[row - 1][j + 1] - T[row - 1][j - 1]) / (2.0 * h);
        Ey[row - 1][j] = -(T[row - 1][j] - T[row - 2][j]) / h;
    }
    for (int i = 1; i < row - 1; i++) {
        Ex[i][0] = -(T[i][1] - T[i][0]) / h;
        Ey[i][0] = -(T[i + 1][0] - T[i - 1][0]) / (2.0 * h);
    }
    for (int i = 1; i < row - 1; i++) {
        Ex[i][col - 1] = -(T[i][col - 1] - T[i][col - 2]) / h;
        Ey[i][col - 1] = -(T[i + 1][col - 1] - T[i - 1][col - 1]) / (2.0 * h);
    }

    // 3. CORNER NODES
    Ex[0][0] = -(T[0][1] - T[0][0]) / h;
    Ey[0][0] = -(T[1][0] - T[0][0]) / h;

    Ex[row - 1][0] = -(T[row - 1][1] - T[row - 1][0]) / h;
    Ey[row - 1][0] = -(T[row - 1][0] - T[row - 2][0]) / h;

    Ex[0][col - 1] = -(T[0][col - 1] - T[0][col - 2]) / h;
    Ey[0][col - 1] = -(T[1][col - 1] - T[0][col - 1]) / h;

    Ex[row - 1][col - 1] = -(T[row - 1][col - 1] - T[row - 1][col - 2]) / h;
    Ey[row - 1][col - 1] = -(T[row - 1][col - 1] - T[row - 2][col - 1]) / h;

    cout << endl << "Electric Field Vector:" << endl;
    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f i+ %.2f j\t", Ex[i][j], Ey[i][j]);
        }
        cout << endl;
    }

    cout << endl << "Please provide the name to save the output on .csv extension: ";
    string filename_input;
    cin >> filename_input;
    cout << endl;

    string csv_file = filename_input + ".csv";
    string image_file = filename_input + ".png";

    ofstream fptr(csv_file);
    fptr << "x,y,V,Ex,Ey,E\n";

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            double em = sqrt(Ex[i][j] * Ex[i][j] + Ey[i][j] * Ey[i][j]);
            fptr << j << "," << i << "," << T[i][j] << "," << Ex[i][j] << "," << Ey[i][j] << "," << em << "\n";
        }
    }
    fptr.close();
    cout << "Output saved to file " << csv_file << " successfully." << endl;

    cout << "\nExecuting map-gen.py plotting script..." << endl;
    string python_cmd = "python map-gen.py " + csv_file;
    
    int ret = system(python_cmd.c_str());
    if (ret != 0) {
        string python3_cmd = "python3 map-gen.py " + csv_file;
        ret = system(python3_cmd.c_str());
    }

    if (ret == 0) {
        cout << "Opening output plot image: " << image_file << endl;
        #if defined(_WIN32)
            string open_cmd = "start " + image_file;
        #elif defined(__APPLE__)
            string open_cmd = "open " + image_file;
        #else
            string open_cmd = "xdg-open " + image_file;
        #endif
        system(open_cmd.c_str());
    }

    return 0;
}