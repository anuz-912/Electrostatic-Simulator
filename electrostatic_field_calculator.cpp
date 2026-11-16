#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <fstream>
#include <string>

using namespace std;

double k(double p, double rm) {
    double b = 1.12; // 1/e0 * 10^-11
    return p * b / rm;
}

int main() {
    int row = 75, col = 75;
    double tp, btm, lf, rg;
    double p, rm;

    vector<vector<double>> T(row, vector<double>(col, 0.0));
    vector<vector<int>> is_fixed(row, vector<int>(col, 0));

    cout << "Boundary conditions (TOP BOTTOM LEFT RIGHT): ";
    cin >> tp >> btm >> lf >> rg;

    cout << "Charge density (10^-11) and relative permittivity: ";
    cin >> p >> rm;

    double mew = k(p, rm);

    // Outer boundaries
    for (int j = 0; j < col; j++) {
        T[0][j] = btm;      is_fixed[0][j] = 1;
        T[row - 1][j] = tp; is_fixed[row - 1][j] = 1;
    }
    for (int i = 1; i < row - 1; i++) {
        T[i][0] = lf;       is_fixed[i][0] = 1;
        T[i][col - 1] = rg; is_fixed[i][col - 1] = 1;
    }

    cout << "\nElectrode geometry:\n"
         << "1. Free Space\n"
         << "2. Parallel Plate Capacitor\n"
         << "3. Circular Conductor\n"
         << "4. Central Rod\n"
         << "5. Dual Electrode\n"
         << "Choice: ";
    int choice;
    cin >> choice;

    switch (choice) {
        case 2: {
            double v_pos, v_neg;
            cout << "Top plate V: "; cin >> v_pos;
            cout << "Bottom plate V: "; cin >> v_neg;
            
            int y_top = row / 2 + 5, y_bot = row / 2 - 5;
            int x_start = col / 2 - 5, x_end = col / 2 + 5;
            for (int j = x_start; j <= x_end; j++) {
                T[y_top][j] = v_pos; is_fixed[y_top][j] = 1;
                T[y_bot][j] = v_neg; is_fixed[y_bot][j] = 1;
            }
            break;
        }
        case 3: {
            double v_center;
            int radius;
            cout << "Conductor V: "; cin >> v_center;
            cout << "Radius (grid units): "; cin >> radius;
            
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
            cout << "Rod V: "; cin >> v_rod;
            int cy = row / 2;
            for (int j = 0; j < col; j++) {
                T[cy][j] = v_rod;
                is_fixed[cy][j] = 1;
            }
            break;
        }
        case 5: {
            double v_left_el, v_right_el;
            cout << "Left electrode V: "; cin >> v_left_el;
            cout << "Right electrode V: "; cin >> v_right_el;

            int cy = row / 2;
            int thickness = max(1, row / 20);
            int length = col / 4;

            for (int i = cy - thickness; i <= cy + thickness; i++) {
                for (int j = 0; j <= length; j++) {
                    T[i][j] = v_left_el;
                    is_fixed[i][j] = 1;
                }
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

    // Gauss-Seidel + SOR
    double tol = 1e-5;
    double omega = 1.75;
    int n = 0;
    double max_change;

    do {
        n++;
        max_change = 0.0;

        for (int i = 1; i < row - 1; i++) {
            for (int j = 1; j < col - 1; j++) {
                if (is_fixed[i][j]) continue;

                double T_old = T[i][j];
                double T_gs = (T[i - 1][j] + T[i + 1][j] + T[i][j - 1] + T[i][j + 1] + mew) / 4.0;
                T[i][j] = (1.0 - omega) * T_old + omega * T_gs;

                double diff = fabs(T[i][j] - T_old);
                if (diff > max_change) max_change = diff;
            }
        }
    } while (max_change > tol);

    cout << "\nConverged in " << n << " iterations.\nPotential matrix:\n";
    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f\t", T[i][j]);
        }
        cout << "\n";
    }

    // E-field calculation (h = 1.0)
    vector<vector<double>> Ex(row, vector<double>(col, 0.0));
    vector<vector<double>> Ey(row, vector<double>(col, 0.0));

    for (int i = 1; i < row - 1; i++) {
        for (int j = 1; j < col - 1; j++) {
            Ex[i][j] = -(T[i][j + 1] - T[i][j - 1]) / 2.0;
            Ey[i][j] = -(T[i + 1][j] - T[i - 1][j]) / 2.0;
        }
    }

    for (int j = 1; j < col - 1; j++) {
        Ex[0][j] = -(T[0][j + 1] - T[0][j - 1]) / 2.0;
        Ey[0][j] = -(T[1][j] - T[0][j]);

        Ex[row - 1][j] = -(T[row - 1][j + 1] - T[row - 1][j - 1]) / 2.0;
        Ey[row - 1][j] = -(T[row - 1][j] - T[row - 2][j]);
    }
    for (int i = 1; i < row - 1; i++) {
        Ex[i][0] = -(T[i][1] - T[i][0]);
        Ey[i][0] = -(T[i + 1][0] - T[i - 1][0]) / 2.0;

        Ex[i][col - 1] = -(T[i][col - 1] - T[i][col - 2]);
        Ey[i][col - 1] = -(T[i + 1][col - 1] - T[i - 1][col - 1]) / 2.0;
    }

    Ex[0][0] = -(T[0][1] - T[0][0]);
    Ey[0][0] = -(T[1][0] - T[0][0]);
    Ex[row - 1][0] = -(T[row - 1][1] - T[row - 1][0]);
    Ey[row - 1][0] = -(T[row - 1][0] - T[row - 2][0]);
    Ex[0][col - 1] = -(T[0][col - 1] - T[0][col - 2]);
    Ey[0][col - 1] = -(T[1][col - 1] - T[0][col - 1]);
    Ex[row - 1][col - 1] = -(T[row - 1][col - 1] - T[row - 1][col - 2]);
    Ey[row - 1][col - 1] = -(T[row - 1][col - 1] - T[row - 2][col - 1]);

    cout << "\nElectric Field Vector:\n";
    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f i+ %.2f j\t", Ex[i][j], Ey[i][j]);
        }
        cout << "\n";
    }

    cout << "\nOutput filename (without .csv): ";
    string fname;
    cin >> fname;

    string csv_file = fname + ".csv";
    ofstream fptr(csv_file);
    fptr << "x,y,V,Ex,Ey,E\n";

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            double em = hypot(Ex[i][j], Ey[i][j]);
            fptr << j << "," << i << "," << T[i][j] << "," << Ex[i][j] << "," << Ey[i][j] << "," << em << "\n";
        }
    }
    fptr.close();
    cout << "Saved " << csv_file << "\n";

    string cmd = "python map-gen.py " + csv_file;
    if (system(cmd.c_str()) == 0) {
        string open_cmd = "start " + fname + ".png";
        system(open_cmd.c_str());
    }

    return 0;
}