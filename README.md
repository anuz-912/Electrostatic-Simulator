# 2D Electrostatic Field Simulator

A 2D numerical solver written in C++ that computes electric potential distributions and electric field vectors across custom boundary conditions and internal electrode geometries, complete with automated Python visualization plotting.

---

## Technical Requirements & Dependencies

To build and run this project, ensure you have the following installed on your system:

### 1. C++ Compiler
* **GCC / G++** (supporting C++11 or higher) or MSVC / Clang.

### 2. Python Environment
* **Python 3.x**

### 3. Required Python Libraries
* `numpy`
* `pandas`
* `matplotlib`

You can install the required Python packages using `pip`:
```bash
pip install numpy pandas matplotlib
```

---

## How to Run

1. **Compile the C++ source file:**
   ```bash
   g++ -O2 electrostatic_field_calculator.cpp -o electrostatic_simulator
   ```

2. **Run the executable:**
   ```bash
   ./electrostatic_simulator
   ```

3. **Follow the on-screen prompts** to input boundary voltages, charge parameters, geometry selection, and output filename. The program will automatically generate the corresponding `.csv` data file and launch `map-gen.py` to create and display the resulting field plot (`.png`).

---

## Physics & Mathematical Concepts

For a detailed explanation of the theoretical background, physics, and mathematical methods used in this simulator, please visit:

[https://anujpanthi.com.np/electrostatic_simulator](https://anujpanthi.com.np/electrostatic_simulator)
