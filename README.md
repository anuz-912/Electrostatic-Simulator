# 2D Electrostatic Field Simulator

A C++ numerical solver and Python visualizer for 2D electrostatic potential and electric field distributions on a 75x75 grid[cite: 2, 3]. The C++ program solves the grid using Gauss-Seidel with Successive Over-Relaxation (SOR), exports the data to a CSV file, and automatically calls `map-gen.py` to generate a heatmap and vector plot[cite: 2, 3].

## Features

* **Electrode Presets:**[cite: 2]
  1. Free Space[cite: 2]
  2. Parallel Plate Capacitor[cite: 2]
  3. Circular Conductor[cite: 2]
  4. Central Rod[cite: 2]
  5. Dual Electrode[cite: 2]
* **Customizable Inputs:** Configure outer boundary potentials (Top, Bottom, Left, Right), charge density, and relative permittivity[cite: 2].
* **Automated Visualization:** Exports grid coordinates, potential ($V$), electric field components ($E_x, E_y$), and field magnitude ($E$) to `.csv`[cite: 2], then renders a 300 DPI `.png` plot using `map-gen.py`[cite: 2, 3].

## Setup & Installation

1. **Install Python dependencies:**
   ```bash
   pip install numpy pandas matplotlib
   ```

2. **Compile the C++ solver:**
   ```bash
   g++ -O3 electrostatic_field_calculator.cpp -o electrostatic_sim
   ```

## Usage

1. **Run the simulator:**
   ```bash
   ./electrostatic_sim
   ```
   *(On Windows, run `electrostatic_field_calculator.exe`)*

2. **Follow the prompts:**
   * Enter boundary potentials (`TOP BOTTOM LEFT RIGHT`)[cite: 2]
   * Enter charge density ($10^{-11}$) and relative permittivity[cite: 2]
   * Choose an electrode preset and enter its parameters[cite: 2]
   * Enter an output filename (without `.csv`)[cite: 2]

3. **Standalone plotting (optional):**
   You can also run the visualizer directly on any generated CSV file:
   ```bash
   python map-gen.py your_output.csv
   ```

## Output Files

* `<filename>.csv`: Grid data with columns `x, y, V, Ex, Ey, E`[cite: 2].
* `<filename>.png`: Potential contour heatmap overlaid with normalized electric field arrows[cite: 2, 3].

---

For related math and science visit https://anujpanthi.com.np/electrostatic-simulator/
