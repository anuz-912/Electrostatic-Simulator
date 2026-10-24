import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Check if filename was passed via command-line argument, else prompt
if len(sys.argv) > 1:
    file_input = sys.argv[1].strip().strip('"\'')
else:
    file_input = input("Please provide the name of your file: ").strip().strip('"\'')

# Validate path
if not os.path.exists(file_input):
    print(f"Error: File '{file_input}' does not exist. Check the path and try again.")
    sys.exit(1)

# Generate output image name based on the CSV filename
base_name, _ = os.path.splitext(file_input)
output_image = f"{base_name}.png"

# 1. Load CSV data
df = pd.read_csv(file_input)

# 2. Extract columns
x = df['x'].values
y = df['y'].values
V = df['V'].values
Ex = df['Ex'].values
Ey = df['Ey'].values

# 3. Create the plot
plt.figure(figsize=(9, 7))

# Contour / Heatmap for Electric Potential (V)
contour = plt.tricontourf(x, y, V, levels=50, cmap='viridis')
plt.colorbar(contour, label='Electric Potential $V$ (V)')

# 4. Normalized field vector downsampling (Stride = 3 for clear visual resolution)
stride = 3
x_sub = x[::stride]
y_sub = y[::stride]
Ex_sub = Ex[::stride]
Ey_sub = Ey[::stride]

E_mag = np.hypot(Ex_sub, Ey_sub)
E_mag_safe = np.where(E_mag == 0, 1.0, E_mag)

Ex_norm = Ex_sub / E_mag_safe
Ey_norm = Ey_sub / E_mag_safe

# Quiver plot with normalized directional arrows
plt.quiver(x_sub, y_sub, Ex_norm, Ey_norm, color='white', alpha=0.75, pivot='mid', scale=30)

# Force axes to span the complete grid range including boundary nodes
x_min, x_max = int(x.min()), int(x.max())
y_min, y_max = int(y.min()), int(y.max())

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

# Force integer ticks at multiples of 5 for x and y axes
plt.xticks(np.arange(x_min, x_max + 1, 5))
plt.yticks(np.arange(y_min, y_max + 1, 5))

plt.xlabel('X position')
plt.ylabel('Y position')
plt.title('Electric Potential Heatmap & Field Vectors ($E$)')
plt.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()

# Save output image with the same base name as the CSV file
plt.savefig(output_image, dpi=300)
print(f"Plot saved successfully as {output_image}")