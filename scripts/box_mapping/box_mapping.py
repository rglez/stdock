# Created by fajardo at 4/1/24
import numpy as np
from prody import parsePDB, getCoords


# todo: convert the whole script into a function with necessary parameters
# todo: convert into functions anything that can be repeated or reused. Find the meaning of the following in Internet: Stay DRY to avoid being WET

def divide_axes(axis_max, axis_min, max_size):
    # todo: add docstrings to all functions
    axis_dim = axis_max - axis_min
    n_cuts = np.ceil(axis_dim / max_size)
    lenght = np.round(axis_dim / n_cuts, 2)
    return int(n_cuts), lenght


from tabulate import tabulate


# Parsing protein structure
receptor_path = '/scripts/03_complexes_explorations/input_files/p97ND1/receptor.pdb'
receptor = parsePDB(
    receptor_path)  # Take into account the location of the PDB file

# Get atom coordinates of the protein
coords = getCoords(receptor)

# Obtain maximum and minimum values along each axis
x_coords = coords[:, 0]
y_coords = coords[:, 1]
z_coords = coords[:, 2]

x_max, x_min = x_coords.max(), x_coords.min()
y_max, y_min = y_coords.max(), y_coords.min()
z_max, z_min = z_coords.max(), z_coords.min()

# Print maximum and minimum values along each axis
# todo: is this a necessary info to print?
print("Maximum value along X-axis:", x_max)
print("Minimum value along X-axis:", x_min)
print("Maximum value along Y-axis:", y_max)
print("Minimum value along Y-axis:", y_min)
print("Maximum value along Z-axis:", z_max)
print("Minimum value along Z-axis:", z_min)

# Calculate box dimensions
max_size = 50
nx, len_x = divide_axes(x_max, x_min, max_size)
ny, len_y = divide_axes(y_max, y_min, max_size)
nz, len_z = divide_axes(z_max, z_min, max_size)
print(f"The dimensions of each box are: x: {len_x}, y: {len_y} and z: {len_z}")

# Calculate the minimum and maximum ranges for each box
range_min = []
range_max = []
# todo: unnecessary redeclaration of variables xm, ym, zm. Why?
xm = x_min
ym = y_min
zm = z_min
# todo: avoid nested loops as much as possible. I dont really get this ...
while zm <= z_max - len_z:
    while ym <= y_max:
        while xm <= x_max - len_x:
            range_min.append(xm)
            range_max.append(xm + len_x)
            xm = xm + len_x

            range_min.append(ym)
            range_max.append(ym + len_y)

            range_min.append(zm)
            range_max.append(zm + len_z)

        ym = ym + len_y
        xm = x_min

    ym = y_min
    zm = zm + len_z

# todo: avoid the use of external dependencies as much as possible. Why tabulate?
# Print the ranges of each box in a tabular form
box = []
for i in range(len(range_min)):
    # todo: you just redefined a reserved variable name (list). NEVER do that.
    list = [range_min[i], range_max[i]]
    box.append(list)
    if (i + 1) % 3 == 0:
        print("box", (i + 1) // 3)
        print(tabulate(box, headers=["range_min", "range_max"]))
        box = []
