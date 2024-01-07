# Created by fajardo at 4/1/24

from prody import parsePDB, getCoords
from tabulate import tabulate


# Parsing protein structure
receptor = parsePDB("/media/fajardo/Fajardo/from_ROY/Docking/ledock/p97ND1.pdb")  # Take into account the location of the PDB file

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
print("Maximum value along X-axis:", x_max)
print("Minimum value along X-axis:", x_min)
print("Maximum value along Y-axis:", y_max)
print("Minimum value along Y-axis:", y_min)
print("Maximum value along Z-axis:", z_max)
print("Minimum value along Z-axis:", z_min)

# Calculate box dimensions
dx = x_max - x_min
dy = y_max - y_min
dz = z_max - z_min

nx = 1

while dx/nx >= 50.0:
    nx = nx + 1
else:
    len_x = dx/nx

ny = 1

while dy/ny >= 50.0:
    ny = ny + 1
else:
    len_y = dy/ny

nz = 1

while dz/nz >= 50.0:
    nz = nz + 1
else:
    len_z = dz/nz

print("The dimensions of each box are: x = %.3f, y = %.3f, z = %.3f" % (len_x, len_y, len_z))

# Calculate the minimum and maximum ranges for each box
range_min = []
range_max = []

xm = x_min
ym = y_min
zm = z_min

while zm <= z_max-len_z:

    while ym <= y_max:

        while xm <= x_max-len_x:
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

#Print the ranges of each box in a tabular form
box = []
for i in range(len(range_min)):
    list = [range_min[i], range_max[i]]
    box.append(list)
    if (i+1)%3==0:
        print("box", (i + 1) // 3)
        print(tabulate(box, headers=["range_min", "range_max"]))
        box = []
