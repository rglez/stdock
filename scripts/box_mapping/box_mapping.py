import numpy as np
from prody import parsePDB, getCoords

def coords():
    receptor_path = input("input your protein path ")
    receptor = parsePDB(receptor_path)
    coords = getCoords(receptor)

    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    z_coords = coords[:, 2]

    x_max, x_min = x_coords.max(), x_coords.min()
    y_max, y_min = y_coords.max(), y_coords.min()
    z_max, z_min = z_coords.max(), z_coords.min()

    max_size = 50

    coords_in_x = [x_min]
    axis_dim_x = np.round(x_max - x_min, 2)
    n_cuts_x = np.ceil(axis_dim_x / max_size)
    lenght_x = np.round(axis_dim_x / n_cuts_x, 2)
    print(f"in the x axis are {n_cuts_x} boxes, each one with a lenght of {lenght_x} Angstrons")
    while x_min < x_max:
        x_min = np.round(x_min + lenght_x, 2)
        coords_in_x.append(x_min)

    coords_in_y = [y_min]
    axis_dim_y = np.round(y_max - y_min, 2)
    n_cuts_y = np.ceil(axis_dim_y / max_size)
    lenght_y = np.round(axis_dim_y / n_cuts_y, 2)
    print(f"in the Y axis are {n_cuts_y} boxes, each one with a lenght of {lenght_y} Angstrons")
    while y_min < y_max:
        y_min = np.round(y_min + lenght_y, 2)
        coords_in_y.append(y_min)

    coords_in_z = [z_min]
    axis_dim_z = np.round(z_max - z_min, 2)
    n_cuts_z = np.ceil(axis_dim_z / max_size)
    lenght_z = np.round(axis_dim_z / n_cuts_z, 2)
    print(f"in the z axis are {n_cuts_z} boxes, each one with a lenght of {lenght_z} Angstrons")
    while z_min < z_max:
        z_min = np.round(z_min + lenght_z, 2)
        coords_in_z.append(z_min)

    return f"coordinates of each vertex in X axis {coords_in_x}, coordinates of each vertex in Y axis {coords_in_y}, coordinates of each vertex in Z axis {coords_in_z}"

print(coords())